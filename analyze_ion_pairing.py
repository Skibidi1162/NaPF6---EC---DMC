#!/usr/bin/env python3
"""Analyze Na-PF6 RDFs, CNs, and SIP/CIP/AGG populations per case.

This command is intentionally independent of electrolyte_generator_v2.py. It
uses the GROMACS installation already required by the project and writes every
result below one generated case directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import statistics
import subprocess
import sys
from typing import Iterable, Sequence


PRIMARY_CUTOFF_NM = 0.5
DEFAULT_SENSITIVITY_CUTOFFS_NM = (0.35, 0.4, 0.405, 0.45, 0.6)
NA_SELECTION = "resname NA and name NA"
P_SELECTION = "resname PF6 and name P"
F_SELECTION = "resname PF6 and name F1 F2 F3 F4 F5 F6"


class AnalysisError(RuntimeError):
    """Raised for an actionable analysis or input error."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate Na-P/Na-F RDFs and CNs and classify Na environments "
            "as SIP/SSIP, CIP, or AGG. Results stay inside the selected case."
        )
    )
    location = parser.add_mutually_exclusive_group(required=True)
    location.add_argument(
        "--case-label",
        help="Generated case label below --cases-root, for example paper_wEC_30_C_1.5M.",
    )
    location.add_argument(
        "--case-dir",
        type=Path,
        help="Path to an existing generated case directory.",
    )
    parser.add_argument(
        "--cases-root",
        type=Path,
        help="Cases directory; defaults to the cases folder beside this script.",
    )
    parser.add_argument(
        "--tpr",
        type=Path,
        help="Override the production TPR discovered inside the case.",
    )
    parser.add_argument(
        "--xtc",
        type=Path,
        help="Override the production XTC discovered inside the case.",
    )
    parser.add_argument(
        "--gmx",
        default="gmx",
        help="GROMACS executable name or path (default: gmx).",
    )
    parser.add_argument(
        "--output-subdir",
        default="analysis_ion_pairing",
        help="Single directory name created inside the case.",
    )
    parser.add_argument(
        "--primary-cutoff-nm",
        type=positive_float,
        default=PRIMARY_CUTOFF_NM,
        help=(
            "Na-P cutoff for the primary SIP/CIP/AGG result "
            f"(default: {PRIMARY_CUTOFF_NM:g} nm)."
        ),
    )
    parser.add_argument(
        "--sensitivity-cutoff-nm",
        type=positive_float,
        action="append",
        help=(
            "Additional Na-P cutoff; repeat as needed. Defaults to "
            "0.35, 0.4, 0.405, 0.45, and 0.6 nm."
        ),
    )
    parser.add_argument(
        "--no-rdf-minimum-cutoff",
        action="store_true",
        help="Do not add the measured first Na-P RDF minimum as a sensitivity cutoff.",
    )
    parser.add_argument(
        "--rdf-bin-nm",
        type=positive_float,
        default=0.002,
        help="RDF bin width in nm (default: 0.002).",
    )
    parser.add_argument(
        "--rdf-max-nm",
        type=positive_float,
        default=1.2,
        help="Maximum RDF distance in nm (default: 1.2).",
    )
    parser.add_argument(
        "--begin-ps",
        type=nonnegative_float,
        help="First trajectory time to analyze, in ps.",
    )
    parser.add_argument(
        "--end-ps",
        type=nonnegative_float,
        help="Last trajectory time to analyze, in ps.",
    )
    parser.add_argument(
        "--dt-ps",
        type=positive_float,
        help="Analyze frames whose times match this interval, in ps.",
    )
    parser.add_argument(
        "--blocks",
        type=positive_int,
        default=5,
        help="Number of contiguous blocks used for standard errors (default: 5).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow known analysis outputs to be overwritten.",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help=(
            "Keep raw GROMACS XVG/index/selection files. By default these "
            "reproducible intermediates are removed after CSV outputs succeed."
        ),
    )
    return parser.parse_args(argv)


def positive_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive finite number")
    return parsed


def nonnegative_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or parsed < 0:
        raise argparse.ArgumentTypeError("value must be a non-negative finite number")
    return parsed


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return parsed


def resolve_case(args: argparse.Namespace) -> tuple[Path, str]:
    script_root = Path(__file__).resolve().parent
    if args.case_label:
        cases_root = (args.cases_root or script_root / "cases").resolve()
        case_dir = (cases_root / args.case_label).resolve()
        if case_dir.parent != cases_root:
            raise AnalysisError("--case-label must name one direct child of --cases-root")
        label = args.case_label
    else:
        case_dir = args.case_dir.expanduser().resolve()
        label = case_dir.name

    if not case_dir.is_dir():
        raise AnalysisError(f"Case directory does not exist: {case_dir}")
    if Path(args.output_subdir).name != args.output_subdir:
        raise AnalysisError("--output-subdir must be a single directory name")
    return case_dir, label


def resolve_input_override(path: Path | None, case_dir: Path) -> Path | None:
    if path is None:
        return None
    if not path.is_absolute():
        path = case_dir / path
    return path.resolve()


def discover_trajectory_pair(
    case_dir: Path,
    label: str,
    tpr_override: Path | None,
    xtc_override: Path | None,
) -> tuple[Path, Path]:
    tpr_override = resolve_input_override(tpr_override, case_dir)
    xtc_override = resolve_input_override(xtc_override, case_dir)
    if (tpr_override is None) != (xtc_override is None):
        raise AnalysisError("Use --tpr and --xtc together")
    if tpr_override is not None and xtc_override is not None:
        require_file(tpr_override, "TPR")
        require_file(xtc_override, "XTC")
        return tpr_override, xtc_override

    expected_stem = f"prod_01_{label}"
    expected_tpr = case_dir / f"{expected_stem}.tpr"
    expected_xtc = case_dir / f"{expected_stem}.xtc"
    if expected_tpr.is_file() and expected_xtc.is_file():
        return expected_tpr.resolve(), expected_xtc.resolve()

    pairs = []
    for tpr in sorted(case_dir.glob("prod_*.tpr")):
        xtc = tpr.with_suffix(".xtc")
        if xtc.is_file():
            pairs.append((tpr.resolve(), xtc.resolve()))
    if len(pairs) == 1:
        return pairs[0]
    if not pairs:
        raise AnalysisError(
            f"No matching production TPR/XTC pair found in {case_dir}. "
            "Run the production stage or provide --tpr and --xtc."
        )
    stems = ", ".join(pair[0].stem for pair in pairs)
    raise AnalysisError(
        f"Multiple production trajectory pairs found ({stems}); use --tpr and --xtc."
    )


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise AnalysisError(f"{label} file does not exist: {path}")


def prepare_output_dir(case_dir: Path, subdir: str, force: bool) -> Path:
    output_dir = (case_dir / subdir).resolve()
    try:
        output_dir.relative_to(case_dir)
    except ValueError as exc:
        raise AnalysisError("Analysis output must remain inside the case directory") from exc
    if output_dir.exists() and any(output_dir.iterdir()) and not force:
        raise AnalysisError(
            f"Analysis output directory is not empty: {output_dir}\n"
            "Use --force to overwrite the named outputs in this directory."
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


class CommandRunner:
    def __init__(self, cwd: Path, log_path: Path):
        self.cwd = cwd
        self.log_path = log_path
        self.environment = os.environ.copy()
        self.environment["GMX_MAXBACKUP"] = "-1"
        self.log_path.write_text("", encoding="utf-8")

    def run(self, command: Sequence[str]) -> str:
        rendered = shlex.join(str(part) for part in command)
        completed = subprocess.run(
            [str(part) for part in command],
            cwd=self.cwd,
            env=self.environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        with self.log_path.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write(f"$ {rendered}\n")
            handle.write(completed.stdout)
            if completed.stdout and not completed.stdout.endswith("\n"):
                handle.write("\n")
            handle.write(f"[exit {completed.returncode}]\n\n")
        if completed.returncode != 0:
            tail = "\n".join(completed.stdout.splitlines()[-20:])
            raise AnalysisError(
                f"Command failed with exit code {completed.returncode}: {rendered}\n{tail}"
            )
        return completed.stdout


def time_arguments(args: argparse.Namespace) -> list[str]:
    result: list[str] = []
    if args.begin_ps is not None:
        result.extend(["-b", format_number(args.begin_ps)])
    if args.end_ps is not None:
        result.extend(["-e", format_number(args.end_ps)])
    if args.dt_ps is not None:
        result.extend(["-dt", format_number(args.dt_ps)])
    return result


def run_rdf(
    runner: CommandRunner,
    gmx: str,
    tpr: Path,
    xtc: Path,
    selection: str,
    rdf_path: Path,
    cn_path: Path,
    args: argparse.Namespace,
) -> None:
    command = [
        gmx,
        "rdf",
        "-f",
        str(xtc),
        "-s",
        str(tpr),
        "-ref",
        NA_SELECTION,
        "-sel",
        selection,
        "-o",
        str(rdf_path),
        "-cn",
        str(cn_path),
        "-bin",
        format_number(args.rdf_bin_nm),
        "-rmax",
        format_number(args.rdf_max_nm),
    ]
    command.extend(time_arguments(args))
    runner.run(command)


def make_static_index(
    runner: CommandRunner,
    gmx: str,
    tpr: Path,
    selection: str,
    output_path: Path,
) -> list[int]:
    runner.run(
        [
            gmx,
            "select",
            "-s",
            str(tpr),
            "-select",
            selection,
            "-on",
            str(output_path),
        ]
    )
    atoms = parse_ndx(output_path)
    if not atoms:
        raise AnalysisError(f"Selection matched no atoms: {selection}")
    return atoms


def parse_ndx(path: Path) -> list[int]:
    atoms: list[int] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("["):
            continue
        atoms.extend(int(value) for value in line.split())
    return atoms


def iter_xvg(path: Path) -> Iterable[list[float]]:
    found = False
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith(("#", "@", "&")):
                continue
            try:
                row = [float(value) for value in line.split()]
            except ValueError as exc:
                raise AnalysisError(
                    f"Could not parse numeric data in {path} at line {line_number}"
                ) from exc
            found = True
            yield row
    if not found:
        raise AnalysisError(f"No numeric data found in {path}")


def parse_xvg(path: Path) -> list[list[float]]:
    return list(iter_xvg(path))


def combine_rdf_cn(rdf_path: Path, cn_path: Path, csv_path: Path) -> list[list[float]]:
    rdf_rows = parse_xvg(rdf_path)
    cn_rows = parse_xvg(cn_path)
    if len(rdf_rows) != len(cn_rows):
        raise AnalysisError(f"RDF and CN lengths differ: {rdf_path} and {cn_path}")
    combined: list[list[float]] = []
    for rdf_row, cn_row in zip(rdf_rows, cn_rows):
        if len(rdf_row) < 2 or len(cn_row) < 2:
            raise AnalysisError("RDF and CN XVG files must each contain at least two columns")
        combined.append([rdf_row[0], rdf_row[1], cn_row[0], cn_row[1]])
    write_csv(
        csv_path,
        ["rdf_r_nm", "g_r", "cn_r_nm", "running_coordination_number"],
        combined,
    )
    return combined


def write_csv(path: Path, headers: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(rows)


def smooth(values: Sequence[float], radius: int = 2) -> list[float]:
    smoothed = []
    for index in range(len(values)):
        start = max(0, index - radius)
        end = min(len(values), index + radius + 1)
        smoothed.append(sum(values[start:end]) / (end - start))
    return smoothed


def find_first_peak_and_minimum(
    combined: Sequence[Sequence[float]],
    peak_start_nm: float = 0.2,
    peak_end_nm: float = 0.5,
    minimum_end_nm: float = 0.75,
) -> dict[str, float | None]:
    radii = [row[0] for row in combined]
    values = smooth([row[1] for row in combined])
    peak_candidates = [
        index
        for index, radius in enumerate(radii)
        if peak_start_nm <= radius <= peak_end_nm
    ]
    if not peak_candidates:
        return {"peak_r_nm": None, "peak_g_r": None, "minimum_r_nm": None, "minimum_g_r": None}
    peak_index = max(peak_candidates, key=lambda index: values[index])
    minimum_index = None
    for index in range(peak_index + 2, len(values) - 1):
        if radii[index] > minimum_end_nm:
            break
        if values[index - 1] > values[index] <= values[index + 1]:
            minimum_index = index
            break
    if minimum_index is None:
        candidates = [
            index
            for index in range(peak_index + 1, len(values))
            if radii[index] <= minimum_end_nm
        ]
        if candidates:
            minimum_index = min(candidates, key=lambda index: values[index])
    return {
        "peak_r_nm": radii[peak_index],
        "peak_g_r": combined[peak_index][1],
        "minimum_r_nm": radii[minimum_index] if minimum_index is not None else None,
        "minimum_g_r": combined[minimum_index][1] if minimum_index is not None else None,
    }


def unique_cutoffs(values: Iterable[float]) -> list[float]:
    result: list[float] = []
    for value in sorted(values):
        if not any(math.isclose(value, existing, abs_tol=1e-9) for existing in result):
            result.append(value)
    return result


def cutoff_token(cutoff: float) -> str:
    return f"{cutoff:.6f}".rstrip("0").rstrip(".").replace(".", "p")


def write_neighbor_selections(
    path: Path,
    na_atoms: Sequence[int],
    cutoffs: Sequence[float],
) -> list[tuple[float, int]]:
    columns: list[tuple[float, int]] = []
    lines: list[str] = []
    for cutoff in cutoffs:
        token = cutoff_token(cutoff)
        for na_atom in na_atoms:
            lines.append(
                f'"cutoff_{token}_na_atom_{na_atom}" '
                f"{P_SELECTION} and within {cutoff:.6f} of atomnr {na_atom};"
            )
            columns.append((cutoff, na_atom))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return columns


def run_neighbor_counts(
    runner: CommandRunner,
    gmx: str,
    tpr: Path,
    xtc: Path,
    selection_path: Path,
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    command = [
        gmx,
        "select",
        "-f",
        str(xtc),
        "-s",
        str(tpr),
        "-sf",
        str(selection_path),
        "-os",
        str(output_path),
        "-xvg",
        "none",
    ]
    command.extend(time_arguments(args))
    runner.run(command)


def classify_neighbor_count(count: int) -> str:
    if count == 0:
        return "SIP/SSIP"
    if count == 1:
        return "CIP"
    return "AGG"


def analyze_neighbor_rows(
    rows: Sequence[Sequence[float]],
    columns: Sequence[tuple[float, int]],
    cutoffs: Sequence[float],
    na_atoms: Sequence[int],
    requested_blocks: int,
) -> tuple[list[list[object]], list[dict[str, object]], list[list[object]]]:
    expected_columns = 1 + len(columns)
    for row in rows:
        if len(row) != expected_columns:
            raise AnalysisError(
                f"Neighbor-count output has {len(row)} columns; expected {expected_columns}"
            )

    offset_by_cutoff: dict[float, int] = {}
    offset = 0
    for cutoff in cutoffs:
        offset_by_cutoff[cutoff] = offset
        offset += len(na_atoms)

    frame_rows: list[list[object]] = []
    summaries: list[dict[str, object]] = []
    histogram_rows: list[list[object]] = []

    for cutoff in cutoffs:
        frame_fractions: dict[str, list[float]] = {
            "SIP/SSIP": [],
            "CIP": [],
            "AGG": [],
        }
        total_by_state = {state: 0 for state in frame_fractions}
        histogram: dict[int, int] = {}
        mean_neighbors_per_frame: list[float] = []
        start = 1 + offset_by_cutoff[cutoff]
        end = start + len(na_atoms)

        for row in rows:
            counts = [integer_count(value) for value in row[start:end]]
            state_counts = {state: 0 for state in frame_fractions}
            for count in counts:
                state = classify_neighbor_count(count)
                state_counts[state] += 1
                total_by_state[state] += 1
                histogram[count] = histogram.get(count, 0) + 1
            mean_neighbors = sum(counts) / len(counts)
            mean_neighbors_per_frame.append(mean_neighbors)
            frame_rows.append(
                [
                    row[0],
                    cutoff,
                    state_counts["SIP/SSIP"],
                    state_counts["CIP"],
                    state_counts["AGG"],
                    state_counts["SIP/SSIP"] / len(counts),
                    state_counts["CIP"] / len(counts),
                    state_counts["AGG"] / len(counts),
                    mean_neighbors,
                ]
            )
            for state in frame_fractions:
                frame_fractions[state].append(state_counts[state] / len(counts))

        total_environments = len(rows) * len(na_atoms)
        block_count = min(requested_blocks, len(rows))
        summary: dict[str, object] = {
            "cutoff_nm": cutoff,
            "frames": len(rows),
            "na_ions": len(na_atoms),
            "na_environments": total_environments,
            "blocks": block_count,
            "mean_pf6_neighbors_per_na": sum(mean_neighbors_per_frame) / len(rows),
        }
        for state, values in frame_fractions.items():
            key = state.lower().replace("/", "_").replace(" ", "_")
            summary[f"{key}_count"] = total_by_state[state]
            summary[f"{key}_fraction"] = total_by_state[state] / total_environments
            summary[f"{key}_block_standard_error"] = block_standard_error(
                values, block_count
            )
        summaries.append(summary)
        for neighbor_count in sorted(histogram):
            count = histogram[neighbor_count]
            histogram_rows.append(
                [cutoff, neighbor_count, count, count / total_environments]
            )
    return frame_rows, summaries, histogram_rows


def analyze_neighbor_file(
    input_path: Path,
    timeseries_path: Path,
    columns: Sequence[tuple[float, int]],
    cutoffs: Sequence[float],
    na_atoms: Sequence[int],
    requested_blocks: int,
) -> tuple[list[dict[str, object]], list[list[object]], int]:
    """Stream neighbor counts so memory does not scale with trajectory length."""
    frame_count = sum(1 for _ in iter_xvg(input_path))
    block_count = min(requested_blocks, frame_count)
    expected_columns = 1 + len(columns)
    states = ("SIP/SSIP", "CIP", "AGG")

    offset_by_cutoff: dict[float, int] = {}
    offset = 0
    for cutoff in cutoffs:
        offset_by_cutoff[cutoff] = offset
        offset += len(na_atoms)

    totals = {
        cutoff: {state: 0 for state in states}
        for cutoff in cutoffs
    }
    histograms: dict[float, dict[int, int]] = {
        cutoff: {} for cutoff in cutoffs
    }
    neighbor_sums = {cutoff: 0.0 for cutoff in cutoffs}
    block_fraction_sums = {
        cutoff: {state: [0.0] * block_count for state in states}
        for cutoff in cutoffs
    }
    block_frames = [0] * block_count

    headers = [
        "time_ps",
        "cutoff_nm",
        "sip_ssip_count",
        "cip_count",
        "agg_count",
        "sip_ssip_fraction",
        "cip_fraction",
        "agg_fraction",
        "mean_pf6_neighbors_per_na",
    ]
    with timeseries_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        for frame_index, row in enumerate(iter_xvg(input_path)):
            if len(row) != expected_columns:
                raise AnalysisError(
                    f"Neighbor-count output has {len(row)} columns; "
                    f"expected {expected_columns}"
                )
            block_index = min(
                frame_index * block_count // frame_count,
                block_count - 1,
            )
            block_frames[block_index] += 1
            for cutoff in cutoffs:
                start = 1 + offset_by_cutoff[cutoff]
                end = start + len(na_atoms)
                counts = [integer_count(value) for value in row[start:end]]
                state_counts = {state: 0 for state in states}
                for count in counts:
                    state = classify_neighbor_count(count)
                    state_counts[state] += 1
                    totals[cutoff][state] += 1
                    histogram = histograms[cutoff]
                    histogram[count] = histogram.get(count, 0) + 1

                mean_neighbors = sum(counts) / len(counts)
                neighbor_sums[cutoff] += mean_neighbors
                fractions = {
                    state: state_counts[state] / len(counts)
                    for state in states
                }
                for state in states:
                    block_fraction_sums[cutoff][state][block_index] += fractions[state]
                writer.writerow(
                    [
                        row[0],
                        cutoff,
                        state_counts["SIP/SSIP"],
                        state_counts["CIP"],
                        state_counts["AGG"],
                        fractions["SIP/SSIP"],
                        fractions["CIP"],
                        fractions["AGG"],
                        mean_neighbors,
                    ]
                )

    total_environments = frame_count * len(na_atoms)
    summaries: list[dict[str, object]] = []
    histogram_rows: list[list[object]] = []
    for cutoff in cutoffs:
        summary: dict[str, object] = {
            "cutoff_nm": cutoff,
            "frames": frame_count,
            "na_ions": len(na_atoms),
            "na_environments": total_environments,
            "blocks": block_count,
            "mean_pf6_neighbors_per_na": neighbor_sums[cutoff] / frame_count,
        }
        for state in states:
            key = state.lower().replace("/", "_").replace(" ", "_")
            block_means = [
                block_fraction_sums[cutoff][state][index] / block_frames[index]
                for index in range(block_count)
                if block_frames[index]
            ]
            error = (
                statistics.stdev(block_means) / math.sqrt(len(block_means))
                if len(block_means) >= 2
                else 0.0
            )
            summary[f"{key}_count"] = totals[cutoff][state]
            summary[f"{key}_fraction"] = (
                totals[cutoff][state] / total_environments
            )
            summary[f"{key}_block_standard_error"] = error
        summaries.append(summary)
        for neighbor_count in sorted(histograms[cutoff]):
            count = histograms[cutoff][neighbor_count]
            histogram_rows.append(
                [cutoff, neighbor_count, count, count / total_environments]
            )
    return summaries, histogram_rows, frame_count


def integer_count(value: float) -> int:
    rounded = int(round(value))
    if not math.isclose(value, rounded, abs_tol=1e-5):
        raise AnalysisError(f"Expected an integer neighbor count; found {value}")
    if rounded < 0:
        raise AnalysisError(f"Neighbor count cannot be negative: {rounded}")
    return rounded


def block_standard_error(values: Sequence[float], block_count: int) -> float:
    if block_count < 2 or len(values) < 2:
        return 0.0
    means = []
    for block in range(block_count):
        start = block * len(values) // block_count
        end = (block + 1) * len(values) // block_count
        if start < end:
            means.append(sum(values[start:end]) / (end - start))
    if len(means) < 2:
        return 0.0
    return statistics.stdev(means) / math.sqrt(len(means))


def interpolate_cn(combined: Sequence[Sequence[float]], radius_nm: float) -> float:
    if radius_nm <= combined[0][2]:
        return combined[0][3]
    for left, right in zip(combined, combined[1:]):
        if left[2] <= radius_nm <= right[2]:
            span = right[2] - left[2]
            if span == 0:
                return left[3]
            fraction = (radius_nm - left[2]) / span
            return left[3] + fraction * (right[3] - left[3])
    return combined[-1][3]


def file_metadata(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "path": str(path),
        "size_bytes": stat.st_size,
        "modified_time_ns": stat.st_mtime_ns,
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def optional_sha256(path: Path) -> str | None:
    return sha256(path) if path.is_file() else None


def analysis_intermediate_paths(output_dir: Path) -> list[Path]:
    return [
        output_dir / name
        for name in (
            "na_atoms.ndx",
            "pf6_phosphorus_atoms.ndx",
            "rdf_Na_P.xvg",
            "coordination_Na_P.xvg",
            "rdf_Na_F.xvg",
            "coordination_Na_F.xvg",
            "na_pf6_neighbor_selections.dat",
            "na_pf6_neighbor_counts.xvg",
        )
    ]


def remove_analysis_intermediates(output_dir: Path) -> dict[str, object]:
    removed: list[str] = []
    reclaimed_bytes = 0
    for path in analysis_intermediate_paths(output_dir):
        if path.is_file():
            reclaimed_bytes += path.stat().st_size
            path.unlink()
            removed.append(path.name)
    return {
        "intermediates_retained": False,
        "removed_files": removed,
        "reclaimed_bytes": reclaimed_bytes,
    }


def write_summary_text(
    path: Path,
    case_label: str,
    tpr: Path,
    xtc: Path,
    extrema: dict[str, float | None],
    summaries: Sequence[dict[str, object]],
    primary_cutoff: float,
) -> None:
    lines = [
        "Na-PF6 RDF, CN, and ion-pairing analysis",
        f"Case label: {case_label}",
        f"TPR: {tpr.name}",
        f"XTC: {xtc.name}",
        "",
        "Na-P RDF features",
        f"First peak (nm): {format_optional(extrema['peak_r_nm'])}",
        f"First minimum (nm): {format_optional(extrema['minimum_r_nm'])}",
        "",
        "Classification: SIP/SSIP = 0 P; CIP = 1 P; AGG = 2 or more P",
    ]
    for summary in summaries:
        marker = " (primary)" if math.isclose(
            float(summary["cutoff_nm"]), primary_cutoff, abs_tol=1e-9
        ) else ""
        lines.extend(
            [
                "",
                f"Cutoff: {float(summary['cutoff_nm']):.6g} nm{marker}",
                f"Frames: {summary['frames']}",
                f"Na ions: {summary['na_ions']}",
                (
                    "Mean PF6 neighbors per Na: "
                    f"{float(summary['mean_pf6_neighbors_per_na']):.6f}"
                ),
                population_line(summary, "sip_ssip", "SIP/SSIP"),
                population_line(summary, "cip", "CIP"),
                population_line(summary, "agg", "AGG"),
            ]
        )
    lines.extend(
        [
            "",
            "Interpretation note: SIP/SSIP here means zero PF6 phosphorus atoms",
            "inside the chosen cutoff; it does not separately identify fully free ions.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def population_line(summary: dict[str, object], key: str, label: str) -> str:
    fraction = float(summary[f"{key}_fraction"])
    error = float(summary[f"{key}_block_standard_error"])
    count = int(summary[f"{key}_count"])
    return f"{label}: {100 * fraction:.3f}% +/- {100 * error:.3f}% block SE ({count})"


def format_optional(value: float | None) -> str:
    return "not detected" if value is None else f"{value:.6g}"


def format_number(value: float) -> str:
    return f"{value:.12g}"


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.begin_ps is not None and args.end_ps is not None:
        if args.end_ps < args.begin_ps:
            raise AnalysisError("--end-ps must be greater than or equal to --begin-ps")

    case_dir, case_label = resolve_case(args)
    tpr, xtc = discover_trajectory_pair(case_dir, case_label, args.tpr, args.xtc)
    output_dir = prepare_output_dir(case_dir, args.output_subdir, args.force)
    runner = CommandRunner(case_dir, output_dir / "commands.log")

    gmx_version = runner.run([args.gmx, "--version"])
    na_atoms = make_static_index(
        runner, args.gmx, tpr, NA_SELECTION, output_dir / "na_atoms.ndx"
    )
    p_atoms = make_static_index(
        runner, args.gmx, tpr, P_SELECTION, output_dir / "pf6_phosphorus_atoms.ndx"
    )
    if len(na_atoms) != len(p_atoms):
        raise AnalysisError(
            f"Expected equal Na and PF6 counts; found {len(na_atoms)} Na and {len(p_atoms)} P"
        )

    na_p_rdf = output_dir / "rdf_Na_P.xvg"
    na_p_cn = output_dir / "coordination_Na_P.xvg"
    na_f_rdf = output_dir / "rdf_Na_F.xvg"
    na_f_cn = output_dir / "coordination_Na_F.xvg"
    run_rdf(runner, args.gmx, tpr, xtc, P_SELECTION, na_p_rdf, na_p_cn, args)
    run_rdf(runner, args.gmx, tpr, xtc, F_SELECTION, na_f_rdf, na_f_cn, args)

    na_p_combined = combine_rdf_cn(
        na_p_rdf, na_p_cn, output_dir / "rdf_coordination_Na_P.csv"
    )
    combine_rdf_cn(na_f_rdf, na_f_cn, output_dir / "rdf_coordination_Na_F.csv")
    extrema = find_first_peak_and_minimum(na_p_combined)

    sensitivity = (
        args.sensitivity_cutoff_nm
        if args.sensitivity_cutoff_nm is not None
        else list(DEFAULT_SENSITIVITY_CUTOFFS_NM)
    )
    requested_cutoffs = [args.primary_cutoff_nm, *sensitivity]
    if not args.no_rdf_minimum_cutoff and extrema["minimum_r_nm"] is not None:
        requested_cutoffs.append(float(extrema["minimum_r_nm"]))
    cutoffs = unique_cutoffs(requested_cutoffs)

    selection_path = output_dir / "na_pf6_neighbor_selections.dat"
    columns = write_neighbor_selections(selection_path, na_atoms, cutoffs)
    neighbor_xvg = output_dir / "na_pf6_neighbor_counts.xvg"
    run_neighbor_counts(
        runner, args.gmx, tpr, xtc, selection_path, neighbor_xvg, args
    )
    summaries, histogram_rows, frame_count = analyze_neighbor_file(
        neighbor_xvg,
        output_dir / "ion_pairing_timeseries.csv",
        columns,
        cutoffs,
        na_atoms,
        args.blocks,
    )
    write_csv(
        output_dir / "ion_pairing_summary.csv",
        [
            "cutoff_nm",
            "frames",
            "na_ions",
            "na_environments",
            "blocks",
            "mean_pf6_neighbors_per_na",
            "sip_ssip_count",
            "sip_ssip_fraction",
            "sip_ssip_block_standard_error",
            "cip_count",
            "cip_fraction",
            "cip_block_standard_error",
            "agg_count",
            "agg_fraction",
            "agg_block_standard_error",
            "na_p_running_cn_at_cutoff",
        ],
        [
            [
                summary["cutoff_nm"],
                summary["frames"],
                summary["na_ions"],
                summary["na_environments"],
                summary["blocks"],
                summary["mean_pf6_neighbors_per_na"],
                summary["sip_ssip_count"],
                summary["sip_ssip_fraction"],
                summary["sip_ssip_block_standard_error"],
                summary["cip_count"],
                summary["cip_fraction"],
                summary["cip_block_standard_error"],
                summary["agg_count"],
                summary["agg_fraction"],
                summary["agg_block_standard_error"],
                interpolate_cn(na_p_combined, float(summary["cutoff_nm"])),
            ]
            for summary in summaries
        ],
    )
    write_csv(
        output_dir / "pf6_neighbor_histogram.csv",
        ["cutoff_nm", "pf6_neighbors", "na_environment_count", "fraction"],
        histogram_rows,
    )

    write_summary_text(
        output_dir / "analysis_summary.txt",
        case_label,
        tpr,
        xtc,
        extrema,
        summaries,
        args.primary_cutoff_nm,
    )
    storage = (
        {
            "intermediates_retained": True,
            "removed_files": [],
            "reclaimed_bytes": 0,
        }
        if args.keep_intermediates
        else remove_analysis_intermediates(output_dir)
    )
    manifest = {
        "command": "analyze_ion_pairing.py",
        "case_label": case_label,
        "case_directory": str(case_dir),
        "output_directory": str(output_dir),
        "inputs": {"tpr": file_metadata(tpr), "xtc": file_metadata(xtc)},
        "atom_selections": {
            "na": NA_SELECTION,
            "pf6_phosphorus": P_SELECTION,
            "pf6_fluorine": F_SELECTION,
        },
        "atom_counts": {"na": len(na_atoms), "pf6_phosphorus": len(p_atoms)},
        "classification": {
            "SIP/SSIP": "zero PF6 phosphorus atoms inside cutoff",
            "CIP": "one PF6 phosphorus atom inside cutoff",
            "AGG": "two or more PF6 phosphorus atoms inside cutoff",
            "primary_cutoff_nm": args.primary_cutoff_nm,
            "analyzed_cutoffs_nm": cutoffs,
        },
        "rdf": {
            "bin_width_nm": args.rdf_bin_nm,
            "maximum_nm": args.rdf_max_nm,
            "na_p_first_peak_and_minimum": extrema,
        },
        "trajectory_window": {
            "begin_ps": args.begin_ps,
            "end_ps": args.end_ps,
            "dt_ps": args.dt_ps,
            "frames_analyzed": frame_count,
        },
        "storage": storage,
        "block_count_requested": args.blocks,
        "gromacs_version_output": gmx_version.strip(),
        "method_references": [
            "https://doi.org/10.26565/2220-637X-2024-43-03",
            "https://doi.org/10.1039/D3EE00864A",
            "https://doi.org/10.1038/s41467-025-63902-4",
            "https://doi.org/10.1080/08927022.2025.2586506",
            "https://doi.org/10.1039/D0CP03639K",
        ],
        "method_note_sha256": optional_sha256(
            Path(__file__).resolve().parent / "ION_PAIRING_METHODS.md"
        ),
        "outputs": sorted(
            {path.name for path in output_dir.iterdir()} | {"analysis_manifest.json"}
        ),
    }
    (output_dir / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n"
    )

    print(f"Analysis complete: {output_dir}")
    primary = next(
        summary
        for summary in summaries
        if math.isclose(
            float(summary["cutoff_nm"]), args.primary_cutoff_nm, abs_tol=1e-9
        )
    )
    print(population_line(primary, "sip_ssip", "SIP/SSIP"))
    print(population_line(primary, "cip", "CIP"))
    print(population_line(primary, "agg", "AGG"))
    if not args.keep_intermediates:
        print(
            "Removed reproducible analysis intermediates: "
            f"{int(storage['reclaimed_bytes']) / (1024 * 1024):.2f} MiB reclaimed"
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AnalysisError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(2)
