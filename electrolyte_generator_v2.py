"""Generate validated NaPF6/EC/DMC Packmol and GROMACS workflow inputs.

Keep this file with the validated PDB, ITP and MDP files, or pass
--template-dir. The original generator is intentionally left unchanged.

Examples:
  python3 electrolyte_generator_v2.py --paper-case
  python3 electrolyte_generator_v2.py --w-ec 30 --concentration 1.5 --density 1.329
  python3 electrolyte_generator_v2.py --n-ec 76 --n-dmc 173 --salt-pairs 30 --density 1.329
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import asdict, dataclass
from pathlib import Path

AVOGADRO = 6.02214076e23
M_EC = 88.062
M_DMC = 90.078
M_NAPF6 = 167.953

TEMPLATES = (
    "custom_atomtypes.itp",
    "na.itp",
    "pf6.itp",
    "ec.itp",
    "dmc.itp",
    "na.pdb",
    "pf6.pdb",
    "ec.pdb",
    "dmc.pdb",
    "em.mdp",
    "prep.mdp",
    "npt_eq.mdp",
    "prod_01.mdp",
)


@dataclass(frozen=True)
class Case:
    mode: str
    requested_w_ec: float | None
    requested_concentration: float | None
    density: float
    n_ec: int
    n_dmc: int
    n_salt: int
    realized_w_ec: float
    realized_concentration: float
    box_angstrom: float
    box_nm: float
    volume_nm3: float
    expected_atoms: int


def from_counts(
    n_ec: int,
    n_dmc: int,
    n_salt: int,
    density: float,
    mode: str,
    requested_w_ec: float | None = None,
    requested_concentration: float | None = None,
) -> Case:
    if density <= 0:
        raise ValueError("Density must be greater than zero.")
    if min(n_ec, n_dmc, n_salt) < 0:
        raise ValueError("Molecule counts must be non-negative.")
    if n_ec + n_dmc == 0:
        raise ValueError("At least one solvent molecule is required.")

    solvent_mass_units = n_ec * M_EC + n_dmc * M_DMC
    total_mass_units = solvent_mass_units + n_salt * M_NAPF6
    mass_g = total_mass_units / AVOGADRO
    volume_cm3 = mass_g / density
    volume_nm3 = volume_cm3 * 1.0e21
    box_angstrom = (volume_cm3 * 1.0e24) ** (1.0 / 3.0)
    realized_w_ec = 100.0 * n_ec * M_EC / solvent_mass_units
    realized_concentration = (
        n_salt / AVOGADRO / (volume_cm3 * 1.0e-3) if n_salt else 0.0
    )

    return Case(
        mode=mode,
        requested_w_ec=requested_w_ec,
        requested_concentration=requested_concentration,
        density=density,
        n_ec=n_ec,
        n_dmc=n_dmc,
        n_salt=n_salt,
        realized_w_ec=realized_w_ec,
        realized_concentration=realized_concentration,
        box_angstrom=box_angstrom,
        box_nm=box_angstrom / 10.0,
        volume_nm3=volume_nm3,
        expected_atoms=8 * n_salt + 10 * n_ec + 12 * n_dmc,
    )


def solve_target(
    w_ec: float,
    concentration: float,
    density: float,
    total_solvents: int,
) -> Case:
    """Include salt mass when converting density to initial box volume."""
    if not 0 <= w_ec <= 100:
        raise ValueError("EC weight percent must be between 0 and 100.")
    if concentration < 0:
        raise ValueError("Concentration must be non-negative.")
    if total_solvents <= 0:
        raise ValueError("Total solvent count must be positive.")

    fraction = w_ec / 100.0
    n_ec = round(
        total_solvents * fraction * M_DMC
        / (M_EC * (1.0 - fraction) + M_DMC * fraction)
    )
    n_dmc = total_solvents - n_ec
    solvent_mass_units = n_ec * M_EC + n_dmc * M_DMC

    denominator = 1000.0 * density - concentration * M_NAPF6
    if denominator <= 0:
        raise ValueError("Concentration and density are physically incompatible.")
    n_salt = max(0, round(concentration * solvent_mass_units / denominator))

    return from_counts(
        n_ec,
        n_dmc,
        n_salt,
        density,
        "target",
        w_ec,
        concentration,
    )


def itp_atoms(text: str) -> list[tuple[str, float]]:
    result: list[tuple[str, float]] = []
    inside = False
    for raw in text.splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            inside = line.lower().replace(" ", "") == "[atoms]"
            continue
        if inside:
            fields = line.split()
            if fields[0].isdigit() and len(fields) >= 8:
                result.append((fields[4], float(fields[6])))
    return result


def pdb_names(text: str) -> list[str]:
    return [
        line[12:16].strip()
        for line in text.splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]


def load_templates(directory: Path) -> dict[str, str]:
    missing = [name for name in TEMPLATES if not (directory / name).is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing templates: " + ", ".join(missing)
            + f"\nTemplate directory: {directory}"
        )

    files = {
        name: (directory / name).read_text(encoding="utf-8")
        .replace("\r\n", "\n")
        .replace("\r", "\n")
        for name in TEMPLATES
    }

    # Fix residue metadata and match the PF6 topology equilibrium bond length.
    files["dmc.pdb"] = files["dmc.pdb"].replace(" UNK ", " DMC ")
    files["pf6.pdb"] = files["pf6.pdb"].replace("1.580", "1.606")

    # Strip box/model records from one-molecule Packmol templates.
    for name in ("na.pdb", "pf6.pdb", "ec.pdb", "dmc.pdb"):
        lines = [
            line
            for line in files[name].splitlines()
            if line.startswith(("ATOM", "HETATM"))
        ]
        files[name] = "\n".join(lines) + "\nEND\n"

    expected_charges = {"na": 1.0, "pf6": -1.0, "ec": 0.0, "dmc": 0.0}
    for stem, expected_charge in expected_charges.items():
        atoms = itp_atoms(files[f"{stem}.itp"])
        if not atoms:
            raise RuntimeError(f"No atoms found in {stem}.itp.")
        top_names = [name for name, _ in atoms]
        coord_names = pdb_names(files[f"{stem}.pdb"])
        if top_names != coord_names:
            raise RuntimeError(
                f"{stem} atom-order mismatch:\nITP {top_names}\nPDB {coord_names}"
            )
        charge = sum(value for _, value in atoms)
        if not math.isclose(charge, expected_charge, abs_tol=1.0e-5):
            raise RuntimeError(
                f"{stem} charge is {charge:.8f}; expected {expected_charge:.1f}."
            )
    return files


def safe_label(text: str) -> str:
    cleaned = re.sub(r"[^0-9A-Za-z_.-]", "_", text)
    if not cleaned:
        raise ValueError("Case label contains no usable characters.")
    return cleaned


def default_label(case: Case) -> str:
    if case.mode == "paper":
        return "paper_wEC_30_C_1.5M"
    if case.mode == "target":
        return (
            f"wEC_{case.requested_w_ec:g}_"
            f"C_{case.requested_concentration:g}M"
        )
    return f"exact_EC_{case.n_ec}_DMC_{case.n_dmc}_salt_{case.n_salt}"


def case_filenames(label: str) -> dict[str, str]:
    """Return names for files whose contents or results depend on the case."""
    return {
        "packmol_input": f"packmol_{label}.inp",
        "topology": f"system_{label}.top",
        "workflow": f"run_workflow_{label}.sh",
        "summary": f"generation_summary_{label}.txt",
        "readme": f"README_{label}.md",
        "manifest": f"manifest_{label}.json",
        "initial_pdb": f"initial_box_{label}.pdb",
        "initial_gro": f"initial_box_{label}.gro",
        "packmol_log": f"packmol_{label}.log",
        "em": f"em_{label}",
        "prep": f"prep_{label}",
        "npt": f"npt_eq_{label}",
        "prod": f"prod_01_{label}",
        "mdout_em": f"mdout_em_{label}.mdp",
        "mdout_prep": f"mdout_prep_{label}.mdp",
        "mdout_npt": f"mdout_npt_{label}.mdp",
        "mdout_prod": f"mdout_prod_{label}.mdp",
    }


def packmol_input(case: Case, label: str) -> str:
    names = case_filenames(label)
    structures = (
        ("na.pdb", case.n_salt),
        ("pf6.pdb", case.n_salt),
        ("ec.pdb", case.n_ec),
        ("dmc.pdb", case.n_dmc),
    )
    blocks = []
    for filename, count in structures:
        if count == 0:
            continue
        blocks.append(
            f"structure {filename}\n"
            f"    number {count}\n"
            "    inside box 0.0 0.0 0.0 "
            f"{case.box_angstrom:.5f} {case.box_angstrom:.5f} "
            f"{case.box_angstrom:.5f}\n"
            "end structure"
        )
    return (
        "# Generated by electrolyte_generator_v2.py\n"
        "# Packmol/PDB distance unit: angstrom\n"
        "tolerance 2.0\n"
        "filetype pdb\n"
        f"output {names['initial_pdb']}\n\n"
        + "\n\n".join(blocks) + "\n"
    )


def system_top(case: Case) -> str:
    return f"""; Generated by electrolyte_generator_v2.py
#include "oplsaa.ff/forcefield.itp"
#include "custom_atomtypes.itp"
#include "na.itp"
#include "pf6.itp"
#include "ec.itp"
#include "dmc.itp"

[ system ]
NaPF6 electrolyte in EC/DMC

[ molecules ]
; Compound   number
NA_M         {case.n_salt}
PF6_M        {case.n_salt}
EC_M         {case.n_ec}
DMC_M        {case.n_dmc}
"""


def workflow_script(case: Case, label: str) -> str:
    names = case_filenames(label)
    return f"""#!/usr/bin/env bash
set -euo pipefail

# Always run inside this case directory, even when the script is invoked from
# the project root or another working directory. This keeps every generated
# Packmol, GROMACS, log, trajectory, and analysis file with its own case.
CASE_DIR="$(cd -- "$(dirname -- "${{BASH_SOURCE[0]}}")" && pwd)"
cd "$CASE_DIR"

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 packmol|em|prep|npt|prod|analysis"
    exit 2
fi

stage="$1"
CASE_LABEL="{label}"
EXPECTED_ATOMS={case.expected_atoms}
BOX_NM={case.box_nm:.8f}
PACKMOL_INPUT="{names['packmol_input']}"
TOP="{names['topology']}"
INITIAL_PDB="{names['initial_pdb']}"
INITIAL_GRO="{names['initial_gro']}"
PACKMOL_LOG="{names['packmol_log']}"
EM="{names['em']}"
PREP="{names['prep']}"
NPT="{names['npt']}"
PROD="{names['prod']}"
MDOUT_EM="{names['mdout_em']}"
MDOUT_PREP="{names['mdout_prep']}"
MDOUT_NPT="{names['mdout_npt']}"
MDOUT_PROD="{names['mdout_prod']}"

case "$stage" in
packmol)
    command -v packmol >/dev/null
    command -v gmx >/dev/null
    packmol < "$PACKMOL_INPUT" | tee "$PACKMOL_LOG"
    actual_atoms=$(grep -cE '^(ATOM|HETATM)' "$INITIAL_PDB" || true)
    [ "$actual_atoms" = "$EXPECTED_ATOMS" ] || {{
        echo "Expected $EXPECTED_ATOMS atoms; found $actual_atoms." >&2
        exit 1
    }}
    gmx editconf -f "$INITIAL_PDB" -o "$INITIAL_GRO" \\
        -box "$BOX_NM" "$BOX_NM" "$BOX_NM"
    ;;
em)
    [ -f "$INITIAL_GRO" ]
    gmx grompp -f em.mdp -c "$INITIAL_GRO" -p "$TOP" \\
        -o "$EM.tpr" -po "$MDOUT_EM"
    nice -n 10 gmx mdrun -deffnm "$EM" -v -ntmpi 1 -ntomp 8
    ;;
prep)
    [ -f "$EM.gro" ]
    gmx grompp -f prep.mdp -c "$EM.gro" -p "$TOP" \\
        -o "$PREP.tpr" -po "$MDOUT_PREP" -maxwarn 1
    nice -n 10 gmx mdrun -deffnm "$PREP" -v -ntmpi 1 -ntomp 8
    ;;
npt)
    [ -f "$PREP.gro" ] && [ -f "$PREP.cpt" ]
    gmx grompp -f npt_eq.mdp -c "$PREP.gro" -t "$PREP.cpt" \\
        -p "$TOP" -o "$NPT.tpr" -po "$MDOUT_NPT" -maxwarn 2
    nice -n 10 gmx mdrun -deffnm "$NPT" -v -ntmpi 1 -ntomp 8
    ;;
prod)
    [ -f "$NPT.gro" ] && [ -f "$NPT.cpt" ]
    gmx grompp -f prod_01.mdp -c "$NPT.gro" -t "$NPT.cpt" \\
        -p "$TOP" -o "$PROD.tpr" -po "$MDOUT_PROD" -maxwarn 2
    nice -n 10 gmx mdrun -deffnm "$PROD" -v -ntmpi 1 -ntomp 8
    ;;
analysis)
    [ -f "$PROD.xtc" ] && [ -f "$PROD.tpr" ]
    gmx rdf -f "$PROD.xtc" -s "$PROD.tpr" \\
        -ref 'resname NA and name NA' \\
        -sel 'resname EC and name O00' \\
        -o "rdf_Na_ECcarbonylO_${{CASE_LABEL}}.xvg" \\
        -cn "coordination_Na_ECcarbonylO_${{CASE_LABEL}}.xvg"
    gmx rdf -f "$PROD.xtc" -s "$PROD.tpr" \\
        -ref 'resname NA and name NA' \\
        -sel 'resname DMC and name O03' \\
        -o "rdf_Na_DMCcarbonylO_${{CASE_LABEL}}.xvg" \\
        -cn "coordination_Na_DMCcarbonylO_${{CASE_LABEL}}.xvg"
    ;;
*)
    echo "Unknown stage: $stage"
    exit 2
    ;;
esac
"""


def summary(case: Case) -> str:
    requested_w = (
        f"{case.requested_w_ec:g}" if case.requested_w_ec is not None
        else "not specified"
    )
    requested_c = (
        f"{case.requested_concentration:g}"
        if case.requested_concentration is not None else "not specified"
    )
    return (
        f"Mode: {case.mode}\n"
        f"Requested EC wt%: {requested_w}\n"
        f"Requested concentration (mol/L): {requested_c}\n"
        f"Input density (g/cm^3): {case.density:g}\n"
        f"EC molecules: {case.n_ec}\n"
        f"DMC molecules: {case.n_dmc}\n"
        f"NaPF6 pairs: {case.n_salt}\n"
        f"Realized EC wt%: {case.realized_w_ec:.8f}\n"
        f"Realized initial concentration (mol/L): "
        f"{case.realized_concentration:.8f}\n"
        f"Box length (angstrom): {case.box_angstrom:.8f}\n"
        f"Box length (nm): {case.box_nm:.8f}\n"
        f"Initial volume (nm^3): {case.volume_nm3:.8f}\n"
        f"Expected atoms: {case.expected_atoms}\n"
    )


def readme(case: Case, label: str) -> str:
    names = case_filenames(label)
    return f"""# Generated NaPF6/EC/DMC case

Case label: `{label}`

## Reproducibility summary

```text
{summary(case).rstrip()}
```

Density is used only for the initial box. Recalculate simulated density and
concentration after NPT equilibration.

## Run one stage at a time in WSL

The workflow automatically runs inside this case directory, so its outputs
cannot accidentally be written into the project root or another case.

```bash
./{names['workflow']} packmol
./{names['workflow']} em
./{names['workflow']} prep
./{names['workflow']} npt
./{names['workflow']} prod
./{names['workflow']} analysis
```

Inspect every stage before continuing. Packmol/PDB uses angstrom; GROMACS uses
nm. The 0.5 ns production stage is a workflow demonstration and preliminary
structural test, not converged transport-property sampling. Case-dependent
results contain `{label}` in their filenames; reusable PDB, ITP and MDP
templates retain their standard names.
"""


def write_case(
    case: Case,
    files: dict[str, str],
    output_dir: Path,
    force: bool,
    label: str,
) -> Path:
    output_dir = output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()) and not force:
        raise FileExistsError(
            f"Output directory is not empty: {output_dir}\n"
            "Use another directory or pass --force."
        )
    output_dir.mkdir(parents=True, exist_ok=True)

    names = case_filenames(label)
    generated = dict(files)
    generated.update({
        names["packmol_input"]: packmol_input(case, label),
        names["topology"]: system_top(case),
        names["workflow"]: workflow_script(case, label),
        names["summary"]: summary(case),
        names["readme"]: readme(case, label),
    })
    for name, content in generated.items():
        with (output_dir / name).open(
            "w", encoding="utf-8", newline="\n"
        ) as handle:
            handle.write(content)

    manifest = {
        "generator": "electrolyte_generator_v2.py",
        "case_label": label,
        "case": asdict(case),
        "molecular_weights_g_mol": {
            "EC": M_EC,
            "DMC": M_DMC,
            "NaPF6": M_NAPF6,
        },
        "generated_files": sorted(generated),
    }
    with (output_dir / names["manifest"]).open(
        "w", encoding="utf-8", newline="\n"
    ) as handle:
        handle.write(json.dumps(manifest, indent=2) + "\n")
    try:
        (output_dir / names["workflow"]).chmod(0o755)
    except OSError:
        pass
    return output_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--w-ec", type=float)
    parser.add_argument("--concentration", type=float)
    parser.add_argument("--density", type=float)
    parser.add_argument("--total-solvents", type=int, default=250)
    parser.add_argument("--n-ec", type=int)
    parser.add_argument("--n-dmc", type=int)
    parser.add_argument("--salt-pairs", type=int)
    parser.add_argument("--paper-case", action="store_true")
    parser.add_argument("--case-label")
    parser.add_argument("--template-dir", type=Path)
    parser.add_argument("-o", "--output-dir", type=Path)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def resolve_case(args: argparse.Namespace) -> Case:
    exact = (args.n_ec, args.n_dmc, args.salt_pairs)
    target = (args.w_ec, args.concentration)

    if args.paper_case:
        if any(value is not None for value in exact + target) or args.density:
            raise ValueError("--paper-case cannot be combined with other inputs.")
        return from_counts(76, 173, 30, 1.329, "paper", 30.0, 1.5)

    exact_requested = any(value is not None for value in exact)
    target_requested = any(value is not None for value in target)
    if exact_requested and target_requested:
        raise ValueError("Do not mix exact counts with target composition.")

    if exact_requested:
        if not all(value is not None for value in exact) or args.density is None:
            raise ValueError(
                "Exact mode needs --n-ec, --n-dmc, --salt-pairs and --density."
            )
        return from_counts(
            args.n_ec, args.n_dmc, args.salt_pairs, args.density, "exact"
        )

    if target_requested or args.density is not None:
        if args.w_ec is None or args.concentration is None or args.density is None:
            raise ValueError(
                "Target mode needs --w-ec, --concentration and --density."
            )
        return solve_target(
            args.w_ec, args.concentration, args.density, args.total_solvents
        )

    print("No composition arguments supplied; using the paper test case.")
    return from_counts(76, 173, 30, 1.329, "paper", 30.0, 1.5)


def main() -> None:
    args = parse_args()
    try:
        case = resolve_case(args)
        project_dir = Path(__file__).resolve().parent
        source = (args.template_dir or project_dir).resolve()
        files = load_templates(source)
        label = safe_label(args.case_label) if args.case_label else default_label(case)
        # Keep the default independent of the shell's current directory. An
        # explicitly supplied output directory remains under the user's control.
        destination = args.output_dir or project_dir / "cases" / label
        written = write_case(case, files, destination, args.force, label)
    except (ValueError, RuntimeError, FileNotFoundError, FileExistsError) as exc:
        raise SystemExit(f"Error: {exc}") from exc

    print(summary(case), end="")
    print(f"Generated: {written}")
    print(
        "Next: cd into that directory and run "
        f"./{case_filenames(label)['workflow']} packmol"
    )


if __name__ == "__main__":
    main()
