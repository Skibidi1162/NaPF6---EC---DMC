# NaPF6 in EC/DMC: Packmol and GROMACS workflow

This repository documents a test workflow for constructing and simulating a
NaPF6 electrolyte in an ethylene carbonate/dimethyl carbonate (EC/DMC) solvent.
The first test case uses 30 wt% EC, a nominal NaPF6 concentration of 1.5 mol/L,
and an initial density of 1.329 g/cm3.

The workflow is:

```text
composition inputs
      -> electrolyte_generator_v2.py
      -> Packmol initial structure
      -> GROMACS energy minimization
      -> 10 ps preparation
      -> 1 ns NPT equilibration
      -> 0.5 ns test production
      -> density and RDF/coordination-number analysis
```

> **Scientific scope.** This is a GROMACS translation and a short workflow
> demonstration based on the system composition in the Karazin paper. The
> paper used MDNAES and 100 ns of total production sampling (50 successive
> 2 ns simulations per system). Therefore, a 0.5 ns test can demonstrate that
> the construction, simulation, and structural-analysis pipeline works, but it
> is not a complete reproduction of the paper and is not sufficient for
> converged transport properties.

## 1. Software and working environment

The tested environment is Windows with Ubuntu running under WSL and VS Code
connected to WSL. All commands below must be entered in an **Ubuntu/WSL Bash
terminal**, not in Windows PowerShell.

Required programs:

- Python 3
- Git
- Packmol
- GROMACS
- VS Code with the WSL extension (recommended, but not required)
- Grace/Xmgrace (optional, for plotting `.xvg` files)

Check the command-line programs before continuing:

```bash
python3 --version
git --version
packmol --version
gmx --version
```

If `packmol --version` is not supported by the installed Packmol wrapper, use:

```bash
command -v packmol
```

It must return the location of a runnable `packmol` command. Installing Packmol
through Julia is acceptable, but the command still needs to be available in
the WSL `PATH` because the generated `run_workflow_<case-label>.sh` script calls
`packmol` directly.

To open the current WSL directory in VS Code:

```bash
code .
```

The lower-left corner of VS Code should show that the window is connected to
WSL/Ubuntu. A new VS Code terminal will then also be an Ubuntu terminal.

## 2. Download the project

Downloading only `electrolyte_generator_v2.py` is not enough. The generator
also requires validated molecular-coordinate, force-field, and simulation
parameter templates. Download or clone the **complete repository**.

The preferred method is to clone it into the Linux filesystem under the WSL
home directory:

```bash
mkdir -p ~/projects
cd ~/projects
git clone <repository-url> napf6-ec-dmc
cd napf6-ec-dmc
```

Replace `<repository-url>` with the URL shown by the green **Code** button on
the GitHub repository. If the project has already been cloned, use:

```bash
cd ~/projects/napf6-ec-dmc
git pull
```

Alternatively, GitHub's **Code -> Download ZIP** option can be used. Extract
the entire archive, not an individual Python file.

## 3. Required input templates

Keep these files in the same directory as `electrolyte_generator_v2.py`:

```text
electrolyte_generator_v2.py
custom_atomtypes.itp
na.itp                 na.pdb
pf6.itp                pf6.pdb
ec.itp                 ec.pdb
dmc.itp                dmc.pdb
em.mdp
prep.mdp
npt_eq.mdp
prod_01.mdp
```

Their roles are:

| File type | Purpose |
|---|---|
| `.pdb` | One-molecule coordinate templates that Packmol copies and positions |
| `.itp` | GROMACS molecule topologies: atom types, charges, bonds, angles, and other force-field terms |
| `custom_atomtypes.itp` | Additional atom types required by the molecular topologies |
| `.mdp` | Settings for minimization, preparation, equilibration, and production |
| `.py` | Calculates the composition and creates a self-contained simulation case |

Verify the files:

```bash
cd ~/projects/napf6-ec-dmc
ls -1 electrolyte_generator_v2.py *.pdb *.itp *.mdp
```

Do not casually rename atoms or reorder atoms in the PDB or ITP files. The
coordinate order for each molecule must exactly match its topology order.

## 4. Generate the paper test case

Run:

```bash
cd ~/projects/napf6-ec-dmc
python3 electrolyte_generator_v2.py --paper-case
```

The expected summary is:

```text
Mode: paper
Requested EC wt%: 30
Requested concentration (mol/L): 1.5
Input density (g/cm^3): 1.329
EC molecules: 76
DMC molecules: 173
NaPF6 pairs: 30
Realized EC wt%: 30.04421848
Realized initial concentration (mol/L): 1.45964846
Box length (angstrom): 32.43700089
Box length (nm): 3.24370009
Initial volume (nm^3): 34.12888327
Expected atoms: 3076
```

The requested and realized values differ slightly because a simulation box can
contain only whole molecules. The `--paper-case` option reproduces the integer
counts reported for this composition in the paper: 30 Na+, 30 PF6-, 76 EC,
and 173 DMC.

The generated case is located at:

```text
cases/paper_wEC_30_C_1.5M/
```

This location is anchored beside `electrolyte_generator_v2.py`; it does not
change when the generator is launched from a different terminal directory.

### Other supported input modes

Generate counts for a target composition, concentration, and density:

```bash
python3 electrolyte_generator_v2.py \
    --w-ec 30 \
    --concentration 1.5 \
    --density 1.329
```

The default target calculation uses approximately 250 solvent molecules. Set a
different scale with `--total-solvents`, for example:

```bash
python3 electrolyte_generator_v2.py \
    --w-ec 30 \
    --concentration 1.5 \
    --density 1.329 \
    --total-solvents 500 \
    --case-label test_500_solvents
```

Specify exact integer counts instead:

```bash
python3 electrolyte_generator_v2.py \
    --n-ec 76 \
    --n-dmc 173 \
    --salt-pairs 30 \
    --density 1.329
```

Use `python3 electrolyte_generator_v2.py --help` to list all options.

## 5. General idea of the generator

The generator performs four connected jobs.

### 5.1 Calculate a discrete composition

For target mode, the EC/DMC molecule counts are chosen to approximate the
requested EC mass fraction. The salt count is calculated while including both
solvent mass and salt mass in the density/volume relation. This coupling is
important: adding NaPF6 increases the total mass and therefore changes the
initial volume associated with the target density.

In simplified form:

```text
mass of contents = mass(EC) + mass(DMC) + mass(NaPF6)
initial volume   = mass of contents / input density
box length       = cube root(initial volume)
concentration    = moles of NaPF6 / box volume in litres
```

The molecular weights used by the code are:

| Species | Molecular weight (g/mol) |
|---|---:|
| EC | 88.062 |
| DMC | 90.078 |
| NaPF6 | 167.953 |

### 5.2 Validate the templates

Before writing a case, the generator checks that:

- every required PDB, ITP, and MDP template exists;
- each PDB atom-name order matches the corresponding ITP atom-name order; and
- the molecular charges are +1 for Na, -1 for PF6, and zero for EC and DMC.

It also normalizes the DMC residue name, makes the PF6 template geometry
consistent with its topology bond length, and removes extra box/model records
from the one-molecule PDB templates.

### 5.3 Generate Packmol and GROMACS inputs

The code creates:

- `packmol_<case-label>.inp`, containing the molecule counts and initial cubic
  box in angstrom;
- `system_<case-label>.top`, containing the topology include order and matching
  molecule counts;
- copies of all validated PDB, ITP, and MDP files; and
- `run_workflow_<case-label>.sh`, which runs one simulation stage at a time and
  can also create the standard RDF/coordination outputs.

The workflow script changes into its own case directory before running any
command. Therefore Packmol, GROMACS, log, trajectory, checkpoint, and RDF files
remain in that case folder even if the script is invoked from the project root
or by its absolute path.

Packmol/PDB coordinates use **angstrom**, while GROMACS coordinates use
**nanometres**. The workflow script performs the conversion when it writes
`initial_box_<case-label>.gro`.

### 5.4 Record reproducibility information

Each case also contains:

- `generation_summary_<case-label>.txt`, a human-readable input/output summary;
- `manifest_<case-label>.json`, the same information in a machine-readable
  form; and
- `README_<case-label>.md`, containing case-specific instructions.

All results that depend on composition include the case label. Reusable
molecular and simulation templates retain their ordinary names, such as
`dmc.itp`, `ec.pdb`, and `npt_eq.mdp`.

These files preserve both the requested values and the realized values after
integer rounding.

### 5.5 Summary

1. The Python generator reads the requested composition.

   - In target mode, it calculates EC, DMC, and NaPF6 counts; EC + DMC defaults
     to 250.
   - In exact mode, it uses the molecule counts supplied by the user.
   - In paper mode, it uses predefined counts from the paper.

2. It calculates the initial box size from the total molecular mass and input
   density.
3. It validates the existing PDB molecular structures and matching ITP
   force-field files.
4. It creates `packmol_<case-label>.inp`, which tells Packmol how many copies of
   each molecule to pack and the box dimensions. Packmol later performs the
   actual packing.
5. It creates `system_<case-label>.top`, which tells GROMACS which molecular
   topology files to use and how many molecules are present.
6. It copies the validated PDB, ITP, and MDP templates into the new case
   directory. It does not generate ITP parameters from the molecules.
7. It generates `run_workflow_<case-label>.sh` to run each Packmol/GROMACS and
   standard RDF-analysis stage. The script automatically works inside its own
   case directory, keeping outputs separated from every other case.
8. It records the case information in `README_<case-label>.md`,
   `generation_summary_<case-label>.txt`, and `manifest_<case-label>.json`.

## 6. Inspect the generated case before running it

Enter the case directory and list its contents:

```bash
cd ~/projects/napf6-ec-dmc/cases/paper_wEC_30_C_1.5M
pwd
ls -1
```

Check the GROMACS molecule counts:

```bash
grep -A6 '\[ molecules \]' system_paper_wEC_30_C_1.5M.top
```

Expected counts:

```text
NA_M          30
PF6_M         30
EC_M          76
DMC_M         173
```

Check the Packmol structures, counts, and box:

```bash
grep -E 'structure|number|inside box' packmol_paper_wEC_30_C_1.5M.inp
```

The four molecule counts must match
`system_paper_wEC_30_C_1.5M.top`, and each maximum box coordinate should be
approximately `32.43700` angstrom.

## 7. Run Packmol

Make the workflow script executable if necessary, then run only the Packmol
stage:

```bash
chmod +x run_workflow_paper_wEC_30_C_1.5M.sh
./run_workflow_paper_wEC_30_C_1.5M.sh packmol
```

This stage:

1. runs Packmol with `packmol_paper_wEC_30_C_1.5M.inp`;
2. saves Packmol output in `packmol_paper_wEC_30_C_1.5M.log`;
3. verifies that `initial_box_paper_wEC_30_C_1.5M.pdb` contains exactly 3076
   atoms; and
4. creates `initial_box_paper_wEC_30_C_1.5M.gro` with a 3.24370009 nm cubic
   GROMACS box.

Verify the outputs:

```bash
grep -cE '^(ATOM|HETATM)' initial_box_paper_wEC_30_C_1.5M.pdb
head -n 2 initial_box_paper_wEC_30_C_1.5M.gro
tail -n 1 initial_box_paper_wEC_30_C_1.5M.gro
```

For this case, the first command and the second line of the GRO file
should both report `3076`. The final line should contain three box lengths near
`3.24370009` nm.

Packmol creates a non-overlapping initial arrangement; it does not equilibrate
the liquid or establish the final simulated density.

## 8. Energy minimization

Run:

```bash
./run_workflow_paper_wEC_30_C_1.5M.sh em
```

This creates `em_paper_wEC_30_C_1.5M.tpr`, then runs steepest-descent
minimization and writes files with the case-labeled `em_paper_wEC_30_C_1.5M`
prefix.

Check the end of the log:

```bash
tail -n 30 em_paper_wEC_30_C_1.5M.log
```

The important outcome is successful convergence without a fatal error. In the
first test, minimization converged in 254 steps with a maximum force below the
1000 kJ mol-1 nm-1 stopping criterion. Exact values can change because Packmol
uses a random initial arrangement.

## 9. Ten-picosecond preparation

Run:

```bash
./run_workflow_paper_wEC_30_C_1.5M.sh prep
```

This stage starts from the labeled minimized structure, generates velocities at
298.15 K, and performs the 10 ps preparation run. It creates labeled `prep`
structure and checkpoint files for the next stage.

The MDP intentionally follows the translated test protocol. GROMACS warns that
the Berendsen thermostat does not generate the correct canonical kinetic-energy
distribution. The script permits the one expected warning for this stage. Do
not add a larger `-maxwarn` value to hide a new or unexplained warning.

Check completion:

```bash
tail -n 20 prep_paper_wEC_30_C_1.5M.log
```

## 10. One-nanosecond NPT equilibration

Run:

```bash
./run_workflow_paper_wEC_30_C_1.5M.sh npt
```

This stage continues from the labeled preparation structure and checkpoint and
allows the box volume and density to relax at 298.15 K and 1 bar. The script
uses eight CPU threads and a lower scheduling priority (`nice -n 10`) to reduce
interference with normal PC use.

The GROMACS translation uses multiple time stepping: a 0.25 fs base step and a
factor of 8 for the slower forces, equivalent to a 2 fs outer step. Settings
such as `nstcalcenergy` must be multiples of 8 when this MTS factor is used.

Check completion:

```bash
tail -n 20 npt_eq_paper_wEC_30_C_1.5M.log
```

Extract temperature, pressure, volume, and density from the equilibrated part
of the trajectory. For example, to analyze the final 500 ps of a 1 ns run:

```bash
gmx energy \
    -f npt_eq_paper_wEC_30_C_1.5M.edr \
    -b 500 \
    -o npt_equilibrated_properties_paper_wEC_30_C_1.5M.xvg
```

At the interactive prompt, select:

```text
Temperature Pressure Volume Density
```

Then press Enter on an empty line to finish the selection.

Pressure in a small liquid box fluctuates strongly; the instantaneous pressure
does not need to stay near 1 bar. Judge equilibration mainly from the absence
of systematic drift in temperature, volume, and density over the chosen
analysis interval.

The first completed test gave an average density of about 1271 kg/m3 and an
average cubic box length of about 3.291 nm over its final 500 ps. Those values
are a reference from one run, not hard-coded pass/fail requirements.

### Initial density versus simulated density

The input density is an experimental or literature estimate used only to make a
reasonable **initial** Packmol box. For the paper case, 1.329 g/cm3 corresponds
to the paper's tabulated/interpolated experimental density for the nominal
30:70 wt%, 1.5 M system. During NPT equilibration, GROMACS changes the box
volume. Therefore, report the equilibrated simulated density and recompute the
realized concentration from the equilibrated volume rather than assuming the
input value remained fixed.

## 11. The 0.5 ns test production run

Run:

```bash
./run_workflow_paper_wEC_30_C_1.5M.sh prod
```

This continues from the labeled NPT structure and checkpoint and writes files
with the prefix `prod_01_paper_wEC_30_C_1.5M`, including the production
trajectory.

Check completion and performance:

```bash
tail -n 30 prod_01_paper_wEC_30_C_1.5M.log
```

The 0.5 ns trajectory is intended for a preliminary structural result such as
an RDF and first-shell coordination number. Do not use it to claim converged
viscosity, conductivity, or diffusion coefficients.

## 12. Preliminary Na-carbonyl oxygen RDF analysis

In the validated topology, the EC carbonyl oxygen is named `O00` and the DMC
carbonyl oxygen is named `O03`. Calculate their RDFs and running coordination
numbers separately.

Run the automated analysis stage:

```bash
./run_workflow_paper_wEC_30_C_1.5M.sh analysis
```

It reads `prod_01_paper_wEC_30_C_1.5M.xtc` and creates:

- `rdf_Na_ECcarbonylO_paper_wEC_30_C_1.5M.xvg`
- `coordination_Na_ECcarbonylO_paper_wEC_30_C_1.5M.xvg`
- `rdf_Na_DMCcarbonylO_paper_wEC_30_C_1.5M.xvg`
- `coordination_Na_DMCcarbonylO_paper_wEC_30_C_1.5M.xvg`

The script uses Na as the reference, EC atom `O00` for the EC selection, and
DMC atom `O03` for the DMC selection. Keeping the case label in every analysis
filename prevents Xmgrace from accidentally opening data from another case.

The RDF, `g(r)`, measures how strongly a selected atom is distributed around
Na+ relative to the bulk. The `-cn` output is the **running/cumulative
coordination number**. It is expected to keep increasing after later shells are
included. The first-shell coordination number is the value of that curve at
the first minimum of the corresponding RDF, not its value at the largest
plotted distance.

In the first short test, the Na-carbonyl-O first peak was near 0.235 nm and the
first minimum/plateau was near 0.36 nm. At 0.36 nm, the preliminary EC and DMC
contributions were approximately 1.53 and 2.38, respectively, for a total of
about 3.91. These are short-trajectory observations to be tested for stability,
not final reproduced values.

## 13. Plot with Xmgrace

Install Grace if needed:

```bash
sudo apt update
sudo apt install grace
```

Open an RDF:

```bash
xmgrace rdf_Na_ECcarbonylO_paper_wEC_30_C_1.5M.xvg
```

Open a running coordination-number file:

```bash
xmgrace coordination_Na_ECcarbonylO_paper_wEC_30_C_1.5M.xvg
```

Useful labels are:

- RDF x-axis: `r (nm)`
- RDF y-axis: `g(r)`
- RDF title: `Na+–carbonyl O radial distribution`
- coordination x-axis: `r (nm)`
- coordination y-axis: `Running coordination number, N(r)`

To save an editable Grace project, use **File -> Save As**. In the save dialog,
enter the complete filename in the bottom **Selection** field, for example:

```text
/home/nguyen/projects/napf6-ec-dmc/cases/paper_wEC_30_C_1.5M/rdf_Na_ECcarbonylO_paper_wEC_30_C_1.5M.agr
```

The top **Filter** field is a filename filter, not the save destination. Putting
the output path there can produce the “not a regular file” error.

Export the saved Grace project to PNG from the terminal:

```bash
gracebat \
    -hdevice PNG \
    -printfile rdf_Na_ECcarbonylO_paper_wEC_30_C_1.5M.png \
    rdf_Na_ECcarbonylO_paper_wEC_30_C_1.5M.agr
```

## 14. Files to preserve for reproducibility

For every reported system, preserve at least:

- generator version or Git commit;
- generator command and `generation_summary_<case-label>.txt`;
- `manifest_<case-label>.json`;
- PDB, ITP, TOP, and MDP inputs;
- the case-labeled Packmol input and log;
- `run_workflow_<case-label>.sh`;
- GROMACS `.log` files;
- final structures/checkpoints needed to continue the run;
- analysis commands, `.xvg` data, and exported figures; and
- GROMACS and Packmol versions.

Large binary simulation outputs can exceed GitHub's practical limits. Do not
commit every trajectory and checkpoint to a normal Git repository. A typical
`.gitignore` can exclude the largest/restartable binary outputs inside generated
case directories while retaining the text logs, inputs, and selected
structures:

```gitignore
cases/**/*.xtc
cases/**/*.trr
cases/**/*.edr
cases/**/*.tpr
cases/**/*.cpt
```

Store excluded reproducibility data in an appropriate research-data archive or
Git LFS and document its location. Do not delete the only copy merely because
it is excluded from Git.

## 15. Recording changes with GitHub

After editing this guide or an input file:

```bash
git status
git diff
git add README.md electrolyte_generator_v2.py
git commit -m "Document electrolyte generation workflow"
git push
```

Add other validated source files to the `git add` command as appropriate. Read
`git diff` before committing so that the commit records only the intended
change.

For each new scientific test, record the exact generator command, case label,
software versions, changes to MDP files, run length, and analysis interval in
this document or in a case-specific notes file.

## 16. Troubleshooting

### `Error: Missing templates`

The generator cannot find one or more required PDB/ITP/MDP files. Download the
complete repository and run the generator from its root, or provide the
template directory explicitly:

```bash
python3 electrolyte_generator_v2.py \
    --paper-case \
    --template-dir /path/to/templates
```

### `Output directory is not empty`

The generator avoids overwriting an existing case. Use a new `--case-label` or
inspect the existing directory. Avoid `--force` for a case that already has
simulation results, because mixing regenerated inputs with older results can
destroy the provenance of that run.

### `packmol: command not found`

Packmol is not available to Bash. Check:

```bash
command -v packmol
```

Fix the Julia Packmol wrapper or WSL `PATH` before continuing. Do not run the
Packmol input by double-clicking it; enter the generated case directory and run
the case-labeled workflow script with the `packmol` stage.

### Packmol reports an input or structure error

Confirm that the terminal is in the generated case directory and that
`na.pdb`, `pf6.pdb`, `ec.pdb`, and `dmc.pdb` are present. Packmol paths inside
the case-labeled Packmol input are relative to the current directory.

### GROMACS reports `No such file or directory` for an ITP

Run GROMACS from inside the generated case directory. Confirm that every
included file named near the top of `system_<case-label>.top` exists there. The
standard GROMACS `oplsaa.ff` directory is supplied by the GROMACS installation;
the other ITP files are supplied by this project.

### GROMACS stops because of warnings

Read each warning. The workflow script permits only the known warning counts
associated with the translated Berendsen/MTS settings. A different warning may
indicate an actual input inconsistency and should not be bypassed by increasing
`-maxwarn`.

### Units look inconsistent

Packmol and the PDB templates use angstrom. GROMACS uses nanometres. For the
paper case, `32.437 angstrom = 3.2437 nm`. The generated workflow performs this
factor-of-ten conversion.

## 17. Reference and software documentation

- A. V. Borovyk, Y. V. Kolesnik, and O. M. Kalugin, *Structure and transport
  properties of NaPF6 solutions in ethylene carbonate/dimethyl carbonate
  mixtures for sodium-ion batteries: MD simulation*, 2024,
  [DOI 10.26565/2220-637X-2024-43-03](https://doi.org/10.26565/2220-637X-2024-43-03).
- [Packmol user guide](https://m3g.github.io/packmol/userguide.shtml)
- [GROMACS installation guide](https://manual.gromacs.org/documentation/current/install-guide/index.html)
- [VS Code: Developing in WSL](https://code.visualstudio.com/docs/remote/wsl)
