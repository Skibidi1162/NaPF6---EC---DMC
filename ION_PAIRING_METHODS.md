# Literature basis for NaPF6 ion-pairing analysis

This note fixes the scientific definitions used by the standalone trajectory
analysis command in this repository. It is intentionally separate from the
simulation generator.

## Primary references

1. A. V. Borovyk, Y. V. Kolesnik, and O. M. Kalugin, "Structure and
   Transport Properties of NaPF6 Solutions in Mixtures of Ethylene Carbonate
   with Dimethyl Carbonate for Sodium-Ion Batteries: MD Simulation,"
   *Kharkiv University Bulletin. Chemical Series* 43 (2024).
   DOI: https://doi.org/10.26565/2220-637X-2024-43-03

   Method adopted here: radial distribution functions (RDFs) are integrated
   to obtain running coordination numbers (CNs); the first RDF minimum defines
   a coordination-shell boundary. The phosphorus atom is used as the
   molecular coordination center of PF6-. This matches the atom naming in
   `pf6.itp`, where every PF6- residue has one atom named `P`.

2. D. A. Rakov et al., "The impact of electrode conductivity on electrolyte
   interfacial structuring and its implications on the Na0/+ electrochemical
   performance," *Energy & Environmental Science* 16 (2023) 3919-3931.
   DOI: https://doi.org/10.1039/D3EE00864A
   Supporting information:
   https://www.rsc.org/suppdata/d3/ee/d3ee00864a/d3ee00864a1.pdf

   Method adopted here: for 1.0 M NaPF6 in EC:DMC (1:1 by volume), the study
   calculates Na-P RDFs and CNs using the phosphorus atom of PF6-. Its reported
   first-shell analysis evaluates the Na coordination environment at 0.405 nm.
   The present command therefore reports ion populations at 0.405 nm as an
   EC/DMC-specific reference sensitivity calculation.

3. H. Lee, J.-H. Shim, and S. Lee, "Solvation structure and ion transport in
   mixed ethylene carbonate and dimethyl carbonate electrolytes for sodium-ion
   batteries: a molecular dynamics study," *Molecular Simulation* 51 (2025)
   1110-1123.
   DOI: https://doi.org/10.1080/08927022.2025.2586506

   Method adopted here: Na-P RDFs are used to monitor direct Na+-PF6-
   association across EC/DMC compositions, and Na-P coordination numbers are
   evaluated with a unified 0.45 nm first-shell boundary. The paper also tests
   0.35 nm, the smallest observed first-minimum position, and reports that the
   smaller boundary changes absolute coordination numbers only marginally
   without changing composition-dependent ordering. The command therefore
   includes both 0.35 and 0.45 nm as EC/DMC-specific sensitivity cutoffs. Its
   primary SIP/CIP/AGG definition remains the phosphorus-count method in
   reference 4 below.

   The paper additionally uses Na-Na RDFs as a qualitative salt-dissociation
   diagnostic. That cation-cation observable is not substituted for the Na-P
   neighbor counts used to define SIP/CIP/AGG in this project.

4. "Directing selective solvent presentations at electrochemical interfaces
   to enable initially anode-free sodium metal batteries," *Nature
   Communications* 16, 8265 (2025).
   DOI: https://doi.org/10.1038/s41467-025-63902-4

   Method adopted here: the number of PF6- phosphorus atoms within 0.5 nm of
   every Na+ is counted over the production trajectory. A Na+ environment is
   classified as:

   - SIP/SSIP: zero P atoms inside 0.5 nm;
   - CIP: one P atom inside 0.5 nm;
   - AGG: two or more P atoms inside 0.5 nm.

   The paper also checks 0.4 and 0.6 nm cutoffs. The command keeps the published
   0.5 nm cutoff as its primary result, makes it configurable, and reports
   sensitivity results at 0.35, 0.4, 0.405, 0.45, and 0.6 nm by default.

5. D. Monti et al., "Towards standard electrolytes for sodium-ion batteries:
   physical properties, ion solvation and ion-pairing in alkyl carbonate
   solvents," *Physical Chemistry Chemical Physics* 22 (2020) 22768-22777.
   DOI: https://doi.org/10.1039/D0CP03639K

   Context adopted here: NaPF6 ion association and Na+ solvation in EC, DMC,
   and their mixtures depend on salt concentration and solvent composition.
   Results must therefore remain case-specific and must retain the generated
   case label in every output path.

## Operational analysis specification

- Primary RDF and molecular CN: `resname NA and name NA` around
  `resname PF6 and name P`.
- Supporting atom-level RDF and CN: Na around all six PF6 fluorine atoms.
  The Na-F CN is an atom count and is not interpreted as a PF6- molecule count.
- RDF and CN calculations use `gmx rdf`, periodic boundary conditions, and the
  production TPR/XTC belonging to one generated case.
- SIP/CIP/AGG counts use dynamic GROMACS selections with periodic boundary
  conditions. For every trajectory frame and every Na+, the command counts P
  atoms inside each configured cutoff.
- The primary population fractions use the published 0.5 nm rule. Alternative
  cutoffs are sensitivity checks, not silently substituted definitions.
- Fixed sensitivity cutoffs retain their literature roles: 0.405 nm is the
  Rakov EC:DMC (1:1) first-shell distance; 0.45 nm is the Lee composition-wide
  EC/DMC boundary; and 0.35, 0.4, and 0.6 nm test cutoff robustness.
- The first minimum after the first Na-P RDF peak is detected and recorded as
  a data-derived coordination-shell estimate. It is reported separately from
  the fixed literature cutoffs.
- All outputs are written below
  `cases/<case-label>/analysis_ion_pairing/`. The project root and other case
  directories are never used for trajectory-analysis output.
- Fractions are calculated over all analyzed Na environments (number of
  analyzed frames multiplied by the number of Na+ ions). Frame-level
  populations are also retained so convergence and block uncertainty can be
  assessed.

## Interpretation limits

The SIP/SSIP label in the adopted population rule means that no PF6-
phosphorus atom lies inside the specified Na+ solvation-shell cutoff. It does
not distinguish a truly solvent-separated pair from a completely dissociated
ion. That distinction would require a second-shell or explicit solvent-bridge
definition and is outside the primary method reproduced here.

The existing 0.5 ns production trajectory is suitable for testing this
analysis pipeline, but final population claims require convergence checks and
substantially longer or replicated production sampling.
