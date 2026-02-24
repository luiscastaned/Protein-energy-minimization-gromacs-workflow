# Protein-energy-minimization-gromacs-workflow
Energy minimization workflow of a ligand-free protein using GROMACS for molecular docking and MD preparation.
📌 Overview

This repository contains a complete and simplified workflow for performing energy minimization of a ligand-free protein using GROMACS.

The commands shown here are intentionally synthetized and standardized for clarity and reproducibility.
For simplicity, the only required input file is named:
protein.pdb
This unified naming convention allows easier replication of the workflow without modifying multiple file names.

The minimized structure can be used for:

Molecular docking

Molecular dynamics simulations

Structural visualization

Further computational refinement

🛠 Software Requirements

GROMACS 2022+

Linux / Ubuntu / WSL

Basic terminal knowledge

Check your installation:

gmx --version
📂 Workflow Steps
1️⃣ Generate Topology
gmx pdb2gmx -f protein.pdb -o protein.gro -water tip3p
🔎 What this command does

Converts protein.pdb into GROMACS format.

Generates topology file (topol.top).

Adds hydrogens.

Assigns a force field.

Defines water model (TIP3P).

🧠 Why this step is important

This prepares the protein structure for molecular simulations.

2️⃣ Define the Simulation Box
gmx editconf -f protein.gro -o protein_box.gro -c -d 1.0 -bt cubic
🔎 What this command does

Centers the protein in the box (-c).

Adds 1.0 nm distance from the protein to box edges (-d 1.0).

Creates a cubic box (-bt cubic).

🧠 Why this step is important

Prevents artificial interactions due to periodic boundary conditions.

3️⃣ Solvate the System
gmx solvate -cp protein_box.gro -cs spc216.gro -o protein_solv.gro -p topol.top
🔎 What this command does

Fills the simulation box with water molecules.

Updates topology file with solvent information.

🧠 Why this step is important

Proteins function in aqueous environments; simulations must reflect this.

4️⃣ Create Energy Minimization Parameters

Create a file named:

minim.mdp

With the following content:

integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000
cutoff-scheme = Verlet
nstlist     = 1
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
🔎 Parameter meaning

steep → Steepest descent minimization

emtol → Convergence threshold

nsteps → Maximum number of steps

PME → Accurate electrostatics calculation

pbc → Periodic boundary conditions

5️⃣ Neutralize the System (Add Ions)

Preprocess:

gmx grompp -f minim.mdp -c protein_solv.gro -p topol.top -o ions.tpr -maxwarn 1

Add ions:

gmx genion -s ions.tpr -o protein_ions.gro -p topol.top -pname NA -nname CL -neutral

When prompted, select:

SOL
🔎 What this step does

Neutralizes total system charge.

Replaces water molecules with Na⁺ or Cl⁻ ions.

🧠 Why it matters

Electrostatic stability requires a neutral simulation box.

6️⃣ Run Energy Minimization

Prepare binary:

gmx grompp -f minim.mdp -c protein_ions.gro -p topol.top -o em.tpr

Run minimization:

gmx mdrun -v -deffnm em
🔎 What this does

Executes steepest descent minimization.

Produces output files:

em.gro

em.edr

em.log

7️⃣ Extract the Minimized Structure
gmx trjconv -s em.tpr -f em.gro -o protein_min.pdb

Select:

Protein
📌 Output
protein_min.pdb

This structure is ready for:

Docking studies

Molecular dynamics

Structural visualization

📊 Extract Potential Energy
gmx energy -f em.edr -o potential.xvg

Select:

Potential
Convert to CSV (Optional)
grep -v '^[@#]' potential.xvg | tr -s ' ' ',' > potential.csv

This removes headers and converts the file into a comma-separated format for Excel or plotting software.

📈 Expected Energy Profile

A successful minimization should show:

Rapid decrease in potential energy (steric clash correction)

Energy plateau indicating structural stabilization

📌 Final Output Files

protein_min.pdb

potential.xvg

potential.csv

em.log

📜 License

MIT License

🔬 Citation Metadata

Create a file named CITATION.cff:

cff-version: 1.2.0
title: Protein Energy Minimization using GROMACS
authors:
  - family-names: Castañeda
    given-names: Luis Ernesto
date-released: 2026
🧠 Notes on Simplification

The commands in this repository are presented in a synthetized format to improve clarity and reproducibility.

To avoid confusion, the workflow assumes:

A single input file named protein.pdb

Standard water model (TIP3P)

Default ion naming (NA / CL)

This makes the protocol easier to replicate and ideal for teaching, reproducible workflows, and portfolio presentation
