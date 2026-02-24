# Protein-energy-minimization-gromacs-workflow
Energy minimization workflow of a ligand-free protein using GROMACS for molecular docking and MD preparation.

🧬 Protein Energy Minimization in GROMACS
📌 Overview

This repository contains a complete workflow for performing energy minimization of a ligand-free protein using GROMACS.

Energy minimization is the first essential step before:

Molecular docking

Molecular dynamics (MD) simulations

Structural refinement

Binding pocket analysis

The goal is to remove steric clashes and optimize atomic geometry before further computational analysis.

🛠 Software Requirements

GROMACS 2022+

Linux / Ubuntu / WSL

Basic terminal knowledge

Nano or any text editor

Check your GROMACS installation:

gmx --version
📂 Workflow Steps (Explained)
1️⃣ Generate Topology
gmx pdb2gmx -f protein.pdb -o protein.gro -water tip3p
🔎 What this does:

pdb2gmx converts a PDB file into a GROMACS-compatible structure.

-f protein.pdb → input protein structure.

-o protein.gro → output coordinate file.

-water tip3p → defines the water model.

🧠 Why it matters:

This step:

Assigns a force field

Adds hydrogens

Generates the topology file (topol.top)

Prepares the system for simulation

2️⃣ Define Simulation Box
gmx editconf -f protein.gro -o protein_box.gro -c -d 1.0 -bt cubic
🔎 What this does:

editconf defines the simulation box.

-c → centers the protein.

-d 1.0 → sets 1.0 nm distance from protein to box edge.

-bt cubic → creates a cubic box.

🧠 Why it matters:

Prevents protein atoms from interacting with their own periodic images.

3️⃣ Solvate the System
gmx solvate -cp protein_box.gro -cs spc216.gro -o protein_solv.gro -p topol.top
🔎 What this does:

solvate fills the box with water molecules.

-cp → protein box structure.

-cs → water configuration file.

-p → updates topology file.

🧠 Why it matters:

Biological systems function in aqueous environments.

4️⃣ Create Energy Minimization Parameters

Create file: mdp/minim.mdp

integrator = steep
emtol = 1000.0
emstep = 0.01
nsteps = 50000
cutoff-scheme = Verlet
nstlist = 1
ns_type = grid
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
🔎 Parameter explanation:

integrator = steep → Steepest descent algorithm.

emtol = 1000.0 → Convergence threshold (kJ/mol/nm).

emstep = 0.01 → Step size.

nsteps = 50000 → Maximum steps.

PME → Particle Mesh Ewald (electrostatics).

pbc = xyz → Periodic boundary conditions.

🧠 Why it matters:

Controls how energy minimization is performed.

5️⃣ Add Ions (Neutralization)

Prepare binary input:

gmx grompp -f minim.mdp -c protein_solv.gro -p topol.top -o ions.tpr -maxwarn 1
🔎 What this does:

grompp preprocesses input files.

Generates .tpr file (portable binary run input).

Now neutralize:

gmx genion -s ions.tpr -o protein_ions.gro -p topol.top -pname NA -nname CL -neutral
🔎 What this does:

Replaces water molecules with ions.

-neutral → neutralizes system charge.

Select group SOL when prompted.

🧠 Why it matters:

Simulations require neutral systems for proper electrostatics.

6️⃣ Run Energy Minimization

Prepare run:

gmx grompp -f minim.mdp -c protein_ions.gro -p topol.top -o em.tpr

Run minimization:

gmx mdrun -v -deffnm em
🔎 What this does:

mdrun executes simulation.

-deffnm em → sets default output names.

-v → verbose mode.

🧠 Expected output:

em.gro

em.edr

em.log

7️⃣ Extract Minimized Structure
gmx trjconv -s em.tpr -f em.gro -o protein_min.pdb

Select group: Protein

🔎 What this does:

Converts final minimized structure to .pdb format.

📌 Final Structure:

protein_min.pdb

Ready for:

Docking

MD

Visualization

📊 Extract Potential Energy
gmx energy -f em.edr -o potential.xvg

Select: Potential

Convert to CSV
grep -v '^[@#]' potential.xvg | tr -s ' ' ',' > potential.csv
🔎 What this does:

Removes header lines.

Converts spaces to commas.

Makes file Excel-compatible.

📈 Expected Energy Curve

A successful minimization shows:

Rapid energy decrease → steric clash correction

Energy plateau → structural stabilization

📌 Output Files

protein_min.pdb

potential.xvg

potential.csv

em.log

em.edr

📜 License

MIT License

🔬 Citation

Create file CITATION.cff

cff-version: 1.2.0
title: Protein Energy Minimization using GROMACS
authors:
  - family-names: Castañeda
    given-names: Luis Ernesto
date-released: 2026

This allows GitHub to display citation metadata automatically.

🚀 Optional Improvements

Add PNG energy plot

Add example PDB dataset

Add workflow diagram

Add Zenodo DOI

Add GROMACS version badge

🎯 Project Type

This repository can be oriented as:

📚 Academic project

🧠 Bioinformatics portfolio

🧪 Reproducible computational workflow

📄 Supplementary material for publication
