# Protein-energy-minimization-gromacs-workflow
Energy minimization workflow of a ligand-free protein using GROMACS for molecular docking and MD preparation.

# 🧬 Protein Energy Minimization Workflow using GROMACS

---

## 📌 Overview

This repository provides a **step-by-step workflow** for performing energy minimization of a ligand-free protein using GROMACS.

To simplify reproducibility, this protocol assumes a single standardized input file:

```
protein.pdb
```

All commands are intentionally presented in a synthetized format to improve clarity, reproducibility, and ease of teaching.

The minimized structure can be used for:

- Molecular docking  
- Molecular dynamics simulations  
- Structural refinement  
- Structural visualization  

---

# 🚀 Step-by-Step Workflow

---

## 🔹 Step 1 — Generate Topology

```bash
gmx pdb2gmx -f protein.pdb -o protein.gro -water tip3p
```

### What this step does

- Converts `protein.pdb` into GROMACS format  
- Generates the topology file (`topol.top`)  
- Adds hydrogens  
- Assigns a force field  
- Defines TIP3P water model  

This prepares the structure for simulation.

---

## 🔹 Step 2 — Define the Simulation Box

```bash
gmx editconf -f protein.gro -o protein_box.gro -c -d 1.0 -bt cubic
```

### What this step does

- Centers the protein in the box  
- Adds 1.0 nm distance from protein to box edges  
- Creates a cubic box  

This prevents artificial periodic interactions.

---

## 🔹 Step 3 — Solvate the System

```bash
gmx solvate -cp protein_box.gro -cs spc216.gro -o protein_solv.gro -p topol.top
```

### What this step does

- Fills the box with water molecules  
- Updates the topology file  

Proteins require aqueous environments for realistic simulations.

---

## 🔹 Step 4 — Define Energy Minimization Parameters

Create a file named:

```
minim.mdp
```

With the following content:

```bash
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
```

### What these parameters control

- Steepest descent minimization  
- Convergence tolerance  
- Maximum number of steps  
- Electrostatics treatment (PME)  
- Periodic boundary conditions  

---

## 🔹 Step 5 — Neutralize the System

### Preprocess:

```bash
gmx grompp -f minim.mdp -c protein_solv.gro -p topol.top -o ions.tpr -maxwarn 1
```

### Add ions:

```bash
gmx genion -s ions.tpr -o protein_ions.gro -p topol.top -pname NA -nname CL -neutral
```

When prompted, select:

```
SOL
```

This replaces water molecules with Na⁺/Cl⁻ ions and neutralizes system charge.

---

## 🔹 Step 6 — Run Energy Minimization

### Prepare binary:

```bash
gmx grompp -f minim.mdp -c protein_ions.gro -p topol.top -o em.tpr
```

### Run minimization:

```bash
gmx mdrun -v -deffnm em
```

### Output files generated

- `em.gro`  
- `em.edr`  
- `em.log`  

---

## 🔹 Step 7 — Extract the Minimized Structure

```bash
gmx trjconv -s em.tpr -f em.gro -o protein_min.pdb
```

Select:

```
Protein
```

### Final minimized structure:

```
protein_min.pdb
```

This file is ready for docking or molecular dynamics.

---

# 📊 Extract Potential Energy

```bash
gmx energy -f em.edr -o potential.xvg
```

Select:

```
Potential
```

---

## Optional: Convert to CSV

```bash
grep -v '^[@#]' potential.xvg | tr -s ' ' ',' > potential.csv
```

This removes headers and converts the file for plotting in Excel or other software.

---

# 📈 Expected Energy Profile

A successful minimization should show:

- Rapid energy decrease (steric clash correction)  
- Energy plateau (structural stabilization)  

---

# 📌 Final Output

- `protein_min.pdb`  
- `potential.xvg`  
- `potential.csv`  
- `em.log`  

---

# 📜 License

MIT License

---

# 🔬 Citation

Create `CITATION.cff`:

```yaml
cff-version: 1.2.0
title: Protein Energy Minimization using GROMACS
authors:
  - family-names: Castañeda
    given-names: Luis Ernesto
date-released: 2026
```

---

# 🧠 Notes

This workflow is intentionally simplified:

- Single standardized input file: `protein.pdb`  
- Default water model: TIP3P  
- Standard ion naming: NA / CL  

Designed for clarity, reproducibility, and educational purposes.
