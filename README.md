# Protein-energy-minimization-gromacs-workflow
Energy minimization workflow of a ligand-free protein using GROMACS for molecular docking and MD preparation.
# 🧬 Protein Energy Minimization in GROMACS

---

## 📌 Overview

This repository contains a simplified and reproducible workflow for performing **energy minimization of a ligand-free protein** using **GROMACS**.

The commands are intentionally presented in a **synthetized format** for clarity and ease of use.  
To simplify the workflow and avoid unnecessary renaming, the only required input file is standardized as:

protein.pdb

This unified naming convention ensures easier replication and teaching purposes.

The minimized structure obtained can be used for:

- Molecular docking  
- Molecular dynamics simulations  
- Structural visualization  
- Further computational refinement  

---

## 🛠 Software Requirements

- GROMACS 2022+
- Linux / Ubuntu / WSL
- Basic terminal knowledge

Check your GROMACS installation:

```bash
gmx --version
