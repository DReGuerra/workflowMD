# WorkflowMD
André Guerra \
July, 2022 \
andre.guerra@mail.mcgill.ca  

---
### Description:
This is repository contains an introduction to molecular dynamics (MD) simulations using [PACKMOL](http://leandro.iqm.unicamp.br/m3g/packmol/examples.shtml), [Moltemplate](https://www.moltemplate.org/), and [LAMMPS](https://www.lammps.org/). I have tried to demonstrate a workflow for simulating a simple system found in `project/`. This workflow can be modified and extented as desired for your project. The system molecule coordinates are defined in `1_packmol`, and the PACKMOL software is used to populate a simulation box. Next, the system's data file is generated with Moltemplate in `2_moltemplate`. Finally, the LAMMPS simulation is conducted in `3_lammps`. This tree structure is designed to contain the workflow, segregate files used/produced by different software packages, and ultimately to organize our project space. The bash script `run.sh` executes the workflow.

Below you will see a list of the contents of the repository. The markdown file `main.md` will be the hub where all information regaring this molecular dynamics workshop is found. So, let's continue this over [there](https://github.com/DReGuerra/molecular_dynamics_workshop/blob/main/main.md).

---
## Core Contents
1. `project/` $\rightarrow$ the main sample MD project files to test workflowMD
2. `lammps-install-configs` $\rightarrow$ configurations and installation details for LAMMPS
3. `1_packmol/` $\rightarrow$ files associated with PACKMOL
4. `2_moltemplate/` $\rightarrow$ files associated with Moltemplate
5. `3_lammps/` $\rightarrow$ files associated with LAMMPS
6. `main.md` $\rightarrow$ this file contains the main content of this workshop
7. `run.sh` $\rightarrow$ bash script that runs the MD workflow
8. `run_cvmfs.sh` $\rightarrow$ modified bash script that runs the MD workflow on a cluster with [slurm scheduler](https://slurm.schedmd.com/documentation.html) and the [CVMFS](https://cernvm.cern.ch/fs/) stack (for reference and for my own testing)

## Tree Structure
<pre>
workflowMD/
├── lammps-install-configs/
│   ├── test/
│   │   ├── in.kokkos-test.lammps
│   │   └── log.lammps
│   └── README.md
├── project/
│   ├── 1_packmol/
│   │   ├── in.2pack_H2O.inp
│   │   ├── out_packed.xyz
│   │   └── water.xyz
│   ├── 2_moltemplate/
│   │   └── input_files/
│   │       ├── system_H2O.lt
│   │       ├── system.lt
│   │       └── water.tip4p-ice.lt
│   ├── 3_lammps/
│   │   ├── data/
│   │   ├── fix/
│   │   ├── restart/
│   │   ├── traj/
│   │   ├── log.npt.lammps
│   │   └── run.ag.lammps
│   ├── run_cvmfs.sh
│   └── run.sh
└── README.md
</pre>