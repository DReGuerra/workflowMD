# LAMMPS Installation and Configurations
I have compiled several installation configurations to enhance LAMMPS capabilitties and/or performance. I have tested these installations on a Windows 10 WSL2 Debian. These configurations are transferrable for a native Linux installation, but this is also a validation of their implementation in WSL. This section focuses on the installation configurations for LAMMPS specifically.<br> 
<br>
For readability, I have grouped some bash commands in the same code block. But as a result, if you copy the whole block and paste into terminal, it will not like that. So, you in blocks that contain multiple lines of commands, you should copy them individually line-by-line. Finally commands that were too long, were broken up into multiple lines using `/`.

## Install MPI
MPI: Message Passing Interface is a library designed to facilitate parallelization of tasks/computations (e.g., LAMMPS). It shares data between processes that do not share physical memory (e.g., computations ocated in different nodes). 
```
sudo apt-get update
sudo apt-get install openmpi-bin libopenmpi-dev
```
Verify the installation:
```
which mpirun
```
Should return:
```
usr/bin/mpirun
```
And
```
mpirun --version
```
Should returrn:
```
mpirun (Open MPI) 4.1.4
```

## KOKKOS-CUDA-OMP
1. Install CUDA toolkit.

```
wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb

sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update
sudo apt-get install cuda-toolkit-12-8
```
Add the CUDA Toolkit to the PATH by adding this line to `.bashrc`
```
export PATH=/usr/local/cuda-12.8/bin:$PATH
```
Validate by getting nvcc version:
```
nvcc --version
```
Should return something like this:
```
nvcc: NVIDIA (R) Cuda compiler driver
Copyright (c) 2005-2025 NVIDIA Corporation
Built on Fri_Feb_21_20:23:50_PST_2025
Cuda compilation tools, release 12.8, V12.8.93
Build cuda_12.8.r12.8/compiler.35583870_0
```
2. Tools for Debian
```
sudo apt-get update
sudo apt-get install build-essential cmake git
```
3. Get LAMMPS source code:
```
wget https://lammps.org/tars/lammps-stable.tar.gz
tar xzf lammps-stable.tar.gz
```
4. Configure LAMMPS with KOKKOS and CUDA

First, prepare for installation. Replace `*` with whatever version you downloaded in the previous step.
```
cd lammps-*
mkdir build
cd build
```

Build LAMMPS with KOKKOS (acceleration package), OpenMP (multithreading CPUs), and CUDA (NVIDIA GPU) (for my specific RTX 3070 we use the AMPERE86 architecture).
```
cmake ../cmake \
  -D PKG_KOKKOS=on \
  -D Kokkos_ENABLE_OPENMP=on \
  -D Kokkos_ENABLE_CUDA=on \
  -D Kokkos_ARCH_AMPERE86=on \
  -D CMAKE_BUILD_TYPE=Release \
  -D BUILD_SHARED_LIBS=off \
  -D CMAKE_INSTALL_PREFIX=$HOME/lammps_kokkos_cuda
```
Add extra build configurations if desired, these are sample ones:
```
  -D PKG_EXTRADUMP=on \
  -D PKG_MOLECULE=on \
  -D PKG_RIGID=on 
```
Still in the same directory `build`:
```
make -j4
make install
```
In the configuration above, I specified the location of the install using the flag<br>
`-D CMAKE_INSTALL_PREFIX=$HOME/lammps_kokkos_cuda`<br>

This means that the executable created after the command `make install`, will be placed at `$HOME/lammps_kokkos_cuda/bin`, and it will be automatically named `lmp`. 

5. Change the executable name:
```
cd $HOME/lammps_kokkos_cuda/bin
mv lmp lmp_kokkos_cuda
```
Add this location to the `PATH` by adding the export statemnet at the end of your `.bashrc` file:
```
export PATH=$HOME/lammps_kokkos_cuda/bin:$PATH
```
Now you can call on the executable `lmp_kokkos_cuda` from anywhere.

## Test KOKKOS-CUDA LAMMPS
Create a simple test input script and run as follows. In this repo, these files are found in the directory `test`. 
```
cd ~/
mkdir test
cd test
touch in.kokkos-test.lammps
```
Use a text editor and add the following to the LAMMPS input file:
```
package kokkos
# simple LJ system
atom_style atomic
lattice fcc 0.8442
region box block 0 5 0 5 0 5
create_box 1 box
create_atoms 1 box
mass 1 1.0
velocity all create 1.0 1234

pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0 2.5

fix 1 all nve
run 1000
```
Run the test simulation:
```
lmp_kokkos_cuda -k on g 1 -sf kk < in.kokkos-test.lammps
```
Notice how the command above calls on the executable `lmp_kokkos_cuda` we made available in the `PATH`. It also uses flags to specify how to run the KOKKOS package: `-k on` turns the KOKKOS package on, `g 1` tells KOKKOS to use up to 1 GPU, `-sf kk` tells LAMMPS to automatically convert any possible `pair_style` or other commands to their KOKKOS accelerated versions.<br>

After running the simulation successfully, the `log.lammps` file will output this:
```
LAMMPS (29 Aug 2024 - Update 1)
KOKKOS mode with Kokkos version 4.3.1 is enabled (src/KOKKOS/kokkos.cpp:72)
  will use up to 1 GPU(s) per node
  using 1 OpenMP thread(s) per MPI task
package kokkos
package kokkos
# simple LJ system
atom_style atomic
lattice fcc 0.8442
Lattice spacing in x,y,z = 1.6795962 1.6795962 1.6795962
region box block 0 5 0 5 0 5
create_box 1 box
Created orthogonal box = (0 0 0) to (8.397981 8.397981 8.397981)
  1 by 1 by 1 MPI processor grid
create_atoms 1 box
Created 500 atoms
  using lattice units in orthogonal box = (0 0 0) to (8.397981 8.397981 8.397981)
  create_atoms CPU = 0.004 seconds
mass 1 1.0
velocity all create 1.0 1234

pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0 2.5

fix 1 all nve
run 1000
Generated 0 of 0 mixed pair_coeff terms from geometric mixing rule
Neighbor list info ...
  update: every = 1 steps, delay = 0 steps, check = yes
  max neighbors/atom: 2000, page size: 100000
  master list distance cutoff = 2.8
  ghost atom cutoff = 2.8
  binsize = 2.8, bins = 3 3 3
  1 neighbor lists, perpetual/occasional/extra = 1 0 0
  (1) pair lj/cut/kk, perpetual
      attributes: full, newton off, kokkos_device
      pair build: full/bin/kk/device
      stencil: full/bin/3d
      bin: kk/device
Per MPI rank memory allocation (min/avg/max) = 3.365 | 3.365 | 3.365 Mbytes
   Step          Temp          E_pair         E_mol          TotEng         Press     
         0   1             -6.7733681      0             -5.2763681     -5.3928057    
      1000   0.55667277    -6.1159279      0             -5.2825887     -1.8795406    
Loop time of 0.698243 on 1 procs for 1000 steps with 500 atoms

Performance: 618695.793 tau/day, 1432.166 timesteps/s, 716.083 katom-step/s
80.2% CPU use with 1 MPI tasks x 1 OpenMP threads

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.08718    | 0.08718    | 0.08718    |   0.0 | 12.49
Neigh   | 0.1079     | 0.1079     | 0.1079     |   0.0 | 15.45
Comm    | 0.16071    | 0.16071    | 0.16071    |   0.0 | 23.02
Output  | 0.00016291 | 0.00016291 | 0.00016291 |   0.0 |  0.02
Modify  | 0.15699    | 0.15699    | 0.15699    |   0.0 | 22.48
Other   |            | 0.1853     |            |       | 26.54

Nlocal:            500 ave         500 max         500 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:           1954 ave        1954 max        1954 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:              0 ave           0 max           0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
FullNghs:        37854 ave       37854 max       37854 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 37854
Ave neighs/atom = 75.708
Neighbor list builds = 91
Dangerous builds = 0
Total wall time: 0:00:01
```