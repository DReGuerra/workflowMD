# LAMMPS Installation and Configurations
I have compiled several installation configurations to enhance LAMMPS capabilitties and/or performance. I have tested these installations on a Windows 10 WSL2 Debian. These configurations are transferrable for a native Linux installation, but this is also a validation of their implementation in WSL. This section focuses on the installation configurations for LAMMPS specifically.

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
