# Molecular Dynamics Workshop
André Guerra \
July, 2022 \
andre.guerra@mail.mcgill.ca  

## Contents
1. Overview
2. Preparation
3. Workflow
4. Analyses

## 1. Overview
This markdown file will include the main contents of this workshop. Here, we will discuss any important information and describe some of the main aspects of the molecular dynamics (MD) workflow demonstrated in this repository. The main sections of `main.md` are listed above in the Contents section.

## 2. Preparation 
First, we will install the necessary software for executing our MD workflow. The following list is intended as a guide to get started. It may not be comprehensive for all purposes. Please refer to each software's webpage for details/requirements of installation. We will use [PACKMOL](http://leandro.iqm.unicamp.br/m3g/packmol/examples.shtml), [Moltemplate](https://www.moltemplate.org/), and [LAMMPS](https://www.lammps.org/).

For the time being, I will assume that you are using a Linux machine (or WSL in Windows10) and that you are familiar with terminal commands. In the future, I may incorporate more detailed instructions on setting up your machine for running this workflow.

### A. This repository
First, clone this repository locally where you usually keep repos (e.g., `~/repos/`). You can `cd` into the desired directory and then run the following command.
```
git clone git@github.com:DReGuerra/molecular_dynamics_workshop.git
```

### B. Auxilliary
1. Install miniconda package manager
    Download the executable installer for Linux x86_64.
    ```
    curl -sL \
    "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" > \
    "Miniconda3.sh"
    ```
    Run installer.
    ```
    bash Miniconda3.sh
    ```
2. Install git version control
    ```
    sudo apt install git
    ```
### C. LAMMPS
More details on LAMMPS can be found [here](https://www.lammps.org/).
1. Create a conda environment of LAMMPS 
    ```
    conda create -n lammps
    ```
2. Install LAMMPS
    Activate the conda environment
    ```
    conda activate lammps
    ```
    Install LAMMPS
    ```
    conda install lammps -c conda-forge
    ```
*NOTE:* From this point on, you need to activate the `lammps` conda environment before running LAMMPS.

3. Clone the LAMMPS git repository in a directory called lammps in home (i.e., `~/lammps`)
    ```
    git clone -b stable https://github.com/lammps/lammps.git ~/lammps
    ```
4. Test LAMMPS
    This will run the example script VISCOSITY provided in the LAMMPS repository. It will take a minute or two to completely run.
    ```
    cd ~/lammps/examples/VISCOSITY/ && lmp_serial -in in.wall.2d
    ```
If the example simulation ran correctly, at this point LAMMPS is fully installed. 

### D. Moltemplate
More details on Moltemplate can be found [here](https://www.moltemplate.org/).
1. Clone the Moltemplate git repository in a directory called moltemplate in home (i.e., `~/moltemplate`)
    ```
    git clone https://github.com/jewettaij/moltemplate ~/moltemplate
    ```
2. Update your system
    ```
    sudo apt update && sudo apt upgrade
    ```
3. Build moltemplate
    ```
    cd ~/moltemplate && pip install . --user
    ```
4. Export `.local/bin` to PATH and make it persistent across logins
    ```
    echo 'export PATH=~/.local/bin:$PATH' >> ~/.bashrc && source ~/.bashrc
    ```
5. Test Moltemplate
    ```
    cd ~/moltemplate/examples/all_atom/force_field_OPLSAA/waterSPC_using_OPLSAA/moltemplate_files/ \
    && moltemplate.sh system.lt
    ```
If the script ran without problems, at this point Moltemplate is fully installed.

### E. PACKMOL
More details on PACKMOL can be found [here](http://leandro.iqm.unicamp.br/m3g/packmol/examples.shtml).
1. Clone the git repository in a directory called packmol in home (i.e., `~/packmol`)
    ```
    git clone http://github.com/mcubeg/packmol ~/packmol
    ```
2. Install FORTRAN compiler
    ```
    sudo apt install gfortran
    ```
3. Compile PACKMOL
    ```
    cd ~/packmol && make
    ```
4. Download the PACKMOL example scripts
    ```
    cd ~/packmol && \
    curl -sL \
    "http://leandro.iqm.unicamp.br/m3g/packmol/examples/examples.tar.gz" > \
    examples.tar.gz
    ```
    Extract the examples from the tar archive
    ```
    tar -xvf examples.tar.gz
    ```
5. Export packmol to PATH and make it persistent across logins
    ```
    echo 'export PATH=~/packmol:$PATH' >> ~/.bashrc && source ~/.bashrc 
    ```
6. Test PACKMOL
    ```
    cd ~/packmol/examples/ && packmol < bilayer.inp
    ```
If the script ran without problems, at this point PACKMOL is fully installed.

## 3. Workflow
The molecular dynamics simulation presented here makes use of the `run.sh` script to control the workflow and execute the stages necessary to complete the simulation. 

### A. Script parameters
The `run.sh` bash script has input parameters which are handled in the "0. Prep" section of the script by the `while getops` block. Below we have this block from the script.
```
#  0. Prep
# handle flag(s)
while getopts g:w:d:e:s:r:n:f:p: flag
do
    case "${flag}" in
        g) gas=${OPTARG};;
        w) water=${OPTARG};;
        d) seed=${OPTARG};;
        e) ensemble=${OPTARG};;
        s) steps=${OPTARG};;
        r) restart=${OPTARG};;
        n) restartFreq=${OPTARG};;
        f) restartFile=${OPTARG};;
        p) rep=${OPTARG};;
    esac
done
```
We now list the inputs, flags, description, and expected values for all inputs to the `run.sh`.

| Input        |Flag | Expected values         | Description                                |
|--------------|-----|-------------------------|--------------------------------------------|
| gas          | -g  | [water, ch4, co2]       | is there ch4 or co2 gas in the system?     |
| water        | -w  | [tip4pice, tip4p2005]   | water potential used                       |
| seed         | -d  | [123456, etc.]          | random number generator seed               |
| ensemble     | -e  | [npt, nvt, nve]         | ensemble fix                               |
| steps        | -s  | [1000000, etc.]         | timesteps to run simulation                |
| restart      | -r  | [false, restart]        | is this a restart simulation?              |
| restartFreq  | -n  | [250000, etc.]          | freq to output restart file in timesteps   |
| restartFile  | -f  | ["./dir/FileName.data"] | name of restart file                       |
| rep          | -p  | [false, true]           | is this a replicate simulation?            |

*NOTE:* for our purposes, the following inputs will be always set to their default values.
| Input        |Flag | Default for our purposes      |
|--------------|-----|-------------------------------|
| gas          | -g  | [water]                       |
| water        | -w  | [tip4pice]                    |
| seed         | -d  | [123456]                      |
| ensemble     | -e  | [npt]                         |
| steps        | -s  | [20000]                       |
| restart      | -r  | [false]                       |
| restartFreq  | -n  | [2500]                        |
| restartFile  | -f  | ["./restart/restartEnd.data"] |
| rep          | -p  | [false]                       |

This means that the script call to execute our workflow using `run.sh` looks like:
```
./run.sh -g water -w tip4pice -d 123456 -e npt -s 20000 -r false -n 2500 -f "./restart/restartEnd.data" -p false
```
If this does not run, make sure your `run.sh` is executable with:
```
chmod +x run.sh
```
*NOTE:* Remember to `cd` into the main directory of this repository before running the `run.sh` script.

### B. PACKMOL
The first execution in our workflow is PACKMOL in section "1. PACKMOL" of the `run.sh` script. 

Objective: The goal of PACKMOL is to randomly populate the simbox of our system with the molecules we want and each of them within a tolerance distance from each other. This helps to minimize the repulsive forces between molecules that can arise when using a lattice command to populate the simbox with a regular pattern space between molecules.

In this section of the code, the script verifies if this is a restarted simulation or not using the variable `restart` which was initiated from our input flag `-r`. Since our simulation is not a restart (`$restart = 'false'`), flow moves into the `if` statement block and executes PACKMOL. Notice that we `cd` in and out of the PACKMOL directory to make sure we can access the input files and that the outputs are contained in the right directory.
```
# 1. Execute PACKMOL
echo -------------------------------------------------
echo ">>>PACKMOL: $(date)<<<"
echo ">>>"
if [ $restart = 'false' ]
then
    # cd into packmol dir
    cd 1_packmol/
    # run packmol
    packmol < in_2pack_H2O.inp
    # return to original project directory
    cd ../
fi
```

### C. Moltemplate
The second step in our workflow is to execute Moltemplate.

Objective: The goal of Moltemplate is to assign the force field potential parameters to each molecules in our system.

The first thing to do is to activate our lammps conda environment (which will also be necessary for the next step running LAMMPS). Next, the script once again checks if this is a restarted simulation (see previous section explanation). Since it is not, flow moves into the `if` statement and Moltemplate is executed. Notice that we `cd` into the directory `input_files` before executing Moltemplate. This is where we have most the required input files.

A couple of peculiarities in this block:
1. we need to rename our `system_H2O.lt` file to `system.lt`. This is so that we can properly run the script `cleanup_moltemplate.sh` after running Moltemplate to clean auxilliary files produced. This script expects the output files to be named on a basis of `system.lt`.
2. We use the output from PACKMOL by refering to it using relative reference `../../1_packmol/out_packed.xyz`
3. We move the output files up one directory level to separate them from the input files
```
# source lammps env
conda activate lammps

# 2. Execute MOLTEMPLATE
echo -------------------------------------------------
echo ">>>MOLTEMPLATE: $(date)<<<"
echo ">>>"
if [ $restart = 'false' ]
then
    # go into 2_moltemplate/input_files/
    cd 2_moltemplate/input_files
    # define system identify
    cp system_H2O.lt system.lt
    # run moltemplate giving system.lt as input and the data file generated by packmol
    moltemplate.sh -xyz ../../1_packmol/out_packed.xyz system.lt
    # delete input lammps script produced
    rm system.in
    # remove extra information in settings and data files
    cleanup_moltemplate.sh
    # copy the output files to the lammps directory for execution
    cp system.data ../../3_lammps
    # cp system.data system.in* ../../3_lammps/equil
    # move the output files up to the parent directory
    mv -f system.data system.in* ../
    # delete the temporary files produced by moltemplate
    # these are good for debugging, comment if debugging is necessary
    rm -rf output_ttree/
    # return to original project directory
    cd ../../
fi
```

### D. LAMMPS
In the final step of our workflow we execute the actual molecular dynamics simulation using LAMMPS.

Objective: In this step, we conduct the molecular dynamics algorithm and produce molecular trajectory file (and possibly other output files as desired).

In this block, we execute LAMMPS with `lmp_serial` (you may change this to your needs). We feed the LAMMPS input script using the flag `-in`, and we specify the name of the output log file with the flag `-log`. Additionally, we feed the input parameters to our `run.sh` script as parameters (see previous sections) into our LAMMPS input script. Finally, our simulation will output various restart files (at every `$restartFreq` timesteps) and one final restart file at the end of the simulation, which is named different than the periodic ones. So, at the end of this block, we delete all the intermediate restart files, IFF the final restart file was produced. This is to save space. If the simulation was trunkated, and no final restart file was produce, then the intermediate files are kept.
```
# 3. Execute LAMMPS
echo -------------------------------------------------
echo ">>>LAMMPS: $(date)<<<"
echo ">>>"

cd 3_lammps/
# run 
lmp_serial -in run.ag.in -log log.$ensemble.lammps \
-v gas $gas \
-v wat $water \
-v seed $seed \
-v restart_sim $restart \
-v restart_freq $restartFreq \
-v ensemble $ensemble \
-v restartFile $restartFile \
-v steps $steps \
-v replicate $rep
# remove extra restarts if the restartEnd exists
if [ -f restart/restartEnd_${ensemble}_${steps}step.data ]
then
    rm restart/restart_$ensemble*
fi
# back up to test dir
cd ../
```

## 4. Analyses
The output from the LAMMPS molecular dynamics simulation will be saved in the file `log.npt.lammps` (this is how we named it according to the `-log` flag input to `run.sh`) in the directory `3_lammps/`. There are various ways to quantify the performance of your MD simulation, and the one we performed here runs for 40,000 fs (0.04 ns) which is very short. Here we will talk about some basic indicators that should be build upon for a definitive analysis of the simulation output. In the following sections, we will be discussing sections of the log file output by simulation.

### A. Initialization section
The initialization of our system was completed without errors. We can see in the log file that the correct `pair_style` was used, and the right data files were read in. The simulation box (simbox) was created and populated using the coordinates coming back from our PACKMOL output. I also like to check that the number of atoms, bonds, angles, and molecules are all as expected.

### B. Minimization section
Our minimization has ran and completed without errors. This energy minimum was achieved within the tolerance limits provided. Once the minimization is complete, we can assign the temperature to the system. The temperature assignment in LAMMPS MD is acomplished by setting velocities to each atom. This is done via the `velocity` command using the `seed` value we fed as input and a Gaussian distribution. 

### C. Thermo table
1. Temperature and pressure
    This simulation is ran with the isothermal-isobaric (NPT) ensemble, so both a thermostat and barostat are applied to the system. We would like to see the results reflect the setpoints. It is clear that we quickly begin to trend towards the temperature setpoint. However, the pressure trend is not so clear, especially if looking at the raw pressure value. The barostat controls pressure by adjusting the volume of the simbox. This introduces considerable mechanical disturbances to the system. The result is a relatively unstable raw pressure value, which should be transformed to a moving average so that we can eliminate some of the signal noise and see the long term trend. Note that even in the moving average column, the value still fluctuates. This is normal due to the barostat. Ideally, we would continue this simulation longer, and then restart it in a canonical (NVT) ensemble to remove the barostat and test that the pressure remains stable without mechanical control.
2. Density
    The density is a good indicator if we know what the expected value of the system should be. In our case, we have a simulation of pure water at low temperature. We still expect the density to be around 0.99 or 0.98 g/cm3. This is the value that we are trending to in this simulation.
3. Potential energy
    We want to make sure that the potential energy of the system approaches (and eventually stabilizes to) the value provided by the original parametrization of the water force field potential. In the case of our simulation, we have used the TIP4P/Ice potential, which has a literature value of [-14.60 kcal/mol](https://aip.scitation.org/doi/full/10.1063/1.1931662). We can see that our simulation is beginning to trend towards that value (-13.35 kcal/mol). This indicates that we are not yet well equilibrated, and we should continue the simulation for a longer simulation time. 
4. Performance table
    On line 339, LAMMPS outputs the performance associated with the thermo table above it. This line in particular has good information on how the simulation performed on the hardware you ran it on. The speed of execution can be inferred from: `6.306 ns/day, 3.806 hours/ns, 36.491 timesteps/s`. This can be used to adjust how you assign hardware to your system. This is out of our scope here, but it's interesting to keep in mind. 
