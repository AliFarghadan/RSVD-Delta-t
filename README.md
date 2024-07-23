# $\text{RSVD}-\Delta t$: Randomized Singular Value Decomposition with Time-stepping for large-scale resolvent analysis

Copyright © 2024 The Regents of the University of Michigan

Welcome to $\text{RSVD}-\Delta t$, a novel algorithm designed to address the computational challenges associated with studying coherent structures in large-scale flows using resolvent analysis. This README file provides an overview of the algorithm, its features, and instructions for usage. Below, behold a visual representation showcasing the resolvent response modes of a three-dimensional jet, exemplifying the utility of $\text{RSVD}-\Delta t$.

![ResModes](image.png)

## Overview

Resolvent analysis is a valuable tool for studying coherent structures in turbulent flows. However, its application to inherently three-dimensional flows and large systems has been limited by the computational cost of computing resolvent modes. $\text{RSVD}-\Delta t$ presents a solution to these challenges by combining randomized singular value decomposition (RSVD) with an optimized time-stepping method, resulting in significant reductions in CPU cost and memory requirements.

## Features

* CPU time scales linearly with problem dimension: $\text{RSVD}-\Delta t$ offers linear scalability with the number of discrete degrees of freedom $O(N)$, significantly reducing computational overhead for large systems compared to existing algorithms.
* Memory usage scales linearly with problem dimension: $\text{RSVD}-\Delta t$ minimizes memory usage via streaming Fourier sums, enabling the study of high-dimensional systems.
* Error control strategies: $\text{RSVD}-\Delta t$ incorporates strategies such as power iteration and transient removal to control errors and ensure the reliability of computed resolvent modes.

## Usage

To use $\text{RSVD}-\Delta t$, follow these steps:

1. **Install dependencies**:
	* OpenMPI (or similar MPI package)
	* C++ compiler
	* PETSc and SLEPc packages (our codes are developed on versions PETSc 3.19.4 and SLEPc 3.19)
2. **Download and extract** the $\text{RSVD}-\Delta t$ package.
3. **Navigate** to the directory containing the extracted files.
4. **Build the executable** using the makefile in the package.
5. **Run the executable** with appropriate input parameters to compute resolvent modes.

We've shown an example usage of our code in the [Tutorial](./Tutorial/Ginzburg-Landau).

## Install PETSc and SLEPc

You can follow the instructions from the official websites:

- [PETSc](https://petsc.org/release/install)
- [SLEPc](https://slepc.upv.es/documentation/instal.htm)

A suggested configuration of PETSc is as follows:\
`./configure --with-debugging=0 --with-scalar-type=complex --with-64-bit-indices PETSC_ARCH=complex-opt`\
This configuration ensures computations with complex values, enables the use of 64-bit integer numbers, and enhances speed by disabling debugging. Here are a few notes:

1. Newer versions of PETSc and SLEPc do not change the core principles of the code; however, syntax updates might be required.
2. `PETSC_ARCH=<PETSc-arch-name>` can be chosen differently.
3. If you do not have MPICH installed locally, you can add `--download-mpich` to the configuration options. Please read the [PETSc configuration](https://petsc.org/main/install/install/) guide for more information.
4. If you need additional packages such as MUMPS, you can add `--download-mumps`. In that case, you need to update the makefile accordingly.

## Setting up environment variables for PETSc and SLEPc

After you have completed the installation of PETSc and SLEPc, you need to define the environment variables `PETSC_DIR` and `SLEPC_DIR` to point to the installation directories of these packages. This is necessary for the proper functioning of your applications that rely on PETSc and SLEPc.

### Steps to Define Environment Variables

1. **Identify the Installation Paths**:
   - Determine the directories where PETSc and SLEPc are installed. For example, if PETSc is installed in `/path/to/PETSC` and SLEPc is installed in `/path/to/SLEPC`, you will use these paths in the following steps.

2. **Export the Environment Variables**:
   - Open a terminal and export the environment variables by running the following commands:
     ```sh
     export PETSC_DIR=/path/to/PETSC
     export SLEPC_DIR=/path/to/SLEPC
     ```

3. **Persist the Environment Variables**:
   - To avoid setting these variables every time you open a new terminal, add the export commands to your shell configuration file. If you are using a Linux environment with Bash, you can add these lines to your `~/.bashrc` file:
     ```sh
     echo 'export PETSC_DIR=/path/to/PETSC' >> ~/.bashrc
     echo 'export SLEPC_DIR=/path/to/SLEPC' >> ~/.bashrc
     ```
   - After adding these lines, apply the changes by running:
     ```sh
     source ~/.bashrc
     ```

   - For other environments (e.g., Windows or macOS), you can add the equivalent commands to your respective shell configuration files (e.g., `.bash_profile`, `.zshrc`).

By setting these environment variables, you ensure that the paths to PETSc and SLEPc are correctly defined, allowing you to compile and run your applications without needing to manually set the paths each time.


## Makefile usage

The `makefile` is written to build the executable from the source code. You can use the following commands:

* `make` to build the executable
* `make clean` to remove the executable

The executable name will be `RSVDt`.

Note: When making the executable via `make` command, by default, `PETSC_ARCH=complex-opt` is assumed. If your PETSC_ARCH name is different, you must specify it in your make command: `make PETSC_ARCH=<PETSc-arch-name>`.

## Example Jobfile

To illustrate the usage of the $\text{RSVD}-\Delta t$ algorithm, we provide an example jobfile. This jobfile specifies the configuration for running the algorithm on a computing cluster.

```bash
#!/bin/bash
#SBATCH --job-name=<name>
#SBATCH --nodes=<count>
#SBATCH --ntasks-per-node=<count>
#SBATCH --mem=<memory>
#SBATCH --time=<dd-hh:mm:ss>

module load gcc/<version> openmpi/<version>

cd /path/to/executable/
make PETSC_ARCH=complex-opt # this will compile the source files to create the executable, or do nothing if the executable is already up-to-date

mpiexec RSVDt -inputs variables.yaml  # runs the algorithm for a given set of variables in variables.yaml file
# Note: You may need to use "srun" or "mpirun" instead of "mpiexec" depending on your installation and/or cluster configurations.
```

In this jobfile, `RSVDt` is the executable for our algorithm, and `variables.yaml` contains the prerequisite variables required for running the algorithm.

## Additional resources

### References

* [Scalable resolvent analysis for three-dimensional flows](https://arxiv.org/pdf/2309.04617.pdf), *JCP (under review)*, 2024
* [A randomized time-domain algorithm for efficiently computing resolvent modes](https://arc.aiaa.org/doi/10.2514/6.2021-2896), *AIAA AVIATION*, 2021

### Contact information

Ali Farghadan, University of Michigan, aliii@umich.edu\
Aaron Towne, University of Michigan, towne@umich.edu

### Cite as

```cite
@Article{Farghadan2023scalable,
  title={Scalable resolvent analysis for three-dimensional flows},
  author={Farghadan, A. and Martini, E. and Towne, A.},
  journal={arXiv preprint arXiv:2309.04617},
  year={2023}
}
```
