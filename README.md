# Basic information of $\text{RSVD}-\Delta t$ for users

Copyright © 2024 The Regents of the University of Michigan

Welcome to $\text{RSVD}-\Delta t$, a novel algorithm designed to address the computational challenges associated with studying coherent structures in large-scale flows using resolvent analysis. This README file provides an overview of the algorithm, its features, and instructions for usage.

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

We've shown an example usage of our code in the [Tutorial](./Tutorial).

## Install PETSc and SLEPc

You can follow the instructions from the official websites:

- [PETSc](https://petsc.org/release/install)
- [SLEPc](https://slpec.upv.es/documentation)

A suggested configuration of PETSc is as follows:\
./configure --with-debugging=0 --with-scalar-type=complex --with-64-bit-indices\
This configuration ensures computations with complex values, enables the use of 64-bit integer numbers, and enhances speed by disabling debugging.

Note: Newer versions of PETSc and SLEPc do not change the code principles; however, syntaxes might need to be updated.

## Makefile usage

The `makefile` is written to build the executable from the source code. You can use the following commands:

* `make` to build the executable
* `make clean` to remove the executable

The executable name will be `RSVDt`.

Note: By default, `PETSC_ARCH=complex-opt`. If your PETSC_ARCH name is different, you must specify it in your make command: `make PETSC_ARCH=<PETSc-arch-name>`.

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

mpirun RSVD-delta-t -inputs variables.yaml  # runs the algorithm for a given set of variables in variables.yaml file
```

In this jobfile, `RSVD-delta-t` is the executable for our algorithm, and `variables.yaml` contains the prerequisite variables required for running the algorithm.

## Additional resources

### References

* [Scalable resolvent analysis for three-dimensional flows](https://arxiv.org/pdf/2309.04617.pdf), *JCP (under review)*, 2024
* [A randomized time-domain algorithm for efficiently computing resolvent modes](https://arc.aiaa.org/doi/10.2514/6.2021-2896), *AIAA AVIATION*, 2021

### Contact information

Ali Farghadan, University of Michigan, aliii@umich.edu\
Aaron Towne, University of Michigan, towne@umich.edu
