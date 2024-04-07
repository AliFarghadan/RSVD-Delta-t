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

To use $\text{RSVD}-\Delta t$, you can follow these steps:
 
1. Install the required dependencies:
	+ OpenMPI (or similar MPI package)
	+ C++ compiler
	+ Petsc and Slepc packages (working versions are Petsc 3.18.5, Slepc 3.18)

2. Download and extract the $\text{RSVD}-\Delta t$ package.
3. Navigate to the directory containing the extracted files.
4. Run the main script (RSVDT.sh) with appropriate input parameters to compute resolvent modes.
5. Obtain the output data.

### Install Petsc and Slepc

You can follow the instructions from the official websites:

- [Petsc](https://petsc.org/release/install)
- [Slepc](https://slepc.upv.es/documentation)

A suggested configuration of Petsc is as follows:\
./configure --with-debugging=0 --with-scalar-type=complex --with-64-bit-indices\
This configuration ensures computations with complex values, enables the use of 64-bit integer numbers, and enhances speed by disabling debugging.

Note: Newer versions of Petsc and Slepc do not change the code principles; however, syntaxes might need to be updated.

## Additional resources

### References

* [Scalable resolvent analysis for three-dimensional flows](https://arxiv.org/pdf/2309.04617.pdf), *JCP (under review)*, 2024
* [A randomized time-domain algorithm for efficiently computing resolvent modes](https://arc.aiaa.org/doi/10.2514/6.2021-2896), *AIAA AVIATION*, 2021

### Contact information

Ali Farghadan, University of Michigan, aliii@umich.edu\
Aaron Towne, University of Michigan, towne@umich.edu
