## Tutorial: Computing resolvent modes of the Ginzburg-Landau system using $\text{RSVD}-\Delta t$

---

### Introduction

In this tutorial, we will walk through the process of using the $\text{RSVD}-\Delta t$ algorithm on a Ginzburg-Landau test case. We showcase both the transient run (pre-processing step) and the main algorithm (for computing resolvent modes). We also describe all related input variables, and the logisitics of I/O directories. The linearized operator is provided in [folder](/./) 

###

### Prerequisite modules

Make sure you have the following prerequisite modules installed:

- GCC
- OpenMPI
- PETSc and SLEPc libraries

In fact, the same modules for GCC and OpenMPI are used to compile both Petsc and Slepc. Please refer to [README](/README.md)

### Example Jobfile

To illustrate the usage of the RSVD-delta-t algorithm, we provide an example jobfile. This jobfile specifies the configuration for running the algorithm on a computing cluster.

```bash
#!/bin/bash
#SBATCH --job-name=GL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40g
#SBATCH --time=01:00:00

module load gcc/10.3.0 openmpi/4.1.6

cd /path/to/root_directory/
make PETSC_ARCH=complex-opt

srun -c1 ./RSVD-delta-t -inputs variables.txt
```

In this jobfile, $\text{RSVD}-\Delta t$ is the executable for our algorithm, and variables.txt contains the prerequisite variables required for running the algorithm.

### List of input variables

Here, we provide a list of variables we used for this test case. You need to change the variables for your own project.



### Running the Algorithm
The process of running the RSVD-delta-t algorithm can be divided into two main parts:

Transient SimulationBefore running the actual RSVD-delta-t algorithm, it's often necessary to perform a transient simulation to initialize the system or prepare the data. To run the transient simulation, set the TransRun flag to true in the input file (e.g., variables.txt). This step is crucial for certain types of problems where initialization or preparation is required.
Executing RSVD-delta-t AlgorithmOnce the transient simulation is complete, or if it's not needed for your specific problem, set the TransRun flag to false in the input file. This signals the algorithm to execute the main RSVD-delta-t computation.

### Conclusion
In this tutorial, we've outlined the basic steps for using the RSVD-delta-t algorithm with a Ginzburg-Landau test case. By following these instructions and adapting them to your specific problem, you can harness the power of this algorithm for your scientific or engineering simulations.

Feel free to adjust the parameters and configurations as needed for your particular use case. Experimentation and iteration are often necessary for achieving optimal results.

If you have any further questions or encounter any issues, don't hesitate to reach out for assistance.





