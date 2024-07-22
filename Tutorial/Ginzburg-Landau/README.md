# Tutorial: Computing resolvent modes of the Ginzburg-Landau system using $\text{RSVD}-\Delta t$

---

## Introduction

In this tutorial, we will walk through the process of using the $\text{RSVD}-\Delta t$ algorithm on a Ginzburg-Landau test case. We showcase both the transient run (the pre-processing step) and the main algorithm (for computing resolvent modes). We also describe all related input variables, and the logisitics of I/O directories. The linearized operator is provided in [here](./). This operator is the operator we used in our [paper](https://arxiv.org/pdf/2309.04617.pdf).

## Prerequisite modules

Make sure you have the following prerequisite modules installed:

- GCC
- OpenMPI
- PETSc and SLEPc libraries

Note that the same versions of GCC and OpenMPI are used to compile PETSc and SLEPc that are used to build the executable. Please refer to the [README](/README.md) file for instructions on installing PETSc and SLEPc.

## List of input variables

Here is a list of variables used for this test case. Please note that you will need to modify these variables to suit your own projects.
Important: For boolean flags, the following values are equivalent:
* False, false, and 0 all represent the same value
* True, true, and 1 all represent the same value


```yaml
# Root directory (string) 
RootDir:         /path/to/root/directory/

# Results directory (string)
# The resolvent modes/gains will be saved in RootDir + ResultsDir + ResolventModes_i
ResultsDir:      /path/to/results/

# The linearized operator (matrix - usually very sparse and saved in binary format)
# The operator directory is defined as RootDir + OperatorDir
OperatorDir:     matrices/A_GL

# Number of test vectors (integer)
k:               2

# Number of power iterations (integer)
q:               0

# Base frequency (real)
# The resolvent modes will be computed at freuqency range -Nw/2*w:w:Nw/2*w
w:               0.05

# Number of frequencies to resolve (integer - even number)
Nw:              42

# Convert frequencies to angular frequencies, w_angular = 2*pi*w (boolean)
TwoPI:           false

# Transient length (real)
TransLen:        10

# Time step (real)
dt:              0.003

# Transient removal flag (boolean)
TransRemoval:    true

# Display option (integer 0 <= Display <= 2)
Display:         2

# Discounting flag (boolean)
# Applies discounting for unstable linear systems, A_discounted = A - beta I
Discounting:     false

# beta value for discounting if true (real > 0)
beta:            0.1

# Seeding random number to replicate data if needed (integer)
Rseed:           14

# Runs the transient simulation if true, otherwise RSVD-delta-t will run (boolean)
# All variables below this are relevant when TransRun is true
TransRun:        false

# Number of periods to integrate (integer)
Periods:         3

# Saves the transient outputs if true (boolean)
# The snapshots will be saved in RootDir + ResultsDir + TransientSnapshots_i
TransSave:       false

# Saves the snapshots every mod number (integer)
mod:             500

# Estimates the transient error using our strategy if true (boolean)
TransRemovalEst: false

# Divergence value, the transient simulaton will stop when reached it (real)
DivVal:          1e3

# Convergence value, the transient simulaton will stop when reached it (real)
ConVal:          1e-6
```

## Running the Algorithm
The process of running the $\text{RSVD}-\Delta t$ algorithm can be divided into two main parts:

1. ### Transient Simulation

Before running the actual $\text{RSVD}-\Delta t$ algorithm, it's often necessary to perform a transient simulation to ensure the stability of the system and visualize the transient decay rate.

To run the transient simulation, set the `TransRun` flag to `true` in the input file (`variables.yaml`).

#### Simulation variables

* `Periods` determines the length of transient simulation.
* We define the period length as $T = 2\pi/\omega_{min}$, where $\omega_{min}$ (defined as `w`) is the base frequency.
* `TransRemovalEst`, if `true`, applies the transient removal strategy we developed for slowly decaying systems. It estimates the updated transient residual at the end of each period.
* `DivVal` and `ConVal` are divergence and convergence values, respectively, which stops the simulation if the transient norm reaches either of them.
* `Rseed` is the seeding number.

#### Saving results

* We create a folder in the results directory with a fixed prename "TransientSnapshots_&lt;int&gt;", where &lt;int&gt; is an integer starting from 0. If "TransientSnapshots_i" exists, the code increments the integer until the folder name is unique, ensuring results from different simulations are not overwritten.
* `TransSave`, if `true`, saves the snapshots "q_transient_&lt;int&gt;" every `mod` time step in the results directory. It also saves the norm of snapshots in a vector "q_transient_norms" and the last snapshot "q_transient_last_snapshot".
* In case the `TransSave` is `false`, we only save "q_transient_norms" and "q_transient_last_snapshot".
* In case the `TransRemovalEst` is `true`, the simulation saves initial and updated transient norms to "Initial_transient_norm_period_&lt;int&gt;" and "Updated_transient_norm_period_&lt;int&gt;", respectively, at the end of each period. For instance, "Initial_transient_norm_period_1" and "Updated_transient_norm_period_1" are the norm of snapshots across the frequency range at the end of the first period.

### Note: all variables here have a default value in case you comment out variables, you will get a warning that the default values are used for those variables that are missed in your input list. Here are the default values:

#### Default Values

* `TransRun`: `false`
* `Periods`: `1`
* `TransSave`: `false`
* `TransRemovalEst`: `false`
* `DivVal`: `1e-6`
* `ConVal`: `1e3`
* `Rseed`: `1373`

2. ### $\text{RSVD}-\Delta t$ Algorithm

The executable runs the $\text{RSVD}-\Delta t$ algorithm by default unless `TransRun` is `true`. 

#### Simulation variables

* `RootDir` determines the root directory path for the simulation.
* `ResultsDir` specifies the results directory path where output files will be saved. The `ResultsDir` is the folder directory that needs to exist in the `RootDir`. If it does not exist, we make directory with the path provided. 
* `OperatorDir` defines the directory path for the linearized operator matrix. If the operator is in the `RootDir`, just provide the operator name (*e.g.* A_GL). Otherwise, you need to specify the path to the operator relative to `RootDir` (*e.g.* matrices/A_GL)
* `k` represents the number of test vectors.
* `q` specifies the number of power iterations performed.
* `w` is the base frequency.
* `Nw` determines the number of frequencies to resolve.
* `TwoPI` is a boolean flag that converts frequencies to angular frequencies by multiplying with 2*pi.
* `TransLen` sets the transient length. We usually get an estimation of the `TransLen` from the transient simulation.
* `dt` specifies the time step.
* `TransRemoval` is a boolean flag that toggles the transient removal strategy on or off. When set to `true`, this flag enables the transient removal strategy to be performed at the end of each time step, refining the steady-state snapshots for improved accuracy.
* `Discounting` is a boolean flag that applies discounting for unstable linear systems. `beta` is the beta value used for discounting when `Discounting` is `true`. You must specify a positive real value for `beta`; otherwise, the algorithm will exit with an error message.
* `Rseed` is the seeding random number used to replicate data.
* `Display` is an integer option that controls the level of output display, with three possible values: 0, 1, and 2.
	+ `Display = 0`: Minimal output, with no information displayed.
	+ `Display = 1`: Standard output, displaying:
		- Problem information
		- Elapsed time for each test vector
		- Total elapsed time
		- Estimated remaining time
	+ `Display = 2`: Detailed output, including everything from `Display = 1`, plus:
		- Progress percentage of the first test vector

#### Saving results

* We create a folder in the results directory with a fixed prename "ResolventModes_&lt;int&gt;", where &lt;int&gt; is an integer starting from 0. If "ResolventModes_i" exists, the code increments the integer until the folder name is unique, ensuring results from different simulations are not overwritten.
* Once the computation is complete, `k` response modes (each of size `N` × `Nw`) are saved as "U_hat_k&lt;int&gt;", where &lt;int&gt; represents the &lt;int&gt;th mode. Similarly, forcing modes are saved as "V_hat_k&lt;int&gt;". For instance, "U_hat_k1" and "V_hat_k1" contain the optimal response and forcing modes, respectively, across all frequencies of interest. Finally, gains are saved as a single matrix "S_hat" of size `k` × `Nw`.

Note: not all variables have default values. If a variable is not specified, you will receive a warning or error message.

#### Default values are provided for the following variables:

* `k`: `3`
* `q`: `0`
* `TwoPI`: `false`
* `TransLen`: `100`
* `TransRemoval`: `false`
* `Display`: `2`
* `Discounting`: `false`
* `Rseed`: `1373`

#### The following variables must be specified with no default values:

* `RootDir`
* `ResultsDir`
* `OperatorDir`
* `w`
* `Nw`
* `dt`
* `beta` (only when `Discounting` is `true`)

## Transfer Data Between MATLAB and PETSc/SLEPc Binary Format

You may frequently need to transfer data between PETSc/SLEPc and MATLAB for post-processing results. This section will guide you through the process of converting data formats to ensure compatibility and seamless integration between MATLAB and PETSc/SLEPc.

### Converting Data from MATLAB to PETSc/SLEPc

MATLAB can save data in various formats, but for use with PETSc/SLEPc, we need to ensure the data is saved in a binary format compatible with these libraries. PETSc provides MATLAB functions to facilitate this conversion.

1. **Add PETSc MATLAB Path:**

    Ensure you add the PETSc MATLAB interface to your MATLAB path:
    ```matlab
    addpath('/path/to/PETSc/share/petsc/matlab/');
    ```

2. **Save Data in MATLAB:**

    Create your matrix in MATLAB and save it using PETSc's binary write function. We have provided A_GL in both binary and .mat formats which you can test.
    ```matlab
    PetscBinaryWrite('/path/to/your/matrix/A', A_GL, 'complex', true, 'indices', 'int64');
    ```
    Note that `A` is an example name for the binary saved file. `A_GL` is the variable in MATLAB. Depending on the PETSc architecture that you have compiled, `'complex', true` and `'indices', 'int64'` can be different. Please refer to `PetscBinaryWrite` function for more information.

3. **Load Data in PETSc:**

    Here is a general way in your PETSc/SLEPc environment, you can now load the binary file:
    ```c
    Mat A;
    PetscViewer viewer;
    
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "/path/to/your/matrix/A", FILE_MODE_READ, &viewer);
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetType(A, MATSEQAIJ); // or MATMPIAIJ if parallel (depending on your matrix, you can vary the type)
    MatLoad(A, viewer);
    PetscViewerDestroy(&viewer);
    ```

### Converting Data from PETSc/SLEPc to MATLAB

To transfer data from PETSc/SLEPc to MATLAB, follow these steps:

1. **Save Data in PETSc/SLEPc:**

    Save the matrix data in PETSc's binary format:
    ```c
    PetscViewer viewer;
    
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "/path/to/your/matrix/A", FILE_MODE_WRITE, &viewer);
    MatView(A, viewer);
    PetscViewerDestroy(&viewer);
    ```

2. **Load Data in MATLAB:**

    Read the PETSc binary file in MATLAB:
    ```matlab
    addpath('/path/to/PETSc/share/petsc/matlab/');
    A = PetscBinaryRead('/path/to/your/matrix/A', 'complex', true, 'indices', 'int64');
    ```

## Conclusion
In this tutorial, we've outlined the basic steps for using the $\text{RSVD}-\Delta t$ algorithm with a Ginzburg-Landau test case and transferring data between MATLAB and PETSc/SLEPc. By following these instructions and adapting them to your specific problem, you can harness the power of this algorithm for your scientific or engineering simulations.

Feel free to adjust the variables and configurations as needed for your particular use case. Experimentation and iteration are often necessary for achieving optimal results.

If you have any further questions or encounter any issues, don't hesitate to reach out for assistance.





