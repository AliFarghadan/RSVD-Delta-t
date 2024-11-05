
#ifndef VARIABLES_H
#define VARIABLES_H

typedef struct {
	PetscInt      N;                                      /* problem size (state dimension) */
	PetscInt      Nb;                                     /* input size */
	PetscInt      Nc;                                     /* output size */
	PetscInt      k;                                      /* number of test vectors */
	PetscInt      q;                                      /* number of power iterations */
	PetscInt      Nw;                                     /* number of input/output frequencies to resolve */
	PetscInt      Nw_eff;                                 /* effective number of frequencies (Nw becomes half for real-valued matrices) */
	PetscReal     w;                                      /* base frequency */
	PetscBool     TwoPI;                                  /* base frequency multiplies by 2*pi if true */
	PetscBool     InputForcingFlg;                        /* input forcing from the specified directory if true, otherwise a random forcing */
	PetscInt      RandSeed;                               /* seeding random number to replicate data if needed */
} RSVD_vars;

typedef struct {
	PetscInt      Ns;                                     /* number of time steps within one steady-state period */
	PetscInt      Nt;                                     /* number of time steps for transient */
	PetscInt      ResRatio;                               /* resolution ratio between Nt and Nw (which is an integer) */
	PetscReal     TransientLength;                        /* transient length */
	PetscReal     dt;                                     /* time step */
	PetscBool     RealOperator;                           /* real-valued A if true, otherwise complex-valued */
	PetscBool     TransientRemoval;                       /* remove the transient projection onto on-fly subspace if true */
	PetscBool     DirAdj;                                 /* direct run if true, otherwise adjoint run */
} TS_vars;

typedef struct {
	PetscBool     DiscFlg;                                /* discounting flag */
	PetscReal     beta;                                   /* discounting parameter */
} Discounting;

typedef struct {
	RSVD_vars     RSVD;                                   /* RSVD variables */
	TS_vars       TS;                                     /* time-stepping variables */
	Discounting   Disc;                                   /* discounting variables */
	PetscInt      Display;                                /* 0: None, 1: Partial, 2: Full progress display */
	PetscInt      SaveResultsOpt;                         /* options to save resolvent modes */
} RSVDt_vars;

typedef struct {
	PetscBool     TransRun;                               /* runs the transient simulation if true */
	PetscBool     TransSave;                              /* saves the transient outputs if true */
	PetscBool     TransRemovalEst;                        /* estimates the transient error using our strategy if true */
	PetscBool     TransICFlg;                             /* initial vector from the specified directory if true, otherwise a random initial condition */
	PetscInt      TransSaveMod;                           /* saves the snapshots every mod number */
	PetscInt      TransPeriods;                           /* number of periods to integrate */
	PetscReal     TransDivVal;                            /* divergence value, stops the simulation when reached it */
	PetscReal     TransConVal;                            /* convergence value, stops the simulation when reached it */
} TransRun_vars;

typedef struct {
	Mat           idft;                                   /* inverse discrete Fourier transform matrix */
	Mat           dft;                                    /* discrete Fourier transform matrix */           
} DFT_matrices;

typedef struct {
	Mat           A;                                      /* LNS matrix */
	RSVDt_vars    RSVDt;                                  /* specifically beta value if discounted */
} LNS_vars;

typedef struct {
	Vec           F1,F2,F3;                               /* RK4 input forcing */
	Vec           k1,k2,k3,k4,y_temp;                     /* RK4 temporary vectors */
} TS_matrices;

typedef struct {
	Mat           W_q_sqrt;                               /* output weight matrix */
	Mat           W_q_sqrt_inv;                           /* inverse output weight matrix */
	Mat           W_f_sqrt_inv;                           /* inverse input weight matrix */
	Mat           B;                                      /* input matrix */
	Mat           C;                                      /* output matrix */
	PetscBool     InvInputWeightFlg;                      /* inverse input weight matrix from the specified directory if true, otherwise identity matrix */
	PetscBool     OutputWeightFlg;                        /* output weight matrix from the specified directory if true, otherwise identity matrix */
	PetscBool     InvOutputWeightFlg;                     /* inverse output weight matrix from the specified directory if true, otherwise identity matrix */
	PetscBool     InputMatrixFlg;                         /* input matrix from the specified directory if true, otherwise identity matrix */
	PetscBool     OutputMatrixFlg;                        /* output matrix from the specified directory if true, otherwise identity matrix */
} Weight_matrices;

typedef struct {
	Mat           F_hat,Y_hat;                            /* RSVD matrices */
} RSVD_matrices;

typedef struct {
	Mat           M_tilde_all;                            /* matrix containing all M_tilde for all frequencies */
	Mat           Y_hat_up;                               /* matrix of updated snapshots */
} TS_removal_matrices;

typedef struct {
	Mat           U_hat;                                  /* response resolvent modes */
	Mat           V_hat;                                  /* forcing resolvent modes */
	Mat           S_hat;                                  /* resolvent gains */
} Resolvent_matrices;

typedef struct {
	char          RootDir[PETSC_MAX_PATH_LEN];            /* root directory */
	char          ResultsDir[PETSC_MAX_PATH_LEN];         /* results folder */
	char          OperatorDir[PETSC_MAX_PATH_LEN];        /* LNS operator directory */
	char          IO_dir[PETSC_MAX_PATH_LEN];             /* I/O directory */
	char          FolderDir[PETSC_MAX_PATH_LEN];          /* results folder directory */
	char          InputForcingDir[PETSC_MAX_PATH_LEN];    /* input forcing directory */ 
	char          InvInputWeightDir[PETSC_MAX_PATH_LEN];  /* inverse input weight directory */ 
	char          OutputWeightDir[PETSC_MAX_PATH_LEN];    /* output weight directory */ 
	char          InvOutputWeightDir[PETSC_MAX_PATH_LEN]; /* inverse output weight directory */
	char          InputMatrixDir[PETSC_MAX_PATH_LEN];     /* input matrix directory */
	char          OutputMatrixDir[PETSC_MAX_PATH_LEN];    /* output matrix directory */
	char          TransICDir[PETSC_MAX_PATH_LEN];         /* initial vector directory for the transient simulation */ 
} Directories;

#endif

