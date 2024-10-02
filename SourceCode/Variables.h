
#ifndef VARIABLES_H
#define VARIABLES_H

typedef struct {
	PetscInt        N;                                /* problem size (state dimension) */
	PetscInt        k;                                /* number of test vectors */
	PetscInt        q;                                /* number of power iterations */
	PetscInt        Nw;                               /* number of input/output frequencies to resolve */
	PetscInt        Nw_eff;                           /* effective number of frequencies (Nw becomes half for real-valued matrices) */
	PetscReal       w;                                /* base frequency */
	PetscBool       TwoPI;                            /* base frequency multiplies by 2*pi if true */
	PetscInt        RandSeed;                         /* seeding random number to replicate data if needed */
} RSVD_vars;

typedef struct {
	PetscInt        Ns;                               /* number of time steps within one steady-state period */
	PetscInt        Nt;                               /* number of time steps for transient */
	PetscInt        ResRatio;                         /* resolution ratio between Nt and Nw (which is an integer) */
	PetscReal       TransientLength;                  /* transient length */
	PetscReal       dt;                               /* time step */
	PetscBool       RealOperator;                     /* real-valued A if true, otherwise complex-valued */
	PetscBool       TransientRemoval;                 /* remove the transient projection onto on-fly subspace if true */
	PetscBool       DirAdj;                           /* direct run if true, otherwise adjoint run */
} TS_vars;

typedef struct {
	PetscBool       DiscFlg;                          /* discounting flag */
	PetscReal       beta;                             /* discounting parameter */
} Discounting;

typedef struct {
	RSVD_vars       RSVD;                             /* RSVD variables */
	TS_vars         TS;                               /* time-stepping variables */
	Discounting     Disc;                             /* discounting variables */
	PetscInt        display;                          /* 0: None, 1: Partial, 2: Full progress display */
	PetscInt        SaveResultsOpt;                   /* options to save resolvent modes */
} RSVDt_vars;

typedef struct {
	PetscBool       TransRun;                         /* runs the transient simulation if true */
	PetscBool       TransSave;                        /* saves the transient outputs if true */
	PetscBool       TransRemovalEst;                  /* estimates the transient error using our strategy if true */
	PetscInt        TransSaveMod;                     /* saves the snapshots every mod number */
	PetscInt        TransPeriods;                     /* number of periods to integrate */
	PetscReal       TransDivVal;                      /* divergence value, stops the simulation when reached it */
	PetscReal       TransConVal;                      /* convergence value, stops the simulation when reached it */
} TransRun_vars;

typedef struct {
	Mat             idft;                             /* inverse discrete Fourier transform matrix */
	Mat             dft;                              /* discrete Fourier transform matrix */           
} DFT_matrices;

typedef struct {
	Mat             A;                                /* LNS matrix */
	RSVDt_vars      RSVDt;                            /* specifically beta value if discounted */
} LNS_vars;

typedef struct {
	Vec             F1,F2,F3;                         /* RK4 input forcing */
	Vec             k1,k2,k3,k4,y_temp;               /* RK4 temporary vectors */
} TS_matrices;

typedef struct {
	Mat             W_in,W_out;                       /* weight matrices */
	PetscBool       flg_Win,flg_Wout;                 /* weight flags */
	PetscBool       flg_B,flg_C;                      /* spatial flags */
	Vec             B_vec,C_vec;                      /* spatial diagonal enteries */
} WS_matrices;

typedef struct {
	Mat             F_hat,Y_hat,Q_hat;                /* RSVD matrices */
} RSVD_matrices;

typedef struct {
	Mat             M_tilde_all,Y_hat_up;             /* transient removal matrices */
} TS_removal_matrices;

typedef struct {
	Mat             U_hat;                            /* response resolvent modes */
	Mat             V_hat;                            /* forcing resolvent modes */
	Mat             S_hat;                            /* resolvent gains */
} Resolvent_matrices;

typedef struct {
	char            RootDir[PETSC_MAX_PATH_LEN];      /* root directory */
	char            ResultsDir[PETSC_MAX_PATH_LEN];   /* results folder */
	char            OperatorDir[PETSC_MAX_PATH_LEN];  /* LNS operator directory */
	char            filename[PETSC_MAX_PATH_LEN];     /* filename of the forcing input if given */
	char            IO_dir[PETSC_MAX_PATH_LEN];       /* I/O directory */
	char            FolderDir[PETSC_MAX_PATH_LEN];    /* results folder directory */
} Directories;

#endif

