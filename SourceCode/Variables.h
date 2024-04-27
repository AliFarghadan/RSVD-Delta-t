
#ifndef VARIABLES_H
#define VARIABLES_H

typedef struct {
	PetscInt        N;                                  /* problem size (state dimension) */
	PetscInt        k;                                  /* number of test vectors */
	PetscInt        q;                                  /* number of power iterations */
	PetscInt        Nw;                                 /* number of input/output frequencies to resolve */
	PetscInt        Nw_eff;                             /* effective number of frequencies (becomes half for real-valued matrices) */
	PetscReal       w;                                  /* base frequency */
	PetscBool       twopi_flg;                          /* base frequency multiplies by 2*pi if true */
} RSVD_vars;

typedef struct {
	PetscInt        Ns;                                 /* number of time steps within one steady-state period */
	PetscInt        Nt;                                 /* number of time steps for transient */
	PetscInt        ResRatio;                           /* resolution ratio between Nt and Nw (which is an integer) */
	PetscReal       TransLen;                           /* transient length */
	PetscReal       dt;                                 /* time step */
	PetscBool       flg_real_A;                         /* real-valued A if true */
	PetscBool       TransRemoval;                       /* remove the transient projection onto on-fly subspace if true */
	PetscBool       flg_dir_adj;                        /* direct run if true, otherwise adjoint run */
} TS_vars;

typedef struct {
	PetscBool       flg_disc;                           /* discounting flag */
	PetscReal       Beta;                               /* discounting parameter */
} Discounting;

typedef struct {
	RSVD_vars       RSVD;                               /* RSVD variables */
	TS_vars         TS;                                 /* time-stepping variables */
	Discounting     disc;                               /* discounting variables */
	PetscInt        display;                            /* 0: None, 1: Partial, 2: Full progress display */
} RSVDt_vars;

typedef struct {
	PetscBool       TransientRun;                       /* runs the transient simulation if true */
	PetscBool       TransientSave;                      /* saves the transient outputs if true */
	// PetscBool       TransientRunDir;                   /* runs the forward system if true, otherwise adjoint */
	PetscBool       EffTransTest;                       /* estimates the transient error using our strategy if true */
	PetscInt        Rseed;                              /* seeding random number to replicate data if needed */
	PetscInt        mod;                                /* saves the snapshots every mod number */
	PetscInt        Periods;                            /* number of periods to integrate */
	PetscReal       DivVal;                             /* divergence value, stops the simulation when reached it */
	PetscReal       ConVal;                             /* convergence value, stops the simulation when reached it */
} TransientRun_vars;

typedef struct {
	Mat             idft;                               /* inverse discrete Fourier transform matrix */
	Mat             dft;                                /* discrete Fourier transform matrix */           
} DFT_matrices;

typedef struct {
	Mat             A;                                  /* LNS matrix */
	RSVDt_vars      RSVDt;                              /* specifically beta value if discounted */
} LNS_vars;

typedef struct {
	Vec             F1,F2,F3;                           /* RK4 input forcing */
	Vec             k1,k2,k3,k4,y_temp;                 /* RK4 temporary vectors */
} TS_matrices;

typedef struct {
	Mat             W_in,W_out;                         /* weight matrices */
	PetscBool       flg_Win,flg_Wout;                   /* weight flags */
	PetscBool       flg_B,flg_C;                        /* spatial flags */
	PetscBool       flg_diag;                           /* true if weight matrices are diagonal */
	Vec             W_in_vec,W_out_vec;                 /* weight diagonal enteries */
	Vec             B_vec,C_vec;                        /* spatial diagonal enteries */
} WS_matrices;

typedef struct {
	Mat             F_hat,Y_hat,Q_hat;                  /* RSVD matrices */
} RSVD_matrices;

typedef struct {
	Mat             M_tilde_all,Y_hat_up;               /* transient removal matrices */
} TS_removal_matrices;

typedef struct {
	Mat             U_hat,U_hat_re;                     /* response resolvent modes */
	Mat             V_hat,V_hat_re;                     /* forcing resolvent modes */
	Mat             S_hat;                              /* resolvent gains */
	PetscBool       all_in_one_flg;                     /* save modes altogether if true, otherwise in separate files */
} Resolvent_matrices;

typedef struct {
	char            root_dir[PETSC_MAX_PATH_LEN];       /* root directory */
	char            input[PETSC_MAX_PATH_LEN];          /* input folder */
	char            output[PETSC_MAX_PATH_LEN];         /* output folder */
	char            Mat_filename[PETSC_MAX_PATH_LEN];   /* filename of the LNS operator */
	char            filename[PETSC_MAX_PATH_LEN];       /* filename of the forcing input if given */
	char            IO_dir[PETSC_MAX_PATH_LEN];         /* I/O directory */
	char            FolderDir[PETSC_MAX_PATH_LEN];      /* results folder directory */
} Directories;

#endif