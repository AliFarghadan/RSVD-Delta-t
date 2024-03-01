
static char help[] =   "V24: equivalent to V24. \
						First check: run and validate. \
						Second check: timing is not worse \
						Transient run for various number of vectors: validated \
						Forward action (k = 1): validated \
						QR decomposition (k = 1): validated \
						Adjoint action (k = 1): validated \
						SVD (k = 1): validated \
						Resolvent modes (k = 1): validated \
						Becoming LOVELY! \
						Resolvent modes (k > 1): validated \
						...? : Needs validation \n";

#include <slepcsvd.h>

typedef struct {
	PetscInt        N;                                /* problem size (state dimension) */
	PetscInt        k;                                /* number of test vectors */
	PetscInt        k_trans;                          /* number of test vectors for transient run */
	PetscInt        q;                                /* number of power iterations */
	PetscInt        Nw;                               /* number of input/output frequencies to resolve */
	PetscReal       w_min;                            /* base frequency */
} RSVD_vars;

typedef struct {
	PetscInt        Nt;                               /* number of time steps within one s.s. period */
	PetscInt        Nt_transient;                     /* number of time steps for transient */
	PetscInt        time_ratio;                       /* ratio between Nt and Nw (which is an integer) */
	PetscReal       t_transient;                      /* transient length */
	PetscReal       dt_max;                           /* maximum time step */
	PetscReal       dt;                               /* actual time step */
	PetscBool       flg_real_A;                       /* real-valued A if true */
	PetscBool       flg_eff_trans;                    /* remove the transient projection onto on-fly subspace if true */
	PetscBool       flg_dir_adj;                      /* direct run if true, else adjoint run */
} TS_vars;

typedef struct {
	PetscBool       flg_disc;                         /* discounting flag */
	PetscReal       beta;                             /* discounting parameter */
} Discounting;

typedef struct {
	PetscBool       TransientRun_flg;                 /* Run transient simulation flag */
	PetscBool       MatExpRun_flg;                    /* Run matrix exponential simulation flag */
	PetscBool       TSTimeDomainRun_flg;              /* Run time-stepping and collect time-domain snapshots flag */
} Function_calls;
// Note: No more than one function can be called. 
// If two or more functions are called simultaneously, the first one will be executed only.
// If none of the Function_calls flags are turned on, the RSVD-\Dt will be running.

typedef struct {
	RSVD_vars       RSVD;                             /* RSVD variables */
	TS_vars         TS;                               /* time-stepping variables */
	Discounting     disc;                             /* discounting variables */
} RSVDt_vars;

typedef struct {
	Mat             inv_dft;                          /* inverse discrete Fourier transform matrix */
	Mat             dft;                              /* discrete Fourier transform matrix */			
} DFT_matrices;

typedef struct {
	Mat             A;                                /* LNS matrix */
	RSVDt_vars      RSVDt;                            /* specifically beta value if discounted */
} LNS_params;

typedef struct {
	Vec             F1,F2,F3;                         /* RK4 input forcing */
	Vec             k1,k2,k3,k4,y_temp;               /* RK4 temporary vectors */
} TS_matrices;

typedef struct {
	Mat             W_in,W_out;                       /* weight matrices */
	PetscBool       flg_Win,flg_Wout;                 /* weight flags */
	PetscBool       flg_B,flg_C;                      /* spatial flags */
	PetscBool       flg_diag;                         /* true if weight matrices are diagonal */
	Vec             W_in_vec,W_out_vec;               /* weight diagonal enteries */
	Vec             B_vec,C_vec;                      /* spatial diagonal enteries */
} WS_matrices;

typedef struct {
	Mat             F_hat,Y_hat,Q_hat;                /* RSVD matrices */
} RSVD_matrices;

typedef struct {
	Mat             U_hat,U_hat_re;                   /* response resolvent modes */
	Mat             V_hat,V_hat_re;                   /* forcing resolvent modes */
	Mat             S_hat;                            /* resolvent gains */
} Resolvent_matrices;

typedef struct {
	char            root_dir[PETSC_MAX_PATH_LEN];     /* root directory */
	char            input[PETSC_MAX_PATH_LEN];        /* input folder */
	char            output[PETSC_MAX_PATH_LEN];       /* output folder */
	char            Mat_filename[PETSC_MAX_PATH_LEN]; /* filename of the LNS operator */
	char            filename[PETSC_MAX_PATH_LEN];     /* filename */
	char            file_dir[PETSC_MAX_PATH_LEN];     /* file directory */
} Directories;

static PetscErrorCode AdjointActionRK4(RSVD_matrices*, RSVDt_vars*, LNS_params*, \
						TS_matrices*, DFT_matrices*, WS_matrices*, Directories*);
static PetscErrorCode ApplyWSMats(RSVD_matrices*, RSVDt_vars*, WS_matrices*, PetscBool);
static PetscErrorCode CreateDFTiDFTMats(RSVDt_vars*, DFT_matrices*, LNS_params*);
static PetscErrorCode CreateRandomMat(RSVD_matrices*, RSVDt_vars*, Directories*);
static PetscErrorCode CreateForcingOnFly(Mat, DFT_matrices*, PetscInt, Vec);
static PetscErrorCode CreateTSMats(TS_matrices*, PetscInt);
static PetscErrorCode DestroyTSMats(TS_matrices*);
static PetscErrorCode DFT(Mat, DFT_matrices*, RSVD_matrices*, RSVDt_vars*);
static PetscErrorCode ForwardActionRK4(RSVD_matrices*, RSVDt_vars*, LNS_params*, \
						TS_matrices*, DFT_matrices*, WS_matrices*, Directories*);
static PetscErrorCode LNSStructure(LNS_params*, TS_matrices*, Directories*);
static PetscErrorCode MatExpActionRK4(Mat, Mat, PetscInt, PetscInt, LNS_params*, \
						WS_matrices*, TS_matrices*, RSVDt_vars*, Directories*);
static PetscErrorCode MatExpRunRK4(RSVDt_vars*, LNS_params*, WS_matrices*, TS_matrices*, Directories*);
// static PetscErrorCode MatExpWForcingRK4(RSVD_matrices*, DFT_matrices*, \
// 						LNS_params*, TS_matrices*, RSVDt_vars*, Directories*);
static PetscErrorCode PermuteMat(Mat, RSVDt_vars*);
static PetscErrorCode PowerIterationRK4(RSVD_matrices*, RSVDt_vars*, LNS_params*, \
						TS_matrices*, DFT_matrices*, WS_matrices*, Directories*);
static PetscErrorCode QR(RSVD_matrices*, RSVDt_vars*);
static PetscErrorCode ReadUserInput(RSVDt_vars*, WS_matrices*, LNS_params*, \
									Function_calls*, Directories*);
static PetscErrorCode ReadWSMats(WS_matrices*, Directories*);
static PetscErrorCode RecoverU(Mat, RSVDt_vars*, RSVD_matrices*, Resolvent_matrices*);
static PetscErrorCode ReducedSVD(RSVD_matrices*, RSVDt_vars*, WS_matrices*, Resolvent_matrices*);
static PetscErrorCode ReversePermuteMat(Mat, RSVDt_vars*);
static PetscErrorCode SaveOutputs(Resolvent_matrices*, Directories*);
static PetscErrorCode SaveSnapshots(Vec, PetscInt, RSVDt_vars*, Mat);
static PetscErrorCode SetupTimeFreqGrid(RSVDt_vars*);
static PetscErrorCode StoreQ(RSVD_matrices*, RSVDt_vars*);
static PetscErrorCode TransientRK4(Mat, Mat, PetscInt, PetscInt, PetscInt, LNS_params*, \
						WS_matrices*, TS_matrices*, RSVDt_vars*, PetscBool, Directories*);
static PetscErrorCode TransientRunRK4(RSVDt_vars*, LNS_params*, WS_matrices*, TS_matrices*, Directories*);
static PetscErrorCode TSActionRK4(RSVD_matrices*, DFT_matrices*, \
						LNS_params*, TS_matrices*, RSVDt_vars*, Directories*);
static PetscErrorCode TSRK4(Mat, DFT_matrices*, LNS_params*, \
					TS_matrices*, RSVDt_vars*, PetscInt, Vec);
static PetscErrorCode TSTimeDomainRunRK4(RSVD_matrices*, DFT_matrices*, LNS_params*, \
								WS_matrices*, TS_matrices*, RSVDt_vars*, Directories*);
static PetscErrorCode TSTransRK4(LNS_params*, TS_matrices*, RSVDt_vars*, PetscInt, Vec);

int main(int argc,char **args)
{
	PetscErrorCode         ierr;                      /* Petsc error code */
	Directories            dirs;                      /* I/O directories */
	Function_calls         Funcs;                     /* desired function to run */
	RSVDt_vars             RSVDt;                     /* RSVDt variables */
	LNS_params             LNS_mat;                   /* LNS matrices */
	DFT_matrices           DFT_mat;                   /* DFT and inverse DFT matrices */
	TS_matrices            TS_mat;                    /* time-stepping matrices */
	WS_matrices            WS_mat;                    /* weight and input/output matrices */
	RSVD_matrices          RSVD_mat;                  /* RSVD matrices */
	Resolvent_matrices     Res_mat;                   /* resolvent modes and gains */
	PetscViewer            fd;                        /* Petsc viewer */
	PetscLogDouble         t1, t2;                    /* time measurement variables */

	/*
		initializes Slepc/Petsc
	*/

	ierr = SlepcInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;

	/*
		Reads in user input parameters
	*/

	ierr = ReadUserInput(&RSVDt, &WS_mat, &LNS_mat, &Funcs, &dirs);

	LNS_mat.RSVDt.TS.flg_real_A ? PetscPrintf(PETSC_COMM_WORLD,"**** Real valued matrix in progress ****\n") : \
					PetscPrintf(PETSC_COMM_WORLD,"**** Complex valued matrix in progress ****\n");

	/*
		Reads in the LNS operator (+ discounting if desired)
	*/

	ierr = LNSStructure(&LNS_mat, &TS_mat, &dirs);CHKERRQ(ierr);
	ierr = MatGetSize(LNS_mat.A,&RSVDt.RSVD.N,NULL);CHKERRQ(ierr);

	/*
		Initializes the time-stepping variables
	*/

	ierr = SetupTimeFreqGrid(&RSVDt);CHKERRQ(ierr);

	/*
		Creates the discrete Fourier transform (DFT) and inverse DFT matrices
	*/

	ierr = CreateDFTiDFTMats(&RSVDt, &DFT_mat, &LNS_mat);CHKERRQ(ierr);

	/*
		Generates a random input matrix (or read in if desired)
	*/

	ierr = CreateRandomMat(&RSVD_mat, &RSVDt, &dirs);CHKERRQ(ierr);

	/*
		Read weight and spatial matrices (if applicable)
	*/

	ierr = ReadWSMats(&WS_mat, &dirs);CHKERRQ(ierr);

	/*
		Transient simulation (if desired -- run and exit)
	*/
	
	if (Funcs.TransientRun_flg) {
		ierr = TransientRunRK4(&RSVDt, &LNS_mat, &WS_mat, &TS_mat, &dirs);CHKERRQ(ierr);
		ierr = SlepcFinalize();
		return ierr;
	}

	/*
		Matrix exponential evolution (if desired -- run and exit)
	*/

	if (Funcs.MatExpRun_flg) {
		ierr = MatExpRunRK4(&RSVDt, &LNS_mat, &WS_mat, &TS_mat, &dirs);CHKERRQ(ierr);
		ierr = SlepcFinalize();
		return ierr;
	}

	/*
		Time-stepping action in the time domain (if desired -- run and exit)
	*/

	if (Funcs.TSTimeDomainRun_flg) {
		ierr = TSTimeDomainRunRK4(&RSVD_mat, &DFT_mat, &LNS_mat, &WS_mat, &TS_mat, &RSVDt, &dirs);CHKERRQ(ierr);
		ierr = SlepcFinalize();
		return ierr;
	}

	/************************************************************************* 
	 *************************************************************************
	    ****************  RSVD - $\Delta t$ algorithm  *******************
	    ****************      resolvent analysis       *******************
     *************************************************************************
	**************************************************************************/

	ierr = ForwardActionRK4(&RSVD_mat, &RSVDt, &LNS_mat, &TS_mat, &DFT_mat, &WS_mat, &dirs);CHKERRQ(ierr);

	// ierr = PetscSNPrintf((char*)&dirs.file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs.root_dir,dirs.output,"Y_hat1");CHKERRQ(ierr);
	// ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs.file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	// ierr = MatView(RSVD_mat.Y_hat,fd);CHKERRQ(ierr);

	// ierr = SlepcFinalize();
	// return ierr;

	ierr = PowerIterationRK4(&RSVD_mat, &RSVDt, &LNS_mat, &TS_mat, &DFT_mat, &WS_mat, &dirs);CHKERRQ(ierr);

	ierr = QR(&RSVD_mat, &RSVDt);CHKERRQ(ierr);

	// ierr = PetscSNPrintf((char*)&dirs.file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs.root_dir,dirs.output,"Q_hat1");CHKERRQ(ierr);
	// ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs.file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	// ierr = MatView(RSVD_mat.Y_hat,fd);CHKERRQ(ierr);

	ierr = StoreQ(&RSVD_mat, &RSVDt);CHKERRQ(ierr);

	ierr = AdjointActionRK4(&RSVD_mat, &RSVDt, &LNS_mat, &TS_mat, &DFT_mat, &WS_mat, &dirs);CHKERRQ(ierr);

	// ierr = PetscSNPrintf((char*)&dirs.file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs.root_dir,dirs.output,"Y_hat_before_SVD");CHKERRQ(ierr);
	// ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs.file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	// ierr = MatView(RSVD_mat.Y_hat,fd);CHKERRQ(ierr);

	ierr = ReducedSVD(&RSVD_mat, &RSVDt, &WS_mat, &Res_mat);CHKERRQ(ierr);

	ierr = SaveOutputs(&Res_mat, &dirs);CHKERRQ(ierr);

	ierr = SlepcFinalize();
	return ierr;

}

static PetscErrorCode AdjointActionRK4(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, LNS_params *LNS_mat, \
			TS_matrices *TS_mat, DFT_matrices *DFT_mat, WS_matrices *WS_mat, Directories *dirs)
{

	/*
		Performs time-stepping to approximate R^* \times \hat{F}
	*/
	
	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	RSVDt->TS.flg_dir_adj = 0;

	ierr = ApplyWSMats(RSVD_mat, RSVDt, WS_mat, 1);CHKERRQ(ierr);
	ierr = TSActionRK4(RSVD_mat, DFT_mat, LNS_mat, TS_mat, RSVDt, dirs);CHKERRQ(ierr);
	ierr = ApplyWSMats(RSVD_mat, RSVDt, WS_mat, 0);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode ApplyWSMats(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, WS_matrices *WS_mat, PetscBool before)
{

	/*
		Applies weight, input and output matrices if defined
	*/

	PetscErrorCode        ierr=0;

	PetscFunctionBeginUser;

	if (RSVDt->TS.flg_dir_adj) {
		if (before) {
			if (WS_mat->flg_Wout) {
				ierr = WS_mat->flg_diag ? MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->W_out_vec,NULL) : \
					MatMatMult(WS_mat->W_out,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);
			}
			if (WS_mat->flg_B) ierr = MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->B_vec,NULL);CHKERRQ(ierr);
		} else {
			if (WS_mat->flg_C) ierr = MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->C_vec,NULL);CHKERRQ(ierr);
			if (WS_mat->flg_Win) {
				ierr = WS_mat->flg_diag ? MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->W_in_vec,NULL) : \
					MatMatMult(WS_mat->W_in,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);
			}
		}
	} else {
		if (before) {
			if (WS_mat->flg_Win) {
				ierr = WS_mat->flg_diag ? MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->W_in_vec,NULL) : \
					MatMatMult(WS_mat->W_in,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);
			}
			if (WS_mat->flg_C) ierr = MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->C_vec,NULL);CHKERRQ(ierr);
		} else {
			if (WS_mat->flg_B) ierr = MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->B_vec,NULL);CHKERRQ(ierr);
			if (WS_mat->flg_Wout) {
				ierr = WS_mat->flg_diag ? MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->W_out_vec,NULL) : \
					MatMatMult(WS_mat->W_out,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);
			}
		}
	}

	PetscFunctionReturn(0);
	
}

static PetscErrorCode CreateDFTiDFTMats(RSVDt_vars *RSVDt, DFT_matrices *DFT_mat, LNS_params *LNS_mat)
{

	/*
		Creates the DFT and inverse DFT matrices
	*/	

	PetscErrorCode        ierr;
	PetscInt              iw, it, Nw, Nt, start, end, Nw_eff;
	PetscScalar           alpha, v;
	Vec                   V_temp, J;

	PetscFunctionBeginUser;

	Nw_eff = RSVDt->TS.flg_real_A ? RSVDt->RSVD.Nw/2 : RSVDt->RSVD.Nw;
	ierr = MatCreate(PETSC_COMM_WORLD,&DFT_mat->dft);CHKERRQ(ierr);
	ierr = MatSetType(DFT_mat->dft,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(DFT_mat->dft,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.Nw,Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(DFT_mat->dft);CHKERRQ(ierr); 

	ierr = MatCreate(PETSC_COMM_WORLD,&DFT_mat->inv_dft);CHKERRQ(ierr);
	ierr = MatSetType(DFT_mat->inv_dft,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(DFT_mat->inv_dft,PETSC_DECIDE,PETSC_DECIDE,2*RSVDt->TS.Nt,Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(DFT_mat->inv_dft);CHKERRQ(ierr);
	ierr = MatGetSize(DFT_mat->inv_dft,&Nt,&Nw);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &J);CHKERRQ(ierr);
	ierr = VecSetSizes(J, PETSC_DECIDE, Nt);CHKERRQ(ierr);
	ierr = VecSetUp(J);CHKERRQ(ierr);

	VecGetOwnershipRange(J,&start,&end);
	for (it=start; it<end; it++) {
		v = it;
		VecSetValues(J,1,&it,&v,INSERT_VALUES);
	}

	ierr = VecAssemblyBegin(J);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(J);CHKERRQ(ierr);

	alpha = -2*PETSC_PI*PETSC_i/Nt;

	for (iw=0; iw<Nw; iw++) {
		ierr = MatDenseGetColumnVecWrite(DFT_mat->inv_dft,iw,&V_temp);CHKERRQ(ierr);
		if (RSVDt->TS.flg_real_A) {
			v    = alpha*iw;
		} else if (PetscFmodReal(Nw,2) == 0)  {
			v    = (iw<=Nw/2-1) ? alpha*iw : alpha*((iw-Nw)+Nt);
		} else {
			v    = (iw<=(Nw-1)/2) ? alpha*iw : alpha*((iw-Nw)+Nt);
		}
		ierr = VecSet(V_temp,v);CHKERRQ(ierr);
		ierr = VecPointwiseMult(V_temp,V_temp,J);CHKERRQ(ierr);
		ierr = VecExp(V_temp);CHKERRQ(ierr);
		ierr = VecConjugate(V_temp);CHKERRQ(ierr);
		ierr = VecScale(V_temp, 2./Nt);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(DFT_mat->inv_dft,iw,&V_temp);CHKERRQ(ierr);
	}

	ierr = VecDestroy(&J);CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &J);CHKERRQ(ierr);
	ierr = RSVDt->TS.flg_real_A ? VecSetSizes(J,PETSC_DECIDE,Nw*2) : \
						VecSetSizes(J,PETSC_DECIDE,Nw);CHKERRQ(ierr);
	ierr = VecSetUp(J);CHKERRQ(ierr);

	VecGetOwnershipRange(J,&start,&end);
	for (it=start; it<end; it++) {
		v = it;
		VecSetValues(J,1,&it,&v,INSERT_VALUES);
	}

	ierr = VecAssemblyBegin(J);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(J);CHKERRQ(ierr);

	alpha = RSVDt->TS.flg_real_A ? -2*PETSC_PI*PETSC_i/(2*Nw) : -2*PETSC_PI*PETSC_i/Nw;
	
	for (iw=0; iw<Nw; iw++) {
		ierr = MatDenseGetColumnVecWrite(DFT_mat->dft,iw,&V_temp);CHKERRQ(ierr);
		v    = alpha*iw;
		ierr = VecSet(V_temp,v);CHKERRQ(ierr);
		ierr = VecPointwiseMult(V_temp,V_temp,J);CHKERRQ(ierr);
		ierr = VecExp(V_temp);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(DFT_mat->dft,iw,&V_temp);CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(DFT_mat->dft,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(DFT_mat->dft,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(DFT_mat->inv_dft,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(DFT_mat->inv_dft,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);	
	ierr = MatTranspose(DFT_mat->inv_dft,MAT_INPLACE_MATRIX,&DFT_mat->inv_dft);CHKERRQ(ierr);

	ierr = VecDestroy(&J);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode CreateRandomMat(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, Directories *dirs)
{

	/*
		Generates a random matrix of size N \times kNw (default)
		The input matrix can be read in using -Fin
		The latter option is useful for testing or resuming from a previous power iteration
	*/	

	PetscErrorCode        ierr;
	PetscBool             flg_Fin;
	PetscInt              Nw_eff;
	PetscViewer           fd;

	PetscFunctionBeginUser;

	Nw_eff = RSVDt->TS.flg_real_A ? RSVDt->RSVD.Nw/2 : RSVDt->RSVD.Nw;
	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->Y_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVD_mat->Y_hat,MATDENSE);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL,NULL,"-Fin",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,&flg_Fin);CHKERRQ(ierr);
	if (flg_Fin) {
		ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the input forcing matrix!\n");
		ierr = MatLoad(RSVD_mat->Y_hat,fd);CHKERRQ(ierr);
		ierr = MatGetSize(RSVD_mat->Y_hat,NULL,&RSVDt->RSVD.k);CHKERRQ(ierr);
		RSVDt->RSVD.k /= Nw_eff;
	} else {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Forcing input is generated randomly (uniform distribution)!\n");CHKERRQ(ierr);
		ierr = MatSetSizes(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,Nw_eff*RSVDt->RSVD.k);CHKERRQ(ierr);
		ierr = MatSetUp(RSVD_mat->Y_hat);CHKERRQ(ierr);
		ierr = MatSetRandom(RSVD_mat->Y_hat,NULL);CHKERRQ(ierr);		
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"@ N = %d, @ k = %d, @ Nw = %d\n",(int)RSVDt->RSVD.N,(int)RSVDt->RSVD.k,(int)Nw_eff);CHKERRQ(ierr);		

	PetscFunctionReturn(0);
	
}

static PetscErrorCode CreateForcingOnFly(Mat F_hat, DFT_matrices *DFT_mat, PetscInt jt_cyc, Vec F)
{

	/*
		Creates a temporary forcing snapshot in time from \hat{F} in Fourier space
	*/	

	PetscErrorCode        ierr;
	Vec                   W_col;

	PetscFunctionBeginUser;

	ierr = MatDenseGetColumnVecRead(DFT_mat->inv_dft,jt_cyc,&W_col);CHKERRQ(ierr);
	ierr = MatMult(F_hat,W_col,F);CHKERRQ(ierr);
	ierr = MatDenseRestoreColumnVecRead(DFT_mat->inv_dft,jt_cyc,&W_col);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode CreateTSMats(TS_matrices *TS_mat, PetscInt N)
{

	/*
		Creates the required RK4 time-stepping matrices 
	*/	

	PetscErrorCode      ierr;

	PetscFunctionBeginUser;

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->F1);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->F1,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->F1,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->F2);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->F2,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->F2,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->F3);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->F3,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->F3,PETSC_DECIDE,N);CHKERRQ(ierr); 

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->k1);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->k1,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->k1,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->k2);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->k2,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->k2,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->k3);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->k3,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->k3,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->k4);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->k4,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->k4,PETSC_DECIDE,N);CHKERRQ(ierr);	

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->y_temp);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->y_temp,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->y_temp,PETSC_DECIDE,N);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode DestroyTSMats(TS_matrices *TS_mat)
{

	/*
		Removes the RK4 time-stepping matrices from memory
	*/	

	PetscErrorCode      ierr;

	PetscFunctionBeginUser;

	ierr = VecDestroy(&TS_mat->F1);CHKERRQ(ierr);
	ierr = VecDestroy(&TS_mat->F2);CHKERRQ(ierr);
	ierr = VecDestroy(&TS_mat->F3);CHKERRQ(ierr);	
	ierr = VecDestroy(&TS_mat->k1);CHKERRQ(ierr);
	ierr = VecDestroy(&TS_mat->k2);CHKERRQ(ierr);
	ierr = VecDestroy(&TS_mat->k3);CHKERRQ(ierr);
	ierr = VecDestroy(&TS_mat->k4);CHKERRQ(ierr);
	ierr = VecDestroy(&TS_mat->y_temp);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode DFT(Mat Y_all, DFT_matrices *DFT_mat, RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt)
{

	/*
		Performs DFT on the collected response snapshots to obtain \hat{Y}
	*/	

	PetscErrorCode        ierr;
	Mat                   Y_temp1, Y_temp2;
	PetscInt              ik, Nw, N, k, Nw_eff;

	PetscFunctionBeginUser;

	ierr = MatGetSize(DFT_mat->dft,&Nw,&Nw_eff);CHKERRQ(ierr);
	ierr = MatGetSize(Y_all,&N,&k);CHKERRQ(ierr);
	k   /= Nw;
	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->Y_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVD_mat->Y_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,N,k*Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(RSVD_mat->Y_hat);CHKERRQ(ierr);

	for (ik=0; ik<k; ik++) {
		ierr = MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nw,(ik+1)*Nw,&Y_temp1);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,ik*Nw_eff,(ik+1)*Nw_eff,&Y_temp2);CHKERRQ(ierr);
		ierr = MatMatMult(Y_temp1,DFT_mat->dft,MAT_REUSE_MATRIX,PETSC_DEFAULT,&Y_temp2);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(RSVD_mat->Y_hat,&Y_temp2);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Y_all,&Y_temp1);CHKERRQ(ierr);
	}

	ierr = MatDestroy(&Y_all);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(RSVD_mat->Y_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(RSVD_mat->Y_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatScale(RSVD_mat->Y_hat, RSVDt->TS.time_ratio);CHKERRQ(ierr);
	ierr = ReversePermuteMat(RSVD_mat->Y_hat, RSVDt);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode ForwardActionRK4(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, LNS_params *LNS_mat, \
			TS_matrices *TS_mat, DFT_matrices *DFT_mat, WS_matrices *WS_mat, Directories *dirs)
{

	/*
		Performs time-stepping to approximate $R \times \hat{F}$ 
	*/

	PetscErrorCode        ierr;

	PetscFunctionBeginUser;
		
	RSVDt->TS.flg_dir_adj = 1;

	ierr = ApplyWSMats(RSVD_mat, RSVDt, WS_mat, 1);CHKERRQ(ierr);
	ierr = TSActionRK4(RSVD_mat, DFT_mat, LNS_mat, TS_mat, RSVDt, dirs);CHKERRQ(ierr);
	ierr = ApplyWSMats(RSVD_mat, RSVDt, WS_mat, 0);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode LNSStructure(LNS_params *LNS_mat, TS_matrices *TS_mat, Directories *dirs)
{

	/*
		Reads the LNS operator and applies discounting if desired
	*/

	PetscErrorCode      ierr;
	PetscInt            N;
	PetscLogDouble      t1, t2;
	PetscViewer         fd;

	PetscFunctionBeginUser;

	ierr = PetscTime(&t1);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->Mat_filename);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&LNS_mat->A);CHKERRQ(ierr);
	ierr = MatSetType(LNS_mat->A,MATMPIAIJ);CHKERRQ(ierr);
	ierr = MatLoad(LNS_mat->A,fd);CHKERRQ(ierr);
	ierr = MatGetSize(LNS_mat->A,&N,NULL);CHKERRQ(ierr);

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"The matrix is loaded in %f seconds (N = %d)\n", t2-t1, (int)N);CHKERRQ(ierr);

	/*
		Discounting
	*/

	if (LNS_mat->RSVDt.disc.flg_disc) {
		ierr = MatShift(LNS_mat->A,-LNS_mat->RSVDt.disc.beta);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"**** Discounting with beta = %g ****\n", LNS_mat->RSVDt.disc.beta);CHKERRQ(ierr);
	}

	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}

static PetscErrorCode MatExpActionRK4(Mat Y_in, Mat Y_out, PetscInt Ngap, PetscInt st_ind, LNS_params *LNS_mat, \
					WS_matrices *WS_mat, TS_matrices *TS_mat, RSVDt_vars *RSVDt, Directories *dirs)
{

	PetscErrorCode       ierr;
	Vec                  y0,y_temp;
	PetscInt             N,k,i,ik,rend=st_ind+Ngap,prg_cnt=0;
	PetscViewer          fd;
	PetscLogDouble       t1,t2;

	PetscFunctionBeginUser;

	ierr = RSVDt->TS.flg_dir_adj ? PetscPrintf(PETSC_COMM_WORLD,"\n*** Computing the subspace evolution (direct action)! ***\n") : \
				PetscPrintf(PETSC_COMM_WORLD,"\n*** Computing the subspace evolution (adjoint action)! ***\n");CHKERRQ(ierr);

	/*
		Create and set the size of matrices/vecs
	*/

	ierr = MatGetSize(Y_in,&N,&k);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** time interval = %f , subspace dim = %d ***\n", (double)Ngap*RSVDt->TS.dt, (int)k);CHKERRQ(ierr);

	ierr = MatSetType(Y_out,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Y_out,PETSC_DECIDE,PETSC_DECIDE,N,k);CHKERRQ(ierr);
	ierr = MatSetUp(Y_out);CHKERRQ(ierr);	

	ierr = VecCreate(PETSC_COMM_WORLD,&y0);CHKERRQ(ierr);
	ierr = VecSetType(y0,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(y0,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = CreateTSMats(TS_mat,N);CHKERRQ(ierr);

	/*
		Loop over all test vectors
	*/

	for (ik=0; ik<k; ik++) {

		/*
			Extract IC's for each ik
		*/

		ierr = MatDenseGetColumnVecRead(Y_in,ik,&y_temp);CHKERRQ(ierr);
		ierr = VecCopy(y_temp,y0);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecRead(Y_in,ik,&y_temp);CHKERRQ(ierr);

		/*
			Perform time-stepping
		*/

		ierr = PetscTime(&t1);CHKERRQ(ierr);

		for (i=st_ind; i<rend; i++) {

			ierr = TSTransRK4(LNS_mat,TS_mat,RSVDt,i,y0);CHKERRQ(ierr);

			if ((double)(i-st_ind)/(rend-st_ind) >= 0.01*prg_cnt){
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Progress percentage: %d% \n",(int)prg_cnt);CHKERRQ(ierr);
				ierr = PetscTime(&t2);CHKERRQ(ierr);
				ierr = PetscPrintf(PETSC_COMM_WORLD,"TS took %f seconds (wall-time)\n", t2-t1);CHKERRQ(ierr);		
				prg_cnt += 10;
			}
		}

		ierr = PetscTime(&t2);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Transient took %f seconds! (i = %d)\n", t2-t1, (int)ik+1);CHKERRQ(ierr);

		ierr = MatDenseGetColumnVecWrite(Y_out,ik,&y_temp);CHKERRQ(ierr);
		ierr = VecCopy(y0,y_temp);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(Y_out,ik,&y_temp);CHKERRQ(ierr);	
	}

	/*
		Save the output
	*/	

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"EQ_test");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = MatView(Y_out,fd);CHKERRQ(ierr);

	/*
		Remove local variables
	*/

	ierr = VecDestroy(&y0);CHKERRQ(ierr);
	ierr = DestroyTSMats(TS_mat);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}

static PetscErrorCode MatExpRunRK4(RSVDt_vars *RSVDt, LNS_params *LNS_mat, WS_matrices *WS_mat, TS_matrices *TS_mat, Directories *dirs)
{
	
	PetscErrorCode        ierr;
	Mat                   Y_in, Y_out;
	PetscBool             period_flg, MatExpRun_flg, MatExpRun_dir;
	PetscInt              Ngap, st_ind = 1;
	PetscViewer           fd;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetString(NULL,NULL,"-input_subspace",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-MatExpRun_dir",&MatExpRun_dir,&MatExpRun_dir);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-Ngap_period_flg",&period_flg,&period_flg);CHKERRQ(ierr);

	RSVDt->TS.flg_dir_adj = MatExpRun_dir;
	
	ierr = MatCreate(PETSC_COMM_WORLD,&Y_out);CHKERRQ(ierr);
	Ngap = period_flg ? RSVDt->TS.Nt : RSVDt->TS.time_ratio;
	ierr = MatCreate(PETSC_COMM_WORLD,&Y_in);CHKERRQ(ierr);
	ierr = MatSetType(Y_in,MATDENSE);CHKERRQ(ierr);
	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatLoad(Y_in,fd);CHKERRQ(ierr);
	ierr = MatExpActionRK4(Y_in,Y_out,Ngap,st_ind,LNS_mat,WS_mat,TS_mat,RSVDt,dirs);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

// static PetscErrorCode MatExpWForcingRK4(RSVD_matrices *RSVD_mat, DFT_matrices *DFT_mat, \
// 		LNS_params *LNS_mat, TS_matrices *TS_mat, RSVDt_vars *RSVDt, Directories *dirs)
// {

// 	PetscErrorCode       ierr;
// 	Mat                  Y_all,F_hat_re,F_temp,Y_all_k,F_hat_k,Y_IC;
// 	Vec                  y,y_temp;
// 	PetscInt             N,k,Nt_saved,Nw_eff,i,ik,rend,prg_cnt,pos,mod,st_ind=1;
// 	PetscReal            norm;
// 	PetscBool            flg_IC;
// 	PetscViewer          fd;
// 	PetscLogDouble       t1,t2;

// 	PetscFunctionBeginUser;

// 	ierr = RSVDt->TS.flg_dir_adj ? PetscPrintf(PETSC_COMM_WORLD,"\n*** Direct action in process! ***\n") : \
// 				PetscPrintf(PETSC_COMM_WORLD,"\n*** Adjoint action in process! ***\n");CHKERRQ(ierr);

// 	/*
// 		Define TS parameters
// 	*/

// 	// rend  = RSVDt->TS.time_ratio;
// 	rend  = RSVDt->TS.Nt_transient;
// 	mod   = RSVDt->TS.time_ratio;
// 	Nt_saved = PetscFloorReal(rend/mod) + 1;
// 	rend += st_ind;

// 	/*
// 		Get the size of matrices
// 	*/

// 	ierr = MatGetSize(RSVD_mat->Y_hat,&N,&k);CHKERRQ(ierr);
// 	Nw_eff = RSVDt->TS.flg_real_A ? RSVDt->RSVD.Nw/2 : RSVDt->RSVD.Nw;
// 	k   /= Nw_eff;
// 	ierr = PetscPrintf(PETSC_COMM_WORLD,"Total snapshots to be saved = %d for k = %d (The first snapshot is the I.C.)\n", (int)Nt_saved, (int)k);CHKERRQ(ierr);

// 	/*
// 		Apply solution from previous action as forcing
// 	*/

// 	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->F_hat);CHKERRQ(ierr);
// 	ierr = MatSetType(RSVD_mat->F_hat,MATDENSE);CHKERRQ(ierr);
// 	ierr = MatSetSizes(RSVD_mat->F_hat,PETSC_DECIDE,PETSC_DECIDE,N,Nw_eff*k);CHKERRQ(ierr);
// 	ierr = MatSetUp(RSVD_mat->F_hat);CHKERRQ(ierr);	
// 	ierr = MatCopy(RSVD_mat->Y_hat,RSVD_mat->F_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
// 	ierr = MatDestroy(&RSVD_mat->Y_hat);CHKERRQ(ierr);		

// 	/*
// 		Create all required matrices (including solution, forcing, RHS, norm, etc.)
// 	*/

// 	ierr = VecCreate(PETSC_COMM_WORLD,&y);CHKERRQ(ierr);
// 	ierr = VecSetType(y,VECMPI);CHKERRQ(ierr);
// 	ierr = VecSetSizes(y,PETSC_DECIDE,N);CHKERRQ(ierr);

// 	ierr = CreateTSMats(TS_mat,N);CHKERRQ(ierr);

// 	ierr = MatCreate(PETSC_COMM_WORLD,&Y_all_k);CHKERRQ(ierr);
// 	ierr = MatSetType(Y_all_k,MATDENSE);CHKERRQ(ierr);
// 	ierr = MatSetSizes(Y_all_k,PETSC_DECIDE,PETSC_DECIDE,N,Nt_saved);CHKERRQ(ierr);
// 	ierr = MatSetUp(Y_all_k);CHKERRQ(ierr);

// 	ierr = MatCreate(PETSC_COMM_WORLD,&Y_IC);CHKERRQ(ierr);
// 	ierr = MatSetType(Y_IC,MATDENSE);CHKERRQ(ierr);
// 	ierr = PetscOptionsGetString(NULL,NULL,"-IC",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,&flg_IC);CHKERRQ(ierr);
// 	if (flg_IC) {
// 		ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,dirs->filename);CHKERRQ(ierr);
// 		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);			
// 		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading IC before forced run!\n");CHKERRQ(ierr);
// 		ierr = MatLoad(Y_IC,fd);CHKERRQ(ierr);
// 		ierr = MatGetSize(Y_IC,NULL,&k);CHKERRQ(ierr);
// 		ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of test vectors = %d\n", (int)k);CHKERRQ(ierr);
// 	} else {
// 		ierr = PetscPrintf(PETSC_COMM_WORLD,"Starting from zero IC for all test vectors!\n");CHKERRQ(ierr);
// 		ierr = MatSetSizes(Y_IC,PETSC_DECIDE,PETSC_DECIDE,N,k);CHKERRQ(ierr);		
// 	}
// 	ierr = MatSetUp(Y_IC);CHKERRQ(ierr);

// 	ierr = MatCreate(PETSC_COMM_WORLD,&Y_all);CHKERRQ(ierr);
// 	ierr = MatSetType(Y_all,MATDENSE);CHKERRQ(ierr);
// 	ierr = MatSetSizes(Y_all,PETSC_DECIDE,PETSC_DECIDE,N,Nt_saved*k);CHKERRQ(ierr);
// 	ierr = MatSetUp(Y_all);CHKERRQ(ierr);

// 	/*
// 		Reshape the forcing matrix
// 	*/

// 	ierr = MatCreate(PETSC_COMM_WORLD,&F_hat_re);CHKERRQ(ierr);
// 	ierr = MatSetType(F_hat_re,MATDENSE);CHKERRQ(ierr);
// 	ierr = PermuteMat(RSVD_mat->F_hat, Nw_eff, F_hat_re);CHKERRQ(ierr);

// 	ierr = MatCreate(PETSC_COMM_WORLD,&F_hat_k);CHKERRQ(ierr);
// 	ierr = MatSetType(F_hat_k,MATDENSE);CHKERRQ(ierr);
// 	ierr = MatSetSizes(F_hat_k,PETSC_DECIDE,PETSC_DECIDE,N,Nw_eff);CHKERRQ(ierr);
// 	ierr = MatSetUp(F_hat_k);CHKERRQ(ierr);

// 	/*
// 		Extract the forcing submatrix for each mode and loop-over all modes
// 	*/	

// 	for (ik=0; ik<k; ik++) {

// 		prg_cnt = 0;

// 		ierr = MatDenseGetColumnVecWrite(Y_IC,ik,&y_temp);CHKERRQ(ierr);
// 		ierr = VecCopy(y_temp,y);CHKERRQ(ierr);
// 		ierr = MatDenseRestoreColumnVecWrite(Y_IC,ik,&y_temp);CHKERRQ(ierr);

// 		ierr = MatDenseGetSubMatrix(F_hat_re,PETSC_DECIDE,PETSC_DECIDE,ik*Nw_eff,(ik+1)*Nw_eff,&F_temp);CHKERRQ(ierr);
// 		ierr = MatCopy(F_temp,F_hat_k,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
// 		ierr = MatAssemblyBegin(F_hat_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
// 		ierr = MatAssemblyEnd(F_hat_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
// 		ierr = MatDenseRestoreSubMatrix(F_hat_re,&F_temp);CHKERRQ(ierr);

// 		/*
// 			Create initial forcing
// 		*/

// 		ierr = CreateForcingOnFly(F_hat_k,DFT_mat->inv_dft,0,TS_mat->F3);CHKERRQ(ierr);

// 		/*
// 			Measure time-stepping wall-time
// 		*/

// 		ierr = PetscTime(&t1);CHKERRQ(ierr);

// 		/*
// 			Perform time-stepping (save "initial conditions" instead of "responses")
// 		*/

// 		for (i=st_ind; i<rend; i++) {

// 			if (PetscFmodReal(i-st_ind,mod) == 0) {
// 				pos  = (i-st_ind)/mod;
// 				ierr = MatDenseGetColumnVecWrite(Y_all_k,pos,&y_temp);CHKERRQ(ierr);
// 				ierr = VecCopy(y,y_temp);CHKERRQ(ierr);
// 				ierr = MatDenseRestoreColumnVecWrite(Y_all_k,pos,&y_temp);CHKERRQ(ierr);

// 				ierr = VecNorm(y,NORM_2,&norm);CHKERRQ(ierr);
// 				ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of the solution: %f @ position index = %d\n",(double)norm,(int)PetscFmodReal(pos,Nw_eff)+1);CHKERRQ(ierr);
// 			}

// 			ierr = TSRK4(F_hat_k,DFT_mat,LNS_mat,TS_mat,RSVDt,i,y);CHKERRQ(ierr);

// 			if ((double)(i-st_ind+1)/(rend-st_ind+1) >= 0.01*prg_cnt){
// 				ierr = PetscPrintf(PETSC_COMM_WORLD,"Progress percentage: %d% \n",(int)prg_cnt);CHKERRQ(ierr);
// 				ierr = PetscTime(&t2);CHKERRQ(ierr);
// 				ierr = PetscPrintf(PETSC_COMM_WORLD,"TS took %f seconds (wall-time)\n", t2-t1);CHKERRQ(ierr);		
// 				prg_cnt += 10;
// 			}
// 		}

// 		ierr = PetscTime(&t2);CHKERRQ(ierr);
// 		ierr = PetscPrintf(PETSC_COMM_WORLD,"TS took %f seconds (k = %d)\n", t2-t1, (int)ik);CHKERRQ(ierr);

// 		ierr = MatAssemblyBegin(Y_all_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
// 		ierr = MatAssemblyEnd(Y_all_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
// 		ierr = MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nt_saved,(ik+1)*Nt_saved,&F_temp);CHKERRQ(ierr);
// 		ierr = MatCopy(Y_all_k,F_temp,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
// 		ierr = MatDenseRestoreSubMatrix(Y_all,&F_temp);CHKERRQ(ierr);	

// 	}	

// 	ierr = RSVDt->TS.flg_dir_adj ? PetscPrintf(PETSC_COMM_WORLD,"Direct action DONE!\n") : \
// 				PetscPrintf(PETSC_COMM_WORLD,"Adjoint action DONE!\n");CHKERRQ(ierr);		

// 	/*
// 		Remove the forcing input as it is no longer needed
// 	*/

// 	ierr = VecDestroy(&y);CHKERRQ(ierr);
// 	ierr = MatDestroy(&F_hat_re);CHKERRQ(ierr);
// 	ierr = DestroyTSMats(TS_mat);CHKERRQ(ierr);

// 	/*
// 		Save the output
// 	*/

// 	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"Y_forced_dir");CHKERRQ(ierr);
// 	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
// 	ierr = MatView(Y_all,fd);CHKERRQ(ierr);			

// 	PetscFunctionReturn(0);

// }

static PetscErrorCode PermuteMat(Mat F, RSVDt_vars *RSVDt)
{

	/*
		Permutes a matrix of size N \times kNw to N \times Nwk (in place)
		Nw blocks of N \times k   ---->   k blocks of N \times Nw
	*/

	PetscErrorCode        ierr;
	PetscInt              iw, ik, ii, jj, N, k;
	Mat                   F_permuted;
	Vec                   V1, V2;

	PetscFunctionBeginUser;

	ierr = MatGetSize(F,&N,&k);CHKERRQ(ierr);
	k   /= RSVDt->RSVD.Nw;
	ierr = MatCreate(PETSC_COMM_WORLD,&F_permuted);CHKERRQ(ierr);
	ierr = MatSetType(F_permuted,MATMPIDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(F_permuted,PETSC_DECIDE,PETSC_DECIDE,N,RSVDt->RSVD.Nw*k);CHKERRQ(ierr);
	ierr = MatSetUp(F_permuted);CHKERRQ(ierr);

	for (ik=0; ik<k; ik++) {
		for (iw=0; iw<RSVDt->RSVD.Nw; iw++) {
			ii   = ik                + iw*k;
			jj   = ik*RSVDt->RSVD.Nw + iw;
			ierr = MatDenseGetColumnVecRead(F,ii,&V1);CHKERRQ(ierr);
			ierr = MatDenseGetColumnVecWrite(F_permuted,jj,&V2);CHKERRQ(ierr);
			ierr = VecCopy(V1,V2);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(F_permuted,jj,&V2);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecRead(F,ii,&V1);CHKERRQ(ierr);
		}
	}

	ierr = MatAssemblyBegin(F_permuted,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(F_permuted,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatCopy(F_permuted,F,SAME_NONZERO_PATTERN);CHKERRQ(ierr);

	ierr = MatDestroy(&F_permuted);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode PowerIterationRK4(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, LNS_params *LNS_mat, \
			TS_matrices *TS_mat, DFT_matrices *DFT_mat, WS_matrices *WS_mat, Directories *dirs)
{

	/*
		Performs power iteration for q times 
	*/

	PetscErrorCode        ierr;
	PetscInt              iq, Nw_eff;

	PetscFunctionBeginUser;	

	Nw_eff = RSVDt->TS.flg_real_A ? RSVDt->RSVD.Nw/2 : RSVDt->RSVD.Nw;

	for (iq=0; iq<RSVDt->RSVD.q; iq++) {

		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n ******* Inside power iteration, %d/%d *******",(int)iq+1,(int)RSVDt->RSVD.q);CHKERRQ(ierr);

		ierr = QR(RSVD_mat, RSVDt);CHKERRQ(ierr);
		ierr = AdjointActionRK4(RSVD_mat, RSVDt, LNS_mat, TS_mat, DFT_mat, WS_mat, dirs);CHKERRQ(ierr);
		ierr = QR(RSVD_mat, RSVDt);CHKERRQ(ierr);
		ierr = ForwardActionRK4(RSVD_mat, RSVDt, LNS_mat, TS_mat, DFT_mat, WS_mat, dirs);CHKERRQ(ierr);
		
	}

	PetscFunctionReturn(0);
	
}

static PetscErrorCode QR(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt)
{

	/*
		Performs QR decomposition for all Nw matrices of size N \times k 
	*/	

	PetscErrorCode        ierr;
	PetscInt              iw, Nw_eff;
	Mat                   Y_temp, Q_temp;
	BV                    Q;

	PetscFunctionBeginUser;

	Nw_eff = RSVDt->TS.flg_real_A ? RSVDt->RSVD.Nw/2 : RSVDt->RSVD.Nw;
	for (iw=0; iw<Nw_eff; iw++) {
		ierr = MatDenseGetSubMatrix(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&Y_temp);CHKERRQ(ierr);
		ierr = BVCreateFromMat(Y_temp,&Q);CHKERRQ(ierr);
		ierr = BVSetType(Q,BVVECS);CHKERRQ(ierr);
		ierr = BVOrthogonalize(Q,NULL);CHKERRQ(ierr);
		ierr = BVCreateMat(Q,&Q_temp);CHKERRQ(ierr);
		ierr = MatCopy(Q_temp,Y_temp,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(RSVD_mat->Y_hat,&Y_temp);CHKERRQ(ierr);
		ierr = BVDestroy(&Q);CHKERRQ(ierr);
		ierr = MatDestroy(&Q_temp);CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);

}

static PetscErrorCode ReadUserInput(RSVDt_vars *RSVDt, WS_matrices *WS_mat, LNS_params *LNS_mat, \
									Function_calls *Funcs, Directories *dirs)
{

	/*
		Reads in user input parameters
	*/

	PetscErrorCode        ierr;

	PetscFunctionBeginUser;
	
	ierr = PetscOptionsGetReal(NULL,NULL,"-t_transient",&RSVDt->TS.t_transient,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,NULL,"-dt_max",&RSVDt->TS.dt_max,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,NULL,"-w_min",&RSVDt->RSVD.w_min,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-TransientRun_flg",&Funcs->TransientRun_flg,&Funcs->TransientRun_flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-MatExpRun_flg",&Funcs->MatExpRun_flg,&Funcs->MatExpRun_flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-TSTimeDomainRun_flg",&Funcs->TSTimeDomainRun_flg,&Funcs->TSTimeDomainRun_flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-diagW",&WS_mat->flg_diag,&WS_mat->flg_diag);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-real_A",&RSVDt->TS.flg_real_A,&RSVDt->TS.flg_real_A);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-eff_tran",&RSVDt->TS.flg_eff_trans,&RSVDt->TS.flg_eff_trans);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-discounting",&LNS_mat->RSVDt.disc.flg_disc,&LNS_mat->RSVDt.disc.flg_disc);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,NULL,"-beta",&LNS_mat->RSVDt.disc.beta,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL,"-root_dir",(char*)&dirs->root_dir,PETSC_MAX_PATH_LEN,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL,"-input_folder",(char*)&dirs->input,PETSC_MAX_PATH_LEN,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL,"-output_folder",(char*)&dirs->output,PETSC_MAX_PATH_LEN,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL,"-Mat_filename",(char*)&dirs->Mat_filename,PETSC_MAX_PATH_LEN,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-k",&RSVDt->RSVD.k,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-N",&RSVDt->RSVD.N,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-q",&RSVDt->RSVD.q,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-Nw",&RSVDt->RSVD.Nw,NULL);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}

static PetscErrorCode ReadWSMats(WS_matrices *WS_mat, Directories *dirs)
{

	/*
		Read the weight, input and outputs matrices if given in the jobscript
	*/

	PetscErrorCode        ierr;
	PetscInt              N;
	Mat                   B,C;
	PetscViewer           fd;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetString(NULL,NULL,"-Win",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,&WS_mat->flg_Win);CHKERRQ(ierr);
	if (WS_mat->flg_Win) {
		ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);		
		ierr = MatCreate(PETSC_COMM_WORLD,&WS_mat->W_in);CHKERRQ(ierr);
		ierr = MatSetType(WS_mat->W_in,MATMPIAIJ);CHKERRQ(ierr);		
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the input weight matrix!\n");
		ierr = MatLoad(WS_mat->W_in,fd);CHKERRQ(ierr);
		if (WS_mat->flg_diag) {
			ierr = MatGetSize(WS_mat->W_in,&N,NULL);CHKERRQ(ierr);
			ierr = VecCreate(PETSC_COMM_WORLD,&WS_mat->W_in_vec);CHKERRQ(ierr);
			ierr = VecSetSizes(WS_mat->W_in_vec,PETSC_DECIDE,N);CHKERRQ(ierr);
			ierr = VecSetUp(WS_mat->W_in_vec);CHKERRQ(ierr);
			ierr = MatGetDiagonal(WS_mat->W_in,WS_mat->W_in_vec);CHKERRQ(ierr);
			ierr = MatDestroy(&WS_mat->W_in);CHKERRQ(ierr);
		}		
	}

	ierr = PetscOptionsGetString(NULL,NULL,"-Wout",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,&WS_mat->flg_Wout);CHKERRQ(ierr);	
	if (WS_mat->flg_Wout) {
		ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);	
		ierr = MatCreate(PETSC_COMM_WORLD,&WS_mat->W_out);CHKERRQ(ierr);
		ierr = MatSetType(WS_mat->W_out,MATMPIAIJ);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the output weight matrix!\n");
		ierr = MatLoad(WS_mat->W_out,fd);CHKERRQ(ierr);
		if (WS_mat->flg_diag) {
			ierr = MatGetSize(WS_mat->W_out,&N,NULL);CHKERRQ(ierr);
			ierr = VecCreate(PETSC_COMM_WORLD,&WS_mat->W_out_vec);CHKERRQ(ierr);
			ierr = VecSetSizes(WS_mat->W_out_vec,PETSC_DECIDE,N);CHKERRQ(ierr);
			ierr = VecSetUp(WS_mat->W_out_vec);CHKERRQ(ierr);
			ierr = MatGetDiagonal(WS_mat->W_out,WS_mat->W_out_vec);CHKERRQ(ierr);
			ierr = MatDestroy(&WS_mat->W_out);CHKERRQ(ierr);
		}		
	}

	ierr = PetscOptionsGetString(NULL,NULL,"-B",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,&WS_mat->flg_B);CHKERRQ(ierr);
	if (WS_mat->flg_B) {
		ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);
		ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
		ierr = MatSetType(B,MATMPIAIJ);CHKERRQ(ierr);		
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the forcing spatial filter matrix!\n");
		ierr = MatLoad(B,fd);CHKERRQ(ierr);
		ierr = MatGetSize(B,&N,NULL);CHKERRQ(ierr);
		ierr = VecCreate(PETSC_COMM_WORLD,&WS_mat->B_vec);CHKERRQ(ierr);
		ierr = VecSetSizes(WS_mat->B_vec,PETSC_DECIDE,N);CHKERRQ(ierr);
		ierr = VecSetUp(WS_mat->B_vec);CHKERRQ(ierr);
		ierr = MatGetDiagonal(B,WS_mat->B_vec);CHKERRQ(ierr);
		ierr = MatDestroy(&B);CHKERRQ(ierr);		
	}

	ierr = PetscOptionsGetString(NULL,NULL,"-C",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,&WS_mat->flg_C);CHKERRQ(ierr);
	if (WS_mat->flg_C) {
		ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);
		ierr = MatCreate(PETSC_COMM_WORLD,&C);CHKERRQ(ierr);
		ierr = MatSetType(C,MATMPIAIJ);CHKERRQ(ierr);		
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the response spatial filter matrix!\n");
		ierr = MatLoad(C,fd);CHKERRQ(ierr);
		ierr = MatGetSize(C,&N,NULL);CHKERRQ(ierr);
		ierr = VecCreate(PETSC_COMM_WORLD,&WS_mat->C_vec);CHKERRQ(ierr);
		ierr = VecSetSizes(WS_mat->C_vec,PETSC_DECIDE,N);CHKERRQ(ierr);
		ierr = VecSetUp(WS_mat->C_vec);CHKERRQ(ierr);
		ierr = MatGetDiagonal(C,WS_mat->C_vec);CHKERRQ(ierr);
		ierr = MatDestroy(&C);CHKERRQ(ierr);		
	}

	PetscFunctionReturn(0);
	
}

static PetscErrorCode ReducedSVD(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, WS_matrices *WS_mat, Resolvent_matrices *Res_mat)
{

	/*
		Performs the economy SVD of matrices of size N \times k for Nw frequencies
	*/
	
	PetscErrorCode        ierr;
	PetscInt              iw, ik, Nw_eff;
	PetscReal             sigma;
	Mat                   Y_temp, U_temp, V_temp, U_til;
	Vec                   U, V, Sig;
	SVD                   svd;

	PetscFunctionBeginUser;

	Nw_eff = RSVDt->TS.flg_real_A ? RSVDt->RSVD.Nw/2 : RSVDt->RSVD.Nw;
	ierr = MatCreate(PETSC_COMM_WORLD,&U_til);CHKERRQ(ierr);
	ierr = MatSetType(U_til,MATMPIDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(U_til,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.k,RSVDt->RSVD.k*Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(U_til);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&Res_mat->V_hat);CHKERRQ(ierr);
	ierr = MatSetType(Res_mat->V_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Res_mat->V_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,RSVDt->RSVD.k*Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(Res_mat->V_hat);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&Res_mat->S_hat);CHKERRQ(ierr);
	ierr = MatSetType(Res_mat->S_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Res_mat->S_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.k,Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(Res_mat->S_hat);CHKERRQ(ierr);

	for (iw=0; iw<Nw_eff; iw++) {

		ierr = MatDenseGetSubMatrix(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&Y_temp);CHKERRQ(ierr);
		ierr = SVDCreate(PETSC_COMM_WORLD,&svd);CHKERRQ(ierr);
		ierr = SVDSetOperators(svd,Y_temp,NULL);CHKERRQ(ierr);
		ierr = SVDSetDimensions(svd,RSVDt->RSVD.k,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
		ierr = SVDSolve(svd);CHKERRQ(ierr);

		ierr = MatDenseGetSubMatrix(U_til,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&U_temp);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(Res_mat->V_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&V_temp);CHKERRQ(ierr);
		ierr = MatDenseGetColumnVecWrite(Res_mat->S_hat,iw,&Sig);CHKERRQ(ierr);
		for (ik=0; ik<RSVDt->RSVD.k; ik++){
			ierr = MatDenseGetColumnVecWrite(U_temp,ik,&U);CHKERRQ(ierr);
			ierr = MatDenseGetColumnVecWrite(V_temp,ik,&V);CHKERRQ(ierr);
			ierr = SVDGetSingularTriplet(svd,ik,&sigma,V,U);
			ierr = VecSetValue(Sig,ik,sigma,INSERT_VALUES);
			ierr = MatDenseRestoreColumnVecWrite(U_temp,ik,&U);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(V_temp,ik,&V);CHKERRQ(ierr);
			ierr = VecAssemblyBegin(Sig);CHKERRQ(ierr);
			ierr = VecAssemblyEnd(Sig);CHKERRQ(ierr);
		}
		ierr = MatDenseRestoreColumnVecWrite(Res_mat->S_hat,iw,&Sig);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(U_til,&U_temp);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Res_mat->V_hat,&V_temp);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(RSVD_mat->Y_hat,&Y_temp);CHKERRQ(ierr);
		ierr = SVDDestroy(&svd);CHKERRQ(ierr);
	}

	ierr = MatDestroy(&RSVD_mat->Y_hat);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(Res_mat->V_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Res_mat->V_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatConjugate(Res_mat->V_hat);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(U_til,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(U_til,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatConjugate(U_til);CHKERRQ(ierr);

	ierr = RecoverU(U_til,RSVDt,RSVD_mat,Res_mat);CHKERRQ(ierr);

	if (WS_mat->flg_Wout) {
		if (WS_mat->flg_diag) {
			ierr = MatDiagonalScale(Res_mat->U_hat,WS_mat->W_out_vec,NULL);CHKERRQ(ierr);
			ierr = MatDiagonalScale(Res_mat->V_hat,WS_mat->W_out_vec,NULL);CHKERRQ(ierr);
		}
		else {
			ierr = MatMatMult(WS_mat->W_out,Res_mat->U_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&Res_mat->U_hat);CHKERRQ(ierr);
			ierr = MatMatMult(WS_mat->W_out,Res_mat->V_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&Res_mat->V_hat);CHKERRQ(ierr);
		}
	}

	PetscFunctionReturn(0);

}

static PetscErrorCode RecoverU(Mat U_til, RSVDt_vars *RSVDt, RSVD_matrices *RSVD_mat, Resolvent_matrices *Res_mat)
{

	/*
		Recovers \hat{U}  for all Nw matrices of size N \times k 
	*/

	PetscErrorCode        ierr;
	PetscInt              iw, Nw_eff;
	Mat                   Q_temp, U_temp1, U_temp2;

	PetscFunctionBeginUser;

	Nw_eff = RSVDt->TS.flg_real_A ? RSVDt->RSVD.Nw/2 : RSVDt->RSVD.Nw;
	ierr = MatCreate(PETSC_COMM_WORLD,&Res_mat->U_hat);CHKERRQ(ierr);
	ierr = MatSetType(Res_mat->U_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Res_mat->U_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,RSVDt->RSVD.k*Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(Res_mat->U_hat);CHKERRQ(ierr);

	for (iw=0; iw<Nw_eff; iw++) {
		ierr = MatDenseGetSubMatrix(RSVD_mat->Q_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&Q_temp);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(U_til,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&U_temp1);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(Res_mat->U_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&U_temp2);CHKERRQ(ierr);
		ierr = MatMatMult(Q_temp,U_temp1,MAT_REUSE_MATRIX,PETSC_DEFAULT,&U_temp2);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Res_mat->U_hat,&U_temp2);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(U_til,&U_temp1);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(RSVD_mat->Q_hat,&Q_temp);CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(Res_mat->U_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Res_mat->U_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = MatDestroy(&RSVD_mat->Q_hat);CHKERRQ(ierr);
	ierr = MatDestroy(&U_til);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode ReversePermuteMat(Mat F, RSVDt_vars *RSVDt)
{

	/*
		Permutes a matrix of size N \times Nwk to N \times kNw (in place)
		k blocks of N \times Nw   ---->   Nw blocks of N \times k 
	*/

	PetscErrorCode        ierr;
	PetscInt              iw, ik, ii, jj, N, k;
	Mat                   F_permuted;
	Vec                   V1, V2;

	PetscFunctionBeginUser;

	ierr = MatGetSize(F,&N,&k);CHKERRQ(ierr);
	k   /= RSVDt->RSVD.Nw;
	ierr = MatCreate(PETSC_COMM_WORLD,&F_permuted);CHKERRQ(ierr);
	ierr = MatSetType(F_permuted,MATMPIDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(F_permuted,PETSC_DECIDE,PETSC_DECIDE,N,RSVDt->RSVD.Nw*k);CHKERRQ(ierr);
	ierr = MatSetUp(F_permuted);CHKERRQ(ierr);

	for (ik=0; ik<k; ik++) {
		for (iw=0; iw<RSVDt->RSVD.Nw; iw++) {
			ii   = ik    + iw*k;
			jj   = ik*RSVDt->RSVD.Nw + iw;
			ierr = MatDenseGetColumnVecRead(F,jj,&V1);CHKERRQ(ierr);
			ierr = MatDenseGetColumnVecWrite(F_permuted,ii,&V2);CHKERRQ(ierr);
			ierr = VecCopy(V1,V2);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(F_permuted,ii,&V2);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecRead(F,jj,&V1);CHKERRQ(ierr);
		}
	}
	
	ierr = MatAssemblyBegin(F_permuted,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(F_permuted,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatCopy(F_permuted,F,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	
	ierr = MatDestroy(&F_permuted);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode SaveOutputs(Resolvent_matrices *Res_mat, Directories *dirs)
{

	/*
		Saves the input/output modes and gains in the output directory
	*/
	
	PetscErrorCode        ierr;
	PetscViewer           fd;

	PetscFunctionBeginUser;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Saving onto the disk in progress... \nThe output directory is: ");
	ierr = PetscPrintf(PETSC_COMM_WORLD,"%s%s\n",dirs->root_dir,dirs->output);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"S_hat_clean_V1_RK4");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = MatView(Res_mat->S_hat,fd);CHKERRQ(ierr);
	// ierr = MatView(Res_mat->S_hat,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"V_hat_clean_V1_RK4");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = MatView(Res_mat->V_hat,fd);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"U_hat_clean_V1_RK4");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = MatView(Res_mat->U_hat,fd);CHKERRQ(ierr);	

	ierr = PetscPrintf(PETSC_COMM_WORLD,"DONE :))\n");CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode SaveSnapshots(Vec y, PetscInt i, RSVDt_vars *RSVDt, Mat Y_all)
{

	/*
		Stores the response snapshots at equal time intervals before perfoming DFT
	*/

	PetscErrorCode        ierr;
	PetscInt              pos;
	Vec                   Y_temp;

	// PetscReal norm_real;

	PetscFunctionBeginUser;

	if (i > RSVDt->TS.Nt_transient) {
		pos = i;
		if (PetscFmodReal(pos,RSVDt->TS.time_ratio) == 0) {
			pos  = pos/RSVDt->TS.time_ratio;
			pos  = PetscFmodReal(pos, RSVDt->RSVD.Nw);
			pos  = RSVDt->TS.flg_dir_adj ? pos : RSVDt->RSVD.Nw-pos;
			pos  = PetscFmodReal(pos, RSVDt->RSVD.Nw);
			pos += RSVDt->TS.flg_eff_trans ? 1 : 0;
			ierr = MatDenseGetColumnVecWrite(Y_all,pos,&Y_temp);CHKERRQ(ierr);
			ierr = VecCopy(y,Y_temp);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(Y_all,pos,&Y_temp);CHKERRQ(ierr);
			// ierr = VecNorm(y,NORM_2,&norm_real);CHKERRQ(ierr);
			// ierr = PetscPrintf(PETSC_COMM_WORLD,"norm = %f, pos = %d, i = %d\n",(float)norm_real,(int)pos,(int)i);CHKERRQ(ierr);
		}
	}

	PetscFunctionReturn(0);
	
}

static PetscErrorCode SetupTimeFreqGrid(RSVDt_vars *RSVDt)
{

	/*
		Initializes the time integration variables including time step and length of time integration
	*/

	PetscErrorCode        ierr;
	PetscReal             T_ss,dt_loc,dt_w_comm,Nw;
	PetscInt              Nt_loc,Nt_transient_loc;

	PetscFunctionBeginUser;

	// Nw_Nw_out                = RSVDt->RSVD.Nw * RSVDt->RSVD.Nw_out;
	Nw                       = RSVDt->RSVD.Nw;
	T_ss                     = 2*PETSC_PI/RSVDt->RSVD.w_min;
	dt_w_comm                = T_ss/Nw;
	// dt_w                     = T_ss/RSVDt->RSVD.Nw;
	// dt_w_out                 = T_ss/RSVDt->RSVD.Nw_out;
	dt_loc                   = dt_w_comm/PetscCeilReal(dt_w_comm/RSVDt->TS.dt_max);
	RSVDt->TS.dt             = dt_loc;
	Nt_loc                   = round(T_ss/dt_loc);
	RSVDt->TS.Nt             = Nt_loc;
	Nt_transient_loc         = round(RSVDt->TS.t_transient/dt_loc);
	RSVDt->TS.Nt_transient   = Nt_transient_loc;
	RSVDt->TS.time_ratio     = Nt_loc/RSVDt->RSVD.Nw;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"@ w_min = %g\n@ t_transient = %g\n@ dt = %g\n@ Nt = %d \
							\n@ Nt_transient = %d\n@ time_ratio = %d\n",\
							RSVDt->RSVD.w_min,RSVDt->TS.t_transient, \
							RSVDt->TS.dt,(int)RSVDt->TS.Nt,(int)RSVDt->TS.Nt_transient,\
							(int)RSVDt->TS.time_ratio);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode StoreQ(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt)
{

	/*
		Stores \hat{Y} on \hat{Q} before the final adjoint action
		to recover \hat{U} after the economy SVD
	*/

	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->Q_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVD_mat->Q_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(RSVD_mat->Q_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,RSVDt->RSVD.Nw*RSVDt->RSVD.k);CHKERRQ(ierr); // Nw_eff to be added
	ierr = MatSetUp(RSVD_mat->Q_hat);CHKERRQ(ierr);	
	ierr = MatCopy(RSVD_mat->Y_hat,RSVD_mat->Q_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode TransientRK4(Mat Y_all, Mat trans_norm, PetscInt k_trans, PetscInt mod, PetscInt rseed, \
		LNS_params *LNS_mat, WS_matrices *WS_mat, TS_matrices *TS_mat, RSVDt_vars *RSVDt, PetscBool Y_all_flg, Directories *dirs)
{

	/*
		Solves \dot{q} - Aq = 0 in time starting from a random or given initial condition
	*/

	PetscErrorCode       ierr;
	Mat                  Y_IC_rand,Y_all_k,Y_temp;
	Vec                  y0,norm_vec,y_temp;
	PetscInt             Nt_saved,i,ik,rend,pos;
	PetscReal            norm_real,div_val=1e3,con_val=1e-16;
	PetscScalar          norm;
	PetscBool            flg_div=0,flg_IC;
	PetscRandom          r;
	PetscViewer          fd;
	PetscLogDouble       t1, t2;

	PetscFunctionBeginUser;

	ierr = RSVDt->TS.flg_dir_adj ? PetscPrintf(PETSC_COMM_WORLD,"\n*** Transient run (direct action) in process! ***\n") : \
				PetscPrintf(PETSC_COMM_WORLD,"\n*** Transient run (adjoint action) in process! ***\n");CHKERRQ(ierr);

	/*
		Define TS parameters
	*/

	RSVDt->TS.dt = RSVDt->TS.dt_max;
	rend  = round(RSVDt->TS.t_transient/RSVDt->TS.dt);
	// mod   = mod <= 10 ? mod : RSVDt->TS.time_ratio;
	Nt_saved = PetscFloorReal(rend/mod) + 1;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Total snapshots to be saved = %d\n", (int) Y_all_flg * (int) Nt_saved);CHKERRQ(ierr);

	/*
		Creates the required matrices/vecs
	*/

	if (Y_all_flg) {
		ierr = MatSetType(Y_all,MATDENSE);CHKERRQ(ierr);
		ierr = MatSetSizes(Y_all,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,Nt_saved*k_trans);CHKERRQ(ierr);
		ierr = MatSetUp(Y_all);CHKERRQ(ierr);

		ierr = MatCreate(PETSC_COMM_WORLD,&Y_all_k);CHKERRQ(ierr);
		ierr = MatSetType(Y_all_k,MATDENSE);CHKERRQ(ierr);
		ierr = MatSetSizes(Y_all_k,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,Nt_saved);CHKERRQ(ierr);
		ierr = MatSetUp(Y_all_k);CHKERRQ(ierr);
	}

	ierr = MatSetType(trans_norm,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(trans_norm,PETSC_DECIDE,PETSC_DECIDE,Nt_saved,k_trans);CHKERRQ(ierr);
	ierr = MatSetUp(trans_norm);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&y0);CHKERRQ(ierr);
	ierr = VecSetType(y0,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(y0,PETSC_DECIDE,RSVDt->RSVD.N);CHKERRQ(ierr);	

	if (k_trans > 1) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Random matrix of size %d x %d!\n", (int)RSVDt->RSVD.N, (int)k_trans);CHKERRQ(ierr);
		ierr = PetscRandomCreate(PETSC_COMM_WORLD,&r);CHKERRQ(ierr);
		ierr = PetscRandomSetSeed(r, rseed);CHKERRQ(ierr);
		ierr = PetscRandomSeed(r);CHKERRQ(ierr);		
		ierr = MatCreate(PETSC_COMM_WORLD,&Y_IC_rand);CHKERRQ(ierr);
		ierr = MatSetType(Y_IC_rand,MATDENSE);CHKERRQ(ierr);
		ierr = MatSetSizes(Y_IC_rand,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,k_trans);CHKERRQ(ierr);
		ierr = MatSetUp(Y_IC_rand);CHKERRQ(ierr);
		ierr = MatSetRandom(Y_IC_rand,r);CHKERRQ(ierr);
	} else {	
		ierr = PetscOptionsGetString(NULL,NULL,"-IC",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,&flg_IC);CHKERRQ(ierr);
		if (flg_IC){
			ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,dirs->filename);CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);			
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading a vector!\n");CHKERRQ(ierr);
			ierr = VecLoad(y0,fd);CHKERRQ(ierr);
			if (WS_mat->flg_Win) {
				ierr = VecPointwiseMult(y_temp,WS_mat->W_in_vec,y0);CHKERRQ(ierr);
				ierr = VecDot(y_temp,y0,&norm);CHKERRQ(ierr);
				norm = PetscSqrtReal(PetscRealPart(norm));
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of test vector: %f\n",(double)norm);CHKERRQ(ierr);  
				// ierr = VecScale(y0,1/norm);CHKERRQ(ierr);
			} else {
				// ierr = VecNormalize(y0,NULL);CHKERRQ(ierr);
			}		
		} else {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Random vector!\n");CHKERRQ(ierr);
			ierr = PetscRandomCreate(PETSC_COMM_WORLD,&r);CHKERRQ(ierr);
			ierr = PetscRandomSetSeed(r, rseed);CHKERRQ(ierr);
			ierr = PetscRandomSeed(r);CHKERRQ(ierr);
			ierr = VecSetRandom(y0,r);CHKERRQ(ierr);
			if (WS_mat->flg_Win && WS_mat->flg_diag) {
				ierr = VecPointwiseMult(y_temp,WS_mat->W_in_vec,y0);CHKERRQ(ierr);
				ierr = VecDot(y_temp,y0,&norm);CHKERRQ(ierr);
				norm = PetscSqrtReal(PetscRealPart(norm));
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of test vector: %f\n",(double)norm);CHKERRQ(ierr);  
				ierr = VecScale(y0,1/norm);CHKERRQ(ierr);
			} else if (WS_mat->flg_Win) {
				ierr = MatMult(WS_mat->W_in,y0,y_temp);CHKERRQ(ierr);
				ierr = VecDot(y_temp,y0,&norm);CHKERRQ(ierr);
				norm = PetscSqrtReal(PetscRealPart(norm));
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of test vector: %f\n",(double)norm);CHKERRQ(ierr);  
				ierr = VecScale(y0,1/norm);CHKERRQ(ierr); 
			} else {
				ierr = VecNormalize(y0,NULL);CHKERRQ(ierr);
			}
		}
	}

	ierr = CreateTSMats(TS_mat,RSVDt->RSVD.N);CHKERRQ(ierr);

	/*
		Loops over all test vectors
	*/	

	// ierr = PetscTime(&t1);CHKERRQ(ierr);

	for (ik=0; ik<k_trans; ik++) {

		if (flg_div) break;
		pos = 0;

		/*
			Extracts the random IC for each test vector
		*/

		if (k_trans > 1) {
			ierr = MatDenseGetColumnVecRead(Y_IC_rand,ik,&y_temp);CHKERRQ(ierr);
			ierr = VecCopy(y_temp,y0);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecRead(Y_IC_rand,ik,&y_temp);CHKERRQ(ierr);
		}

		if (WS_mat->flg_Win && WS_mat->flg_diag) {
			ierr = VecPointwiseMult(y_temp,WS_mat->W_in_vec,y0);CHKERRQ(ierr);
			ierr = VecDot(y_temp,y0,&norm);CHKERRQ(ierr);
			norm = PetscSqrtReal(PetscRealPart(norm));
			ierr = VecScale(y0,1/norm);CHKERRQ(ierr);
		} else if (WS_mat->flg_Win) {
			ierr = MatMult(WS_mat->W_in,y0,y_temp);CHKERRQ(ierr);
			ierr = VecDot(y_temp,y0,&norm);CHKERRQ(ierr);
			norm_real = PetscSqrtReal(PetscRealPart(norm));
			ierr = VecScale(y0,1/norm);CHKERRQ(ierr);
		} else{
			ierr = VecNorm(y0,NORM_2,&norm_real);CHKERRQ(ierr);
			norm = norm_real;
			// ierr = VecScale(y0,1/norm);CHKERRQ(ierr);
			ierr = VecNorm(y0,NORM_2,&norm_real);CHKERRQ(ierr);
			norm = norm_real;			
		}

		ierr = VecCreate(PETSC_COMM_WORLD,&norm_vec);CHKERRQ(ierr);
		ierr = VecSetType(norm_vec,VECMPI);CHKERRQ(ierr);
		ierr = VecSetSizes(norm_vec,PETSC_DECIDE,Nt_saved);CHKERRQ(ierr);	
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of test vector after normalization: %f @ t = %d\n",(double)norm,(int)0);CHKERRQ(ierr);           
		ierr = VecSetValues(norm_vec,1,&pos,&norm,INSERT_VALUES);CHKERRQ(ierr);			

		if (Y_all_flg) {
			ierr = MatDenseGetColumnVecWrite(Y_all_k,ik,&y_temp);CHKERRQ(ierr);
			ierr = VecCopy(y0,y_temp);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(Y_all_k,ik,&y_temp);CHKERRQ(ierr);		
		}

		/*
			Performs time stepping
		*/

		ierr = PetscTime(&t1);CHKERRQ(ierr);

		for (i=1; i<=rend; i++) {

			ierr = TSTransRK4(LNS_mat,TS_mat,RSVDt,i,y0);CHKERRQ(ierr);

			if (PetscFmodReal(i,mod) == 0) {
				pos  = i/mod;
				if (Y_all_flg) {
					ierr = MatDenseGetColumnVecWrite(Y_all_k,pos,&y_temp);CHKERRQ(ierr);
					ierr = VecCopy(y0,y_temp);CHKERRQ(ierr);
					ierr = MatDenseRestoreColumnVecWrite(Y_all_k,pos,&y_temp);CHKERRQ(ierr);
				}

				if (WS_mat->flg_Win && WS_mat->flg_diag) {
					ierr = VecPointwiseMult(y_temp,WS_mat->W_in_vec,y0);CHKERRQ(ierr);
					ierr = VecDot(y_temp,y0,&norm);CHKERRQ(ierr);
					norm_real = PetscSqrtReal(PetscRealPart(norm));
				} else if (WS_mat->flg_Win) {
					ierr = MatMult(WS_mat->W_in,y0,y_temp);CHKERRQ(ierr);
					ierr = VecDot(y_temp,y0,&norm);CHKERRQ(ierr);
					norm_real = PetscSqrtReal(PetscRealPart(norm));
				} else{
					ierr = VecNorm(y0,NORM_2,&norm_real);CHKERRQ(ierr);
					norm = norm_real;
				}
				
				if (norm_real > div_val) {
					PetscPrintf(PETSC_COMM_WORLD,"Diverged!\n");CHKERRQ(ierr); 
					break;
					flg_div = 1;
				}

				if (norm_real < con_val) {
					PetscPrintf(PETSC_COMM_WORLD,"Converged!\n");CHKERRQ(ierr); 
					break;
				}
				
				norm = norm_real;
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of test vector: %f @ t = %f\n",(double)norm,i*RSVDt->TS.dt);CHKERRQ(ierr);
				ierr = VecSetValues(norm_vec,1,&pos,&norm,INSERT_VALUES);CHKERRQ(ierr);
			}
		}

		ierr = PetscTime(&t2);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"TS took %f seconds (k = %d)\n", t2-t1, (int)ik);CHKERRQ(ierr);

		ierr = VecAssemblyBegin(norm_vec);CHKERRQ(ierr);
		ierr = VecAssemblyEnd(norm_vec);CHKERRQ(ierr);
		ierr = MatDenseGetColumnVecWrite(trans_norm,ik,&y_temp);CHKERRQ(ierr);
		ierr = VecCopy(norm_vec,y_temp);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(trans_norm,ik,&y_temp);CHKERRQ(ierr);

		if (Y_all_flg) {
			ierr = MatAssemblyBegin(Y_all_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
			ierr = MatAssemblyEnd(Y_all_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
			ierr = MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nt_saved,(ik+1)*Nt_saved,&Y_temp);CHKERRQ(ierr);
			ierr = MatCopy(Y_all_k,Y_temp,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
			ierr = MatDenseRestoreSubMatrix(Y_all,&Y_temp);CHKERRQ(ierr);
		}

		ierr = PetscTime(&t2);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Transient took %f seconds!\n", t2-t1);CHKERRQ(ierr);
	
	}

	/*
		Saves the outputs
	*/	

	if (Y_all_flg && !flg_div) {
		ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"Y_all_transient");CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
		ierr = MatView(Y_all,fd);CHKERRQ(ierr);
	}

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"Y_mod_transient_norms");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = MatView(trans_norm,fd);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"Y_transient_last");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = VecView(y0,fd);CHKERRQ(ierr);	

	/*
		Removes variables from memory and return
	*/

	ierr = VecDestroy(&y0);CHKERRQ(ierr);
	ierr = MatDestroy(&Y_all);CHKERRQ(ierr);
	ierr = MatDestroy(&Y_all_k);CHKERRQ(ierr);
	ierr = MatDestroy(&trans_norm);CHKERRQ(ierr);
	ierr = DestroyTSMats(TS_mat);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}

static PetscErrorCode TransientRunRK4(RSVDt_vars *RSVDt, LNS_params *LNS_mat, WS_matrices *WS_mat, TS_matrices *TS_mat, Directories *dirs)
{

	/*
		Transient simulation using RK4
		This function runs only if -TransientRun_flg is True
	*/
	
	PetscErrorCode        ierr;
	Mat                   Y_all, trans_norm;
	PetscBool             TransientSave_flg, TransientRun_dir;
	PetscInt              k_trans, rseed, mod = 10;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetBool(NULL,NULL,"-TransientSave_flg",&TransientSave_flg,&TransientSave_flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-TransientRun_dir",&TransientRun_dir,&TransientRun_dir);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-k_trans",&k_trans,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-rseed",&rseed,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-mod",&mod,NULL);CHKERRQ(ierr);

	RSVDt->TS.flg_dir_adj = TransientRun_dir;

	ierr = MatCreate(PETSC_COMM_WORLD,&Y_all);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&trans_norm);CHKERRQ(ierr);

	/*
		Time stepping of the homogenous system
	*/

	ierr = TransientRK4(Y_all,trans_norm,k_trans,mod,rseed,LNS_mat,WS_mat,TS_mat,RSVDt,TransientSave_flg,dirs);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode TSActionRK4(RSVD_matrices *RSVD_mat, DFT_matrices *DFT_mat, \
		LNS_params *LNS_mat, TS_matrices *TS_mat, RSVDt_vars *RSVDt, Directories *dirs)
{

	/*
		Time stepping of the LNS equations to obtain either action of direct or adjoint resolvent operator
	*/

	PetscErrorCode       ierr;
	KSP                  ksp;
	Mat                  Y_all,F_temp,Y_all_k,F_hat_k,Q_basis,EQ,M_tilde;
	Vec                  y,y_temp,y1,y2,b;
	PetscInt             N,k,Nw,Nw_eff,i,ik,rend,prg_cnt=0,period_ind,pos;
	// PetscReal            norm;
	PetscLogDouble       t1,t2;

	// PetscReal norm_real;
	// PetscViewer     fd;
	// char            output_dir[PETSC_MAX_PATH_LEN];           /* output root directory */
	// char            file_temp[PETSC_MAX_PATH_LEN];            /* temporary file names */
	// char            filename_temp[64];                        /* temporary output file names */
	
	PetscFunctionBeginUser;

	// strcpy(output_dir, "/scratch/towne_root/towne0/aliii/output/");

	ierr = RSVDt->TS.flg_dir_adj ? PetscPrintf(PETSC_COMM_WORLD,"\n*** Direct action in process! ***\n") : \
				PetscPrintf(PETSC_COMM_WORLD,"\n*** Adjoint action in process! ***\n");CHKERRQ(ierr);

	/*
		Defines the time-stepping parameters
	*/

	k     = RSVDt->RSVD.k;
	rend  = RSVDt->TS.Nt_transient + RSVDt->TS.Nt;
	rend += RSVDt->TS.flg_eff_trans ? RSVDt->TS.time_ratio : 0;

	/*
		Gets the size of matrices
	*/

	ierr = MatGetSize(RSVD_mat->Y_hat,&N,&Nw);CHKERRQ(ierr);
	Nw  /= k;
	
	/*
		Uses the solution from the previous action as forcing for the current action
	*/

	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->F_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVD_mat->F_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(RSVD_mat->F_hat,PETSC_DECIDE,PETSC_DECIDE,N,Nw*k);CHKERRQ(ierr);
	ierr = MatSetUp(RSVD_mat->F_hat);CHKERRQ(ierr);	
	ierr = MatCopy(RSVD_mat->Y_hat,RSVD_mat->F_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = MatDestroy(&RSVD_mat->Y_hat);CHKERRQ(ierr);

	/*
		Creates all required matrices (including solution, forcing, RHS, norm, etc.)
	*/

	ierr = VecCreate(PETSC_COMM_WORLD,&y);CHKERRQ(ierr);
	ierr = VecSetType(y,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(y,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = CreateTSMats(TS_mat,N);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&Y_all_k);CHKERRQ(ierr);
	ierr = MatSetType(Y_all_k,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Y_all_k,PETSC_DECIDE,PETSC_DECIDE,N,Nw);CHKERRQ(ierr);
	ierr = MatSetUp(Y_all_k);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&Y_all);CHKERRQ(ierr);
	ierr = MatSetType(Y_all,MATDENSE);CHKERRQ(ierr);
	ierr = !RSVDt->TS.flg_eff_trans ? MatSetSizes(Y_all,PETSC_DECIDE,PETSC_DECIDE,N,Nw*k) : \
						MatSetSizes(Y_all,PETSC_DECIDE,PETSC_DECIDE,N,(Nw+1)*k);CHKERRQ(ierr);
	ierr = MatSetUp(Y_all);CHKERRQ(ierr);

	/*
		Permutes the forcing matrix
	*/

	ierr = PermuteMat(RSVD_mat->F_hat, RSVDt);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&F_hat_k);CHKERRQ(ierr);
	ierr = MatSetType(F_hat_k,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(F_hat_k,PETSC_DECIDE,PETSC_DECIDE,N,Nw);CHKERRQ(ierr);
	ierr = MatSetUp(F_hat_k);CHKERRQ(ierr);

	/*
		Extracts the forcing submatrix for each mode and loop-over all modes
	*/	

	for (ik=0; ik<k; ik++) {

		ierr = VecZeroEntries(y);CHKERRQ(ierr);

		ierr = MatDenseGetSubMatrix(RSVD_mat->F_hat,PETSC_DECIDE,PETSC_DECIDE,ik*Nw,(ik+1)*Nw,&F_temp);CHKERRQ(ierr);
		ierr = MatCopy(F_temp,F_hat_k,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		ierr = MatAssemblyBegin(F_hat_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(F_hat_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(RSVD_mat->F_hat,&F_temp);CHKERRQ(ierr);

		/*
			Creates the initial forcing
		*/

		ierr = CreateForcingOnFly(F_hat_k,DFT_mat,0,TS_mat->F3);CHKERRQ(ierr);

		/*
			Measures the time-stepping wall-time
		*/

		ierr = PetscTime(&t1);CHKERRQ(ierr);

		/*
			Performs time stepping
		*/

		for (i=1; i<=rend; i++) {

			ierr = TSRK4(F_hat_k,DFT_mat,LNS_mat,TS_mat,RSVDt,i,y);CHKERRQ(ierr);

			ierr = SaveSnapshots(y,i,RSVDt,Y_all_k);CHKERRQ(ierr);

			// ierr = VecNorm(y,NORM_2,&norm_real);CHKERRQ(ierr);
			// ierr = PetscPrintf(PETSC_COMM_WORLD,"norm is: %f% \n",(float)norm_real);CHKERRQ(ierr);

			// ierr = DisplayProgress(y,i,RSVDt,Y_all_k);CHKERRQ(ierr);

			if ((double)i/rend >= 0.01*prg_cnt){
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Progress percentage: %d% \n",(int)prg_cnt);CHKERRQ(ierr);
				ierr = PetscTime(&t2);CHKERRQ(ierr);
				ierr = PetscPrintf(PETSC_COMM_WORLD,"TS took %f seconds (wall-time)\n", t2-t1);CHKERRQ(ierr);		
				prg_cnt += 10;
			}

		}

		ierr = PetscTime(&t2);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"TS took %f seconds (k = %d)\n", t2-t1, (int)ik);CHKERRQ(ierr);

		ierr = MatAssemblyBegin(Y_all_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(Y_all_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nw,(ik+1)*Nw,&F_temp);CHKERRQ(ierr);
		ierr = MatCopy(Y_all_k,F_temp,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Y_all,&F_temp);CHKERRQ(ierr);

	}

	// snprintf(filename_temp, sizeof(char) * 64, "Y_all_k_clean");
	// strcpy(file_temp, "");
	// strcat(file_temp,output_dir);
	// strcat(file_temp,filename_temp);
	// ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file_temp,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	// ierr = MatView(Y_all_k,fd);CHKERRQ(ierr);

	// // snprintf(filename_temp, sizeof(char) * 64, "Y_all_clean");
	// // strcpy(file_temp, "");
	// // strcat(file_temp,output_dir);
	// // strcat(file_temp,filename_temp);
	// // ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file_temp,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	// // ierr = MatView(Y_all,fd);CHKERRQ(ierr);

	/*
		Removes the forcing input as it is no longer needed after the repsonse is obtained
	*/

	ierr = VecDestroy(&y);CHKERRQ(ierr);
	ierr = MatDestroy(&RSVD_mat->F_hat);CHKERRQ(ierr);
	ierr = DestroyTSMats(TS_mat);CHKERRQ(ierr);

	/*
		Takes the response to the frequency domain
		Efficient transient removal is perfomed before DFT if desired
	*/

	ierr = DFT(Y_all,DFT_mat,RSVD_mat,RSVDt);CHKERRQ(ierr);

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!     Add this     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	// Efficient transient removal if desired

	// if (!flg_eff_trans) {
	// 	ierr = ManualDFT(Y_all,dft,time_ratio,B_hat,flg_real_A);CHKERRQ(ierr);
	// } else {
	// 	ierr = SubspaceTransientRemoval(A,Y_all,B_hat,w_min,dt,k,time_ratio,dft,flg_real_A,flg_flip);CHKERRQ(ierr);
	// }
	
	ierr = RSVDt->TS.flg_dir_adj ? PetscPrintf(PETSC_COMM_WORLD,"Direct action DONE!\n") : \
				PetscPrintf(PETSC_COMM_WORLD,"Adjoint action DONE!\n");CHKERRQ(ierr);

	PetscFunctionReturn(0);

}

static PetscErrorCode TSRK4(Mat F_hat, DFT_matrices *DFT_mat, LNS_params *LNS_mat, \
					TS_matrices *TS_mat, RSVDt_vars *RSVDt, PetscInt i, Vec y)
{

	/*
		Performs one RK4 iteration after creating forcing terms on fly
	*/
	
	PetscErrorCode        ierr;
	PetscInt              jt_cyc;

	PetscFunctionBeginUser;

	ierr = RSVDt->TS.flg_dir_adj ? MatMultAdd(LNS_mat->A,y,TS_mat->F3,TS_mat->k1) : \
				MatMultHermitianTransposeAdd(LNS_mat->A,y,TS_mat->F3,TS_mat->k1);CHKERRQ(ierr);

	jt_cyc = RSVDt->TS.flg_dir_adj ? PetscFmodReal(2*(i-1)+1,2*RSVDt->TS.Nt) : \
					PetscFmodReal(1e6*RSVDt->TS.Nt-1-(2*(i-1)), 2*RSVDt->TS.Nt);
	ierr = CreateForcingOnFly(F_hat,DFT_mat,jt_cyc,TS_mat->F2);CHKERRQ(ierr);

	jt_cyc = RSVDt->TS.flg_dir_adj ? PetscFmodReal(2*(i-1)+2,2*RSVDt->TS.Nt) : \
					PetscFmodReal(1e6*RSVDt->TS.Nt-1-(2*(i-1)+1), 2*RSVDt->TS.Nt);
	ierr = CreateForcingOnFly(F_hat,DFT_mat,jt_cyc,TS_mat->F3);CHKERRQ(ierr);

	if (RSVDt->TS.flg_dir_adj) {
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k1,y);CHKERRQ(ierr);
		ierr = MatMultAdd(LNS_mat->A,TS_mat->y_temp,TS_mat->F2,TS_mat->k2);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k2,y);CHKERRQ(ierr);
		ierr = MatMultAdd(LNS_mat->A,TS_mat->y_temp,TS_mat->F2,TS_mat->k3);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt,TS_mat->k3,y);CHKERRQ(ierr);
		ierr = MatMultAdd(LNS_mat->A,TS_mat->y_temp,TS_mat->F3,TS_mat->k4);CHKERRQ(ierr);
	} else {
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k1,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTransposeAdd(LNS_mat->A,TS_mat->y_temp,TS_mat->F2,TS_mat->k2);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k2,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTransposeAdd(LNS_mat->A,TS_mat->y_temp,TS_mat->F2,TS_mat->k3);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt,TS_mat->k3,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTransposeAdd(LNS_mat->A,TS_mat->y_temp,TS_mat->F3,TS_mat->k4);CHKERRQ(ierr);	
	}
	ierr = VecAXPY(y,RSVDt->TS.dt/6,TS_mat->k1);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6*2,TS_mat->k2);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6*2,TS_mat->k3);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6,TS_mat->k4);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode TSTimeDomainRunRK4(RSVD_matrices *RSVD_mat, DFT_matrices *DFT_mat, LNS_params *LNS_mat, WS_matrices *WS_mat, \
							TS_matrices *TS_mat, RSVDt_vars *RSVDt, Directories *dirs)
{
	
	PetscErrorCode        ierr;
	Mat                   Y_all, trans_norm;
	PetscBool             TSTimeDomainRun_flg, TSTimeDomainRun_dir;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetBool(NULL,NULL,"-TSTimeDomainRun_dir",&TSTimeDomainRun_dir,&TSTimeDomainRun_dir);CHKERRQ(ierr);

	RSVDt->TS.flg_dir_adj = TSTimeDomainRun_dir;

	if (RSVDt->TS.flg_dir_adj) {	
		if (WS_mat->flg_Wout) {
			ierr = WS_mat->flg_diag ? MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->W_out_vec,NULL) : \
				MatMatMult(WS_mat->W_out,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);

		}
		if (WS_mat->flg_B) ierr = MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->B_vec,NULL);CHKERRQ(ierr);
	} else {
		if (WS_mat->flg_Win) {
			ierr = WS_mat->flg_diag ? MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->W_in_vec,NULL) : \
				MatMatMult(WS_mat->W_in,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);
		}
		if (WS_mat->flg_C) ierr = MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->C_vec,NULL);CHKERRQ(ierr);
	}

	ierr = TSActionRK4(RSVD_mat,DFT_mat,LNS_mat,TS_mat,RSVDt,dirs);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

static PetscErrorCode TSTransRK4(LNS_params *LNS_mat, \
					TS_matrices *TS_mat, RSVDt_vars *RSVDt, PetscInt i, Vec y)
{

	/*
		Performs one RK4 iteration assuming zero forcing 
	*/
	
	PetscErrorCode       ierr;
	PetscInt             jt_cyc;
	PetscLogDouble       t1, t2;

	PetscFunctionBeginUser;

	if (RSVDt->TS.flg_dir_adj) {
		ierr = MatMultAdd(LNS_mat->A,y,TS_mat->F1,TS_mat->k1);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k1,y);CHKERRQ(ierr);
		ierr = MatMult(LNS_mat->A,TS_mat->y_temp,TS_mat->k2);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k2,y);CHKERRQ(ierr);
		ierr = MatMult(LNS_mat->A,TS_mat->y_temp,TS_mat->k3);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt,TS_mat->k3,y);CHKERRQ(ierr);
		ierr = MatMult(LNS_mat->A,TS_mat->y_temp,TS_mat->k4);CHKERRQ(ierr);
	} else{
		ierr = MatMultHermitianTransposeAdd(LNS_mat->A,y,TS_mat->F1,TS_mat->k1);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k1,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTranspose(LNS_mat->A,TS_mat->y_temp,TS_mat->k2);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k2,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTranspose(LNS_mat->A,TS_mat->y_temp,TS_mat->k3);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt,TS_mat->k3,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTranspose(LNS_mat->A,TS_mat->y_temp,TS_mat->k4);CHKERRQ(ierr);	
	}

	ierr = VecAXPY(y,RSVDt->TS.dt/6,TS_mat->k1);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6*2,TS_mat->k2);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6*2,TS_mat->k3);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6,TS_mat->k4);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
	
}




