
/* 	
	Author: Ali Farghadan, 2024
	Produced at the Univeristy of Michigan, Towne Lab
	Reference paper: Scalable resolvent analysis for three-dimensional flows, JCP, 2024
*/

/*

	List of inputs ** Description ************************************ Format

	RSVD-delta-t variables:
	A                 linear (or LNS) operator                         matrix
	B                 input matrix as defined in the reference paper   matrix
	C                 output matrix as defined in the reference paper  matrix
	W_q_sqrt          W_q^(1/2) as defined in the reference paper      matrix
	W_q_sqrt_inv      W_q^(-1/2) as defined in the reference paper     matrix
	W_f_sqrt_inv      W_f^(-1/2) as defined in the reference paper     matrix
	k                 number of test vectors                           integer
	q                 number of power iterations                       integer
	w                 base frequency                                   real 
	TwoPI             base frequency multiplies by 2*pi if true        boolean 
	Nw                number of frequencies to resolve                 integer
	dt                time step                                        real
	RootDir           root directory                                   string
	ResultsDir        results directory (RootDir/ResultsDir)           string
	TransientLength   transient length                                 real
	beta              beta value for discounting (A <-- A - beta I)    real > 0
	TransientRemoval  performs Galerkin transient removal if true      boolean
	RandSeed          seeding random number                            integer
	DiscFlg           applies discounting for unstable linear systems  boolean
	InputForcingFlg   starts from a specified forcing input            boolean
	InputMatrixFlg    applies input matrix                             boolean 
	OutputMatrixFlg   applies output weight matrix                     boolean 
	InputWeightFlg    applies input weight matrix                      boolean 
	OutputWeightFlg   applies output weight matrix                     boolean 
	Display           display options                                  integer
	    case 1) Display = 0: nothing
	    case 2) Display = 1: problem information + elapsed time of each test vector (and total elapsed time) + estimated remaining time
	    case 3) Display = 2: "Display = 1" information + progress percentage of the first test vector (every 10 percent)
	SaveResultsOpt    saving resolvent modes options                   integer
	    case 1) SaveResultsOpt = 1: saves resolvent modes as k  matrices of size N x Nw
	    case 2) SaveResultsOpt = 2: saves resolvent modes as Nw matrices of size N x k

	Transient simulation variables:
	TransRun          runs transient simulation and exits              boolean	
	TransRemovalEst   estimates the transient error if true            boolean
	TransSave         saves the transient outputs if true              boolean
	TransPeriods      number of periods to integrate                   integer
	TransSaveMod      saves the snapshots every "TransSaveMod" number  integer
	TransDivVal       divergence value                                 real
	TransConVal       convergence value                                real
	TransICFlg        starts transient simulation from a specified IC  boolean
	TransICDir        initial vector directory (RootDir/TransICDir)    boolean

	List of outputs * Description ************************************ Format

	U                 response resolvent modes                         matrix
	V                 forcing resolvent modes                          matrix
	Sigma             resolvent gains                                  matrix

*/

/* 	
	List of input libraries and functions
*/

#include <slepcsvd.h>
#include <Variables.h>
#include <PreProcessing.h>
#include <TransientRunRK4.h>
#include <DirectActionRK4.h>
#include <PowerIterationRK4.h>
#include <StoreU.h>
#include <AdjointActionRK4.h>
#include <SVDAllFreqs.h>
#include <SVDAllFreqsBeforeAdjoint.h>
#include <SaveResults.h>

/* 	
	Beginning of the simulation
*/

int main(int argc,char **args)
{

	/* 	
		Defines variable types
	*/

	PetscErrorCode        ierr;                         /* Petsc error code */
	Directories           dirs;                         /* I/O directories */
	RSVDt_vars            RSVDt;                        /* RSVDt variables */
	TransRun_vars         TR_vars;                      /* transient run variables */
	LNS_vars              LNS_mat;                      /* LNS matrix */
	DFT_matrices          DFT_mat;                      /* DFT and inverse DFT matrices */
	TS_removal_matrices   TSR;                          /* transient removal matrices */
	Weight_matrices       Weight_mat;                   /* weight and input/output matrices */
	RSVD_matrices         RSVD_mat;                     /* RSVD matrices */
	Resolvent_matrices    Res_mat;                      /* resolvent modes and gains */

	/*
		Initializes the SLEPc
	*/

	ierr = SlepcInitialize(&argc,&args,(char*)0,NULL); if (ierr) return ierr;

	/*
		Reads user inputs and create required matrices before running the algorithm
	*/
	
	ierr = PreProcessing(&RSVDt, &Weight_mat, &LNS_mat, &RSVD_mat, &Res_mat, &TR_vars, &DFT_mat, &dirs);CHKERRQ(ierr);

	/*
		Transient simulation (if desired -- run and exit)
	*/
	
	if (TR_vars.TransRun) {
		ierr = TransientRunRK4(&TR_vars, &RSVDt, &LNS_mat, &DFT_mat, &dirs);CHKERRQ(ierr);
		ierr = PetscOptionsClear(NULL);CHKERRQ(ierr);
		ierr = SlepcFinalize();
		return ierr;
	}

	/*************************************************************************
		****************  RSVD - $\Delta t$ algorithm  *******************
		****************    for resolvent analysis     *******************
	**************************************************************************/

	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n********************************************\n"
			"*************** RSVD-\\Delta t **************\n********************************************\n");CHKERRQ(ierr);

	ierr = DirectActionRK4(&RSVD_mat, &RSVDt, &LNS_mat, &DFT_mat, &Weight_mat, &dirs, &TSR);CHKERRQ(ierr);

	ierr = PowerIterationRK4(&RSVD_mat, &RSVDt, &LNS_mat, &DFT_mat, &Weight_mat, &dirs, &TSR);CHKERRQ(ierr);

	ierr = SVDAllFreqsBeforeAdjoint(&RSVD_mat, &RSVDt);CHKERRQ(ierr);

	ierr = StoreU(&RSVD_mat, &Res_mat);CHKERRQ(ierr);

	ierr = AdjointActionRK4(&RSVD_mat, &RSVDt, &LNS_mat, &DFT_mat, &Weight_mat, &dirs, &TSR);CHKERRQ(ierr);

	ierr = SVDAllFreqs(&RSVD_mat, &RSVDt, &Weight_mat, &Res_mat);CHKERRQ(ierr);

	ierr = SaveResults(&Res_mat, &RSVDt, &dirs);CHKERRQ(ierr);
	
	ierr = PetscOptionsClear(NULL);CHKERRQ(ierr);
	ierr = SlepcFinalize();
	return ierr;

}


/* 	
	The end!
*/


