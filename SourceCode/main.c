
/* 	
	Author: Ali Farghadan, Sep 2024
	Produced at the Univeristy of Michigan, Towne Lab
	Reference: Scalable resolvent analysis for three-dimensional flows, JCP, 2024
*/

/*

	Inputs ********** Description ************************************ Format

	RSVD-delta-t variables:
	A                 linear (or LNS) operator                         matrix (sparse N x N)
	B                 input matrix as defined in the reference paper   matrix
	C                 output matrix as defined in the reference paper  matrix
	W_q_sqrt          W_q^(1/2) as defined in the reference paper      matrix
	W_f_sqrt_inv      W_f^(-1/2) as defined in the reference paper     matrix
	k                 number of test vectors                           integer
	q                 number of power iterations                       integer
	w                 base frequency                                   real 
	TwoPI             base frequency multiplies by 2*pi if true        boolean 
	Nw                number of input/output frequencies to resolve    integer
	dt                time step                                        real
	RootDir           root directory                                   string
	ResultsDir        results directory (RootDir/ResultsDir)           string
	TransientLength   transient length                                 real
	DiscFlg           applies discounting for unstable linear systems  boolean
	beta              beta value for discounting (A <-- A - beta I)    real > 0
	TransientRemoval  performs Galerkin transient removal if true      boolean
	RandSeed          seeding random number                            integer
	Display           display options                                  integer
	 case 1) Display = 0: nothing
	 case 2) Display = 1: problem information + elapsed time of each test vector (and total elapsed time) + estimated remaining time
	 case 3) Display = 2: "Display = 1" information + progress percentage of the first test vector (every 10 percent)

	SaveResultsOpt    saving resolvent modes options                   integer
	 case 1) SaveResultsOpt = 1: saves resolvent modes as k  matrices of size N x Nw
	 case 2) SaveResultsOpt = 2: saves resolvent modes as Nw matrices of size N x k

	Transient simulation variables:
	TransRun          run transient simulation and exit                boolean	
	TransRemovalEst   estimates the transient error if true            boolean
	TransSave         saves the transient outputs if true              boolean
	TransPeriods      number of periods to integrate                   integer
	TransSaveMod      saves the snapshots every "TransSaveMod" number  integer
	TransDivVal       divergence value                                 real
	TransConVal       convergence value                                real

	Outputs ********* Description ************************************ Format

	U                 response resolvent modes                         matrix (k matrices of size N x Nw)
	V                 forcing resolvent modes                          matrix (k matrices of size N x Nw)
	Sigma             resolvent gains                                  matrix (one matrix of size k x Nw)

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
#include <QRAllFreqs.h>
#include <StoreQ.h>
#include <AdjointActionRK4.h>
#include <SVDAllFreqs.h>
#include <SaveResults.h>

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
	WS_matrices           WS_mat;                       /* weight and input/output matrices */
	RSVD_matrices         RSVD_mat;                     /* RSVD matrices */
	Resolvent_matrices    Res_mat;                      /* resolvent modes and gains */

	/*
		Initializes the SLEPc
	*/

	ierr = SlepcInitialize(&argc,&args,(char*)0,NULL); if (ierr) return ierr;

	/*
		Reads user inputs and create required matrices before running the algorithm
	*/
	
	ierr = PreProcessing(&RSVDt, &WS_mat, &LNS_mat, &RSVD_mat, &Res_mat, &TR_vars, &DFT_mat, &dirs);CHKERRQ(ierr);

	/*
		Transient simulation (if desired -- run and exit)
	*/
	
	if (TR_vars.TransRun) {
		ierr = TransientRunRK4(&TR_vars, &RSVDt, &LNS_mat, &DFT_mat, &dirs);CHKERRQ(ierr);
		ierr = SlepcFinalize();
		return ierr;
	}

	/*************************************************************************
		****************  RSVD - $\Delta t$ algorithm  *******************
		****************    for resolvent analysis     *******************
	**************************************************************************/

	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n********************************************\n"
			"*************** RSVD-\\Delta t **************\n********************************************\n");CHKERRQ(ierr);

	ierr = DirectActionRK4(&RSVD_mat, &RSVDt, &LNS_mat, &DFT_mat, &WS_mat, &dirs, &TSR);CHKERRQ(ierr);

	ierr = PowerIterationRK4(&RSVD_mat, &RSVDt, &LNS_mat, &DFT_mat, &WS_mat, &dirs, &TSR);CHKERRQ(ierr);

	ierr = QRAllFreqs(&RSVD_mat, &RSVDt);CHKERRQ(ierr);

	ierr = StoreQ(&RSVD_mat, &RSVDt);CHKERRQ(ierr);

	ierr = AdjointActionRK4(&RSVD_mat, &RSVDt, &LNS_mat, &DFT_mat, &WS_mat, &dirs, &TSR);CHKERRQ(ierr);

	ierr = SVDAllFreqs(&RSVD_mat, &RSVDt, &WS_mat, &Res_mat);CHKERRQ(ierr);

	ierr = SaveResults(&Res_mat, &RSVDt, &dirs);CHKERRQ(ierr);

	ierr = SlepcFinalize();
	return ierr;

}





