
/* 	
	Author: Ali Farghadan, 2024
	Produced at the Univeristy of Michigan, Towne Lab
	Reference paper: Scalable resolvent analysis for three-dimensional flows, JCP, 2024
*/

/*

	List of inputs *** Description ************************************ Format

	RSVD-delta-t variables:
	A                  linear (or LNS) operator                         matrix
	B                  input matrix as defined in the reference paper   matrix
	C                  output matrix as defined in the reference paper  matrix
	W_q_sqrt           W_q^(1/2) as defined in the reference paper      matrix
	W_q_sqrt_inv       W_q^(-1/2) as defined in the reference paper     matrix
	W_f_sqrt_inv       W_f^(-1/2) as defined in the reference paper     matrix
	k                  number of test vectors                           integer
	q                  number of power iterations                       integer
	w                  base frequency                                   real 
	TwoPI              base frequency multiplies by 2*pi if true        boolean 
	Nw                 number of frequencies to resolve                 integer
	dt                 time step                                        real
	RootDir            root directory                                   string
	ResultsDir         results directory (RootDir/ResultsDir)           string
	TransientLength    transient length                                 real
	beta               beta value for discounting (A <-- A - beta I)    real > 0
	TransientRemoval   performs Galerkin transient removal if true      boolean
	RandSeed           seeding random number                            integer
	DiscFlg            applies discounting for unstable linear systems  boolean
	InputForcingFlg    starts from a specified forcing input            boolean
	InputMatrixFlg     applies input matrix                             boolean 
	OutputMatrixFlg    applies output matrix                            boolean 
	InputWeightFlg     applies input weight matrix                      boolean
	InvInputWeightFlg  applies inverse input weight matrix              boolean 
	InvOutputWeightFlg applies inverse output weight matrix             boolean 
	Display            display options                                  integer
	    case 1) Display = 0: nothing
	    case 2) Display = 1: problem information + elapsed time of each test vector (and total elapsed time) + estimated remaining time
	    case 3) Display = 2: "Display = 1" information + progress percentage of the first test vector (every 10 percent)
	SaveResultsOpt     saving resolvent modes options                   integer
	    case 1) SaveResultsOpt = 1: saves resolvent modes as k  matrices of size N x Nw
	    case 2) SaveResultsOpt = 2: saves resolvent modes as Nw matrices of size N x k

	Transient simulation variables:
	TransRun           runs transient simulation and exits              boolean	
	TransRemovalEst    estimates the transient error if true            boolean
	TransSave          saves the transient outputs if true              boolean
	TransPeriods       number of periods to integrate                   integer
	TransSaveMod       saves the snapshots every "TransSaveMod" number  integer
	TransDivVal        divergence value                                 real
	TransConVal        convergence value                                real
	TransICFlg         starts transient simulation from a specified IC  boolean
	TransICDir         initial vector directory (RootDir/TransICDir)    boolean

	List of outputs ** Description ************************************ Format

	U                  response resolvent modes                         matrix
	V                  forcing resolvent modes                          matrix
	Sigma              resolvent gains                                  matrix

*/

/* 	
	List of input libraries and functions
*/

#include <slepcsys.h>
#include <Variables.h>
#include <PreProcessing.h>
#include <TransientRunRK4.h>
#include <DirectActionRK4.h>
#include <PowerIterationRK4.h>
#include <AdjointActionRK4.h>
#include <SVDAllFreqs4Response.h>
#include <SVDAllFreqs4Forcing.h>

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
	TransRun_vars         TR;                           /* transient run variables */
	LNS_vars              LNS;                          /* LNS matrix */
	DFT_matrices          DFT;                          /* DFT and inverse DFT matrices */
	TS_removal_matrices   TSR;                          /* transient removal matrices */
	Weight_matrices       Weight;                       /* weight and input/output matrices */
	RSVD_matrices         RSVD;                         /* RSVD matrices */
	Resolvent_matrices    Res;                          /* resolvent modes and gains */
	PetscLogDouble        t1, t2;                       /* time measurement variables for simulation elapsed time */

	/*
		Initializes the SLEPc
	*/

	ierr = SlepcInitialize(&argc,&args,(char*)0,NULL); if (ierr) return ierr;
	ierr = PetscTime(&t1);CHKERRQ(ierr);

	/*
		Reads user inputs and create required matrices before running the algorithm
	*/
	
	ierr = PreProcessing(&RSVDt, &Weight, &LNS, &RSVD, &TR, &DFT, &dirs);CHKERRQ(ierr);

	/*
		Transient simulation (if desired -- run and exit)
	*/
	
	if (TR.TransRun) {
		ierr = TransientRunRK4(&TR, &RSVDt, &LNS, &DFT, &dirs);CHKERRQ(ierr);
		ierr = PetscOptionsClear(NULL);CHKERRQ(ierr);
		ierr = SlepcFinalize();
		return ierr;
	}

	/*************************************************************************
		****************  RSVD - $\Delta t$ algorithm  *******************
		****************    for resolvent analysis     *******************
	**************************************************************************/

	if (RSVDt.Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n********************************************\n"
			"*************** RSVD-\\Delta t **************\n********************************************\n");CHKERRQ(ierr);

	ierr = DirectActionRK4(&RSVD, &RSVDt, &LNS, &DFT, &Weight, &dirs, &TSR);CHKERRQ(ierr);

	ierr = PowerIterationRK4(&RSVD, &RSVDt, &LNS, &DFT, &Weight, &dirs, &TSR);CHKERRQ(ierr);

	ierr = SVDAllFreqs4Response(&RSVD, &RSVDt, &Weight, &Res, &dirs);CHKERRQ(ierr);

	ierr = AdjointActionRK4(&RSVD, &RSVDt, &LNS, &DFT, &Weight, &dirs, &TSR);CHKERRQ(ierr);

	ierr = SVDAllFreqs4Forcing(&RSVD, &RSVDt, &Weight, &Res, &dirs);CHKERRQ(ierr);

	/*
		Prints out the elapsed time and exits
	*/

	ierr          = PetscTime(&t2);CHKERRQ(ierr);
	PetscInt hh   = (t2-t1)/3600;
	PetscInt mm   = (t2-t1-3600*hh)/60;
	PetscInt ss   = t2-t1-3600*hh-mm*60;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Done :))\n\nEntire simulation elapsed time = %02d:%02d:%02d\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);
	
	ierr = PetscOptionsClear(NULL);CHKERRQ(ierr);
	ierr = SlepcFinalize();
	return ierr;

}


/* 	
	The end!
*/


