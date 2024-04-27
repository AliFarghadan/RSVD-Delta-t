
/* 	
	Author: Ali Farghadan, 2020-2024
	Produced at the Univeristy of Michigan, Towne Lab
*/

/*

	Inputs ****** Description ************************************ Format

	A             LNS operator                                     matrix (sparse N x N)
	k             number of test vectors                           integer
	q             number of power iterations                       integer
	w             base frequency                                   real 
	twopi_flg     base frequency multiplies by 2*pi                boolean 
	Nw            number of input/output frequencies to resolve    integer
	real_A        real-valued if true, otherwise complex           boolean
	dt            time step                                        real
	TransLen      transient length                                 real
	Discounting   applies discounting for unstable linear systems  boolean
	Beta          Beta value for discounting (A_dis = A - Beta I)  real > 0
	TransRemoval  efficient transient removal                      boolean
	TransRun      run transient simulation and exit                boolean
	root_dir      root directory                                   string
	input_dir     root_dir/input  (default is root dir)            string
	results_dir   root_dir/output (default is root dir)            string
	display       display options                                  integer
		case 1) display = 0: nothing
		case 2) display = 1: problem information + total elapsed time of each action
		case 3) display = 2: display 1 + elapsed time of each test vector + estimated remaining time
		case 4) display = 3: display 2 + progress percentage of the first test vector

	Outputs ****** Description ************************************ Format

	U              response resolvent modes                         matrix (NkNw)
	V              forcing resolvent modes                          matrix (NkNw)
	Sigma          resolvent gains                                  matrix (kNw)

*/

#include <slepcsvd.h>
#include "Variables.h"
#include "PreProcessing.h"
#include "TransientRunRK4.h"
#include "ForwardActionRK4.h"
#include "PowerIterationRK4.h"
#include "ReducedQR.h"
#include "StoreQ.h"
#include "AdjointActionRK4.h"
#include "ReducedSVD.h"
#include "SaveResults.h"

int main(int argc,char **args)
{
	PetscErrorCode        ierr;                         /* Petsc error code */
	Directories           dirs;                         /* I/O directories */
	RSVDt_vars            RSVDt;                        /* RSVDt variables */
	TransientRun_vars     TR_vars;                      /* transient run variables */
	LNS_vars              LNS_mat;                      /* LNS matrices */
	DFT_matrices          DFT_mat;                      /* DFT and inverse DFT matrices */
	TS_removal_matrices   TSR;                          /* transient removal matrices */
	WS_matrices           WS_mat;                       /* weight and input/output matrices */
	RSVD_matrices         RSVD_mat;                     /* RSVD matrices */
	Resolvent_matrices    Res_mat;                      /* resolvent modes and gains */

	ierr = SlepcInitialize(&argc,&args,(char*)0,NULL); if (ierr) return ierr;

	/*
		Reads user inputs and create required matrices before running the algorithm
	*/
	
	ierr = PreProcessing(&RSVDt, &WS_mat, &LNS_mat, &RSVD_mat, &Res_mat, &TR_vars, &DFT_mat, &dirs);CHKERRQ(ierr);

	/*
		Transient simulation (if desired -- run and exit)
	*/
	
	if (TR_vars.TransientRun) {
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

	ierr = ForwardActionRK4(&RSVD_mat, &RSVDt, &LNS_mat, &DFT_mat, &WS_mat, &dirs, &TSR);CHKERRQ(ierr);

	ierr = PowerIterationRK4(&RSVD_mat, &RSVDt, &LNS_mat, &DFT_mat, &WS_mat, &dirs, &TSR);CHKERRQ(ierr);

	ierr = ReducedQR(&RSVD_mat, &RSVDt);CHKERRQ(ierr);

	ierr = StoreQ(&RSVD_mat, &RSVDt);CHKERRQ(ierr);

	ierr = AdjointActionRK4(&RSVD_mat, &RSVDt, &LNS_mat, &DFT_mat, &WS_mat, &dirs, &TSR);CHKERRQ(ierr);

	ierr = ReducedSVD(&RSVD_mat, &RSVDt, &WS_mat, &Res_mat);CHKERRQ(ierr);

	ierr = SaveResults(&Res_mat, &RSVDt, &dirs);CHKERRQ(ierr);

	ierr = SlepcFinalize();
	return ierr;

}





