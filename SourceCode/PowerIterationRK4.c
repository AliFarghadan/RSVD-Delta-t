
#include <petscmat.h>
#include "Variables.h"
#include "ReducedQR.h"
#include "AdjointActionRK4.h"
#include "ForwardActionRK4.h"

PetscErrorCode PowerIterationRK4(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, LNS_vars *LNS_mat, \
			DFT_matrices *DFT_mat, WS_matrices *WS_mat, Directories *dirs, TS_removal_matrices *TSR)
{
	/*
		Performs power iteration for q times 
	*/

	PetscErrorCode        ierr;
	PetscInt              iq;

	PetscFunctionBeginUser;

	for (iq=0; iq<RSVDt->RSVD.q; iq++) {

		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n******** Inside power iteration, %d/%d *******\n",(int)iq+1,(int)RSVDt->RSVD.q);CHKERRQ(ierr);

		ierr = ReducedQR(RSVD_mat, RSVDt);CHKERRQ(ierr);
		ierr = AdjointActionRK4(RSVD_mat, RSVDt, LNS_mat, DFT_mat, WS_mat, dirs, TSR);CHKERRQ(ierr);
		ierr = ReducedQR(RSVD_mat, RSVDt);CHKERRQ(ierr);
		ierr = ForwardActionRK4(RSVD_mat, RSVDt, LNS_mat, DFT_mat, WS_mat, dirs, TSR);CHKERRQ(ierr);
		
	}

	if (RSVDt->RSVD.q > 0)
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n******** Power iteration DONE! *************\n");CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}


