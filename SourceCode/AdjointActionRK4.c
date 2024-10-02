
#include <petscmat.h>
#include "Variables.h"
#include "ApplyWSMats.h"
#include "TSActionRK4.h"

PetscErrorCode AdjointActionRK4(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, LNS_vars *LNS_mat, \
			DFT_matrices *DFT_mat, WS_matrices *WS_mat, Directories *dirs, TS_removal_matrices *TSR)
{
	/*
		Performs time-stepping to approximate R' \times \hat{F}, where (.)' indicates complex conjugate transpose
		For a modified resolvent operator, it computes B' * W_out' * R' * W_in' * C' * \times \hat{F} 
		In the latter case, the weight and input/output matrices are given as inputs
	*/
	
	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	/*
		1: direct action, 0: adjoint action
	*/
	
	RSVDt->TS.DirAdj = 0;

	ierr = ApplyWSMats(RSVD_mat, RSVDt, WS_mat, 1);CHKERRQ(ierr);
	ierr = TSActionRK4(RSVD_mat, DFT_mat, LNS_mat, RSVDt, dirs, TSR);CHKERRQ(ierr);
	ierr = ApplyWSMats(RSVD_mat, RSVDt, WS_mat, 0);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}



