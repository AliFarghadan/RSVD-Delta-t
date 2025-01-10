
#include <petscmat.h>
#include <Variables.h>
#include <ApplyWeightMats.h>
#include <TSActionRK4.h>

PetscErrorCode AdjointActionRK4(RSVD_matrices *RSVD, RSVDt_vars *RSVDt, LNS_vars *LNS, \
			DFT_matrices *DFT, Weight_matrices *Weight, Directories *dirs, TS_removal_matrices *TSR)
{
	/*
		Performs time-stepping to approximate R' \times \hat{F}, where (.)' indicates complex conjugate transpose
		For a modified resolvent operator, it computes B' * W_f_sqrt_inv' * R' * W_q_sqrt' * C' \times \hat{F} 
		In the latter case, the weight and input/output matrices are given as inputs
	*/
	
	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	/*
		1: direct action, 0: adjoint action
	*/
	
	RSVDt->TS.DirAdj = 0;

	ierr = ApplyWeightMats(RSVD, RSVDt, Weight, 1);CHKERRQ(ierr);
	ierr = TSActionRK4(RSVD, DFT, LNS, RSVDt, dirs, TSR);CHKERRQ(ierr);
	ierr = ApplyWeightMats(RSVD, RSVDt, Weight, 0);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}



