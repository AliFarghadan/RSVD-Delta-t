
#include <petscmat.h>
#include "Variables.h"
#include "ApplyWSMats.h"
#include "TSActionRK4.h"

PetscErrorCode ForwardActionRK4(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, LNS_vars *LNS_mat, \
			DFT_matrices *DFT_mat, WS_matrices *WS_mat, Directories *dirs, TS_removal_matrices *TSR)
{
	/*
		Performs time-stepping to approximate $R \times \hat{F}$ 
		For a modified resolvent operator, it computes C * W_in * R * W_out * B * \times \hat{F} 
		In the latter case, the weight and input/output matrices are given as inputs
	*/

	PetscErrorCode        ierr;

	PetscFunctionBeginUser;
		
	RSVDt->TS.flg_dir_adj = 1;

	ierr = ApplyWSMats(RSVD_mat, RSVDt, WS_mat, 1);CHKERRQ(ierr);
	ierr = TSActionRK4(RSVD_mat, DFT_mat, LNS_mat, RSVDt, dirs, TSR);CHKERRQ(ierr);
	ierr = ApplyWSMats(RSVD_mat, RSVDt, WS_mat, 0);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}



