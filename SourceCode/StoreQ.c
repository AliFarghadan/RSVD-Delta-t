
#include <petscmat.h>
#include "Variables.h"

PetscErrorCode StoreQ(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt)
{
	/*
		Stores \hat{Y} on \hat{Q} before the final adjoint action
		to recover \hat{U} after the economy SVD
	*/

	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	ierr = MatDuplicate(RSVD_mat->Y_hat,MAT_COPY_VALUES,&RSVD_mat->Q_hat);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}




