
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode StoreU(RSVD_matrices *RSVD_mat, Resolvent_matrices *Res_mat)
{
	/*
		Stores \hat{Y} on \hat{U}
	*/

	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	ierr = MatDuplicate(RSVD_mat->Y_hat,MAT_COPY_VALUES,&Res_mat->U_hat);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}




