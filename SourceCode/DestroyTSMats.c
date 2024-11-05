
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode DestroyTSMats(TS_matrices *TS_mat)
{
	/*
		Removes the RK4 time-stepping matrices from memory
	*/  

	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	ierr = VecDestroy(&TS_mat->F1);CHKERRQ(ierr);
	ierr = VecDestroy(&TS_mat->F2);CHKERRQ(ierr);
	ierr = VecDestroy(&TS_mat->F3);CHKERRQ(ierr);
	ierr = VecDestroy(&TS_mat->k1);CHKERRQ(ierr);
	ierr = VecDestroy(&TS_mat->k2);CHKERRQ(ierr);
	ierr = VecDestroy(&TS_mat->k3);CHKERRQ(ierr);
	ierr = VecDestroy(&TS_mat->k4);CHKERRQ(ierr);
	ierr = VecDestroy(&TS_mat->y_temp);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}




