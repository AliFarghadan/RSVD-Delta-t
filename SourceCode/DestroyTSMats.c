
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode DestroyTSMats(TS_matrices *TS)
{
	/*
		Removes the RK4 time-stepping matrices from memory
	*/  

	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	ierr = VecDestroy(&TS->F1);CHKERRQ(ierr);
	ierr = VecDestroy(&TS->F2);CHKERRQ(ierr);
	ierr = VecDestroy(&TS->F3);CHKERRQ(ierr);
	ierr = VecDestroy(&TS->k1);CHKERRQ(ierr);
	ierr = VecDestroy(&TS->k2);CHKERRQ(ierr);
	ierr = VecDestroy(&TS->k3);CHKERRQ(ierr);
	ierr = VecDestroy(&TS->k4);CHKERRQ(ierr);
	ierr = VecDestroy(&TS->y_temp);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}




