
#include <RecoverU.h>

PetscErrorCode RecoverU(Mat U_til, Mat Q_hat, Mat *U_hat)
{
	
	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	ierr = MatMatMult(Q_hat,U_til,MAT_INITIAL_MATRIX,PETSC_DEFAULT,U_hat);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}
