
#include <CreateForcingOnFly.h>

PetscErrorCode CreateForcingOnFly(Mat F_hat_re, Mat inv_dft, PetscInt jt_cyc, Vec F)
{

	PetscErrorCode        ierr;
	Vec                   W_read;

	PetscFunctionBeginUser;

	ierr = MatDenseGetColumnVecRead(inv_dft,jt_cyc,&W_read);CHKERRQ(ierr);
	ierr = MatMult(F_hat_re,W_read,F);CHKERRQ(ierr);
	ierr = MatDenseRestoreColumnVecRead(inv_dft,jt_cyc,&W_read);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

