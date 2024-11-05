
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode CreateForcingOnFly(Mat F_hat, DFT_matrices *DFT_mat, PetscInt jt_cyc, Vec F)
{

	/*
		Creates a temporary forcing snapshot in time from \hat{F} in Fourier space
	*/  

	PetscErrorCode        ierr;
	Vec                   W_col;

	PetscFunctionBeginUser;

	ierr = MatDenseGetColumnVecRead(DFT_mat->idft,jt_cyc,&W_col);CHKERRQ(ierr);
	ierr = MatMult(F_hat,W_col,F);CHKERRQ(ierr);
	ierr = MatDenseRestoreColumnVecRead(DFT_mat->idft,jt_cyc,&W_col);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

