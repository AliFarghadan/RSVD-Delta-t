
#include <petscmat.h>
#include "Variables.h"
#include "ReversePermuteMat.h"

PetscErrorCode DFT(Mat Y_all, DFT_matrices *DFT_mat, RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt)
{
	/*
		Performs DFT on the collected response snapshots to obtain \hat{Y}
	*/  

	PetscErrorCode        ierr;
	Mat                   Y_temp1, Y_temp2;
	PetscInt              ik, Nw, N, k;

	PetscFunctionBeginUser;

	ierr = MatGetSize(DFT_mat->dft,&Nw,NULL);CHKERRQ(ierr);
	ierr = MatGetSize(Y_all,&N,&k);CHKERRQ(ierr);
	k   /= Nw;
	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->Y_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVD_mat->Y_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,N,k*RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(RSVD_mat->Y_hat);CHKERRQ(ierr);

	for (ik=0; ik<k; ik++) {
		ierr = MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nw,(ik+1)*Nw,&Y_temp1);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,ik*RSVDt->RSVD.Nw_eff,(ik+1)*RSVDt->RSVD.Nw_eff,&Y_temp2);CHKERRQ(ierr);
		ierr = MatMatMult(Y_temp1,DFT_mat->dft,MAT_REUSE_MATRIX,PETSC_DEFAULT,&Y_temp2);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(RSVD_mat->Y_hat,&Y_temp2);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Y_all,&Y_temp1);CHKERRQ(ierr);
	}

	ierr = MatDestroy(&Y_all);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(RSVD_mat->Y_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(RSVD_mat->Y_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatScale(RSVD_mat->Y_hat, RSVDt->TS.ResRatio);CHKERRQ(ierr);
	ierr = ReversePermuteMat(RSVD_mat->Y_hat, RSVDt);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}


