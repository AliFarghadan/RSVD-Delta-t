
#include <petscmat.h>
#include <Variables.h>
#include <ReversePermuteMat.h>

PetscErrorCode DiscreteFourierTransform(Mat Y_all, DFT_matrices *DFT, RSVD_matrices *RSVD, RSVDt_vars *RSVDt)
{
	/*
		Performs DFT on the collected response snapshots to obtain \hat{Y}
	*/  

	PetscErrorCode        ierr;
	Mat                   Y1, Y2;
	PetscInt              ik, Nw, N, k;

	PetscFunctionBeginUser;

	ierr = MatGetSize(DFT->dft,&Nw,NULL);CHKERRQ(ierr);
	ierr = MatGetSize(Y_all,&N,&k);CHKERRQ(ierr);
	k   /= Nw;
	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD->Y_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVD->Y_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(RSVD->Y_hat,PETSC_DECIDE,PETSC_DECIDE,N,k*RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(RSVD->Y_hat);CHKERRQ(ierr);

	for (ik=0; ik<k; ik++) {
		ierr = MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nw,(ik+1)*Nw,&Y1);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(RSVD->Y_hat,PETSC_DECIDE,PETSC_DECIDE,ik*RSVDt->RSVD.Nw_eff,(ik+1)*RSVDt->RSVD.Nw_eff,&Y2);CHKERRQ(ierr);
		ierr = MatMatMult(Y1,DFT->dft,MAT_REUSE_MATRIX,PETSC_DEFAULT,&Y2);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(RSVD->Y_hat,&Y2);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Y_all,&Y1);CHKERRQ(ierr);
	}

	ierr = MatDestroy(&Y_all);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(RSVD->Y_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(RSVD->Y_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatScale(RSVD->Y_hat, RSVDt->TS.ResRatio);CHKERRQ(ierr);
	ierr = ReversePermuteMat(RSVD->Y_hat, RSVDt);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}


