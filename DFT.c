
#include <DFT.h>

PetscErrorCode DFT(Mat Y_all, Mat dft, PetscInt time_ratio, Mat Y_hat, PetscBool flg_real_A)
{

	PetscErrorCode        ierr;
	Mat                   Y_temp1, Y_temp2, Y_hat_re;
	PetscInt              ik, Nw, N, k, Nw_eff;

	PetscFunctionBeginUser;

	ierr = MatGetSize(dft,&Nw,NULL);CHKERRQ(ierr);
	Nw_eff  = (flg_real_A) ? Nw/2 : Nw;
	ierr = MatGetSize(Y_all,&N,&k);CHKERRQ(ierr);
	k   /= Nw;
	ierr = MatCreate(PETSC_COMM_WORLD,&Y_hat_re);CHKERRQ(ierr);
	ierr = MatSetType(Y_hat_re,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Y_hat_re,PETSC_DECIDE,PETSC_DECIDE,N,k*Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(Y_hat_re);CHKERRQ(ierr);

	for (ik=0; ik<k; ik++) {
		ierr = MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nw,(ik+1)*Nw,&Y_temp1);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(Y_hat_re,PETSC_DECIDE,PETSC_DECIDE,ik*Nw_eff,(ik+1)*Nw_eff,&Y_temp2);CHKERRQ(ierr);
		ierr = MatMatMult(Y_temp1,dft,MAT_REUSE_MATRIX,PETSC_DEFAULT,&Y_temp2);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Y_hat_re,&Y_temp2);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Y_all,&Y_temp1);CHKERRQ(ierr);
	}
	ierr = MatAssemblyBegin(Y_hat_re,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Y_hat_re,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatScale(Y_hat_re, time_ratio);CHKERRQ(ierr);
	ierr = ReversePermuteForcing(Y_hat_re, Nw_eff, Y_hat);CHKERRQ(ierr);

	ierr = MatDestroy(&Y_all);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

