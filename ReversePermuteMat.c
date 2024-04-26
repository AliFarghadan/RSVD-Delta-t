
#include <ReversePermuteMat.h>

PetscErrorCode ReversePermuteMat(Mat F_hat, PetscInt Nw, Mat F_hat_re)
{

	PetscErrorCode        ierr;
	PetscInt              iw, ik, ii, jj, N, k;
	Vec                   V1, V2;

	PetscFunctionBeginUser;

	ierr = MatGetSize(F_hat,&N,&k);CHKERRQ(ierr);
	k   /= Nw;
	ierr = MatSetSizes(F_hat_re,PETSC_DECIDE,PETSC_DECIDE,N,Nw*k);CHKERRQ(ierr);
	ierr = MatSetUp(F_hat_re);CHKERRQ(ierr);

	for (ik=0; ik<k; ik++) {
		for (iw=0; iw<Nw; iw++) {
			ii   = ik    + iw*k;
			jj   = ik*Nw + iw;
			ierr = MatDenseGetColumnVecRead(F_hat,jj,&V1);CHKERRQ(ierr);
			ierr = MatDenseGetColumnVecWrite(F_hat_re,ii,&V2);CHKERRQ(ierr);
			ierr = VecCopy(V1,V2);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(F_hat_re,ii,&V2);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecRead(F_hat,jj,&V1);CHKERRQ(ierr);
		}
	}
	
	ierr = MatAssemblyBegin(F_hat_re,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(F_hat_re,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	ierr = MatDestroy(&F_hat);CHKERRQ(ierr);
	ierr = VecDestroy(&V1);CHKERRQ(ierr);
	ierr = VecDestroy(&V2);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}
