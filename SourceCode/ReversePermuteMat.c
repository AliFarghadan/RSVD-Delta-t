
#include <petscmat.h>
#include "Variables.h"

PetscErrorCode ReversePermuteMat(Mat F, RSVDt_vars *RSVDt)
{
	/*
		Permutes a matrix of size N \times Nwk to N \times kNw (in place)
		k blocks of N \times Nw   ---->   Nw blocks of N \times k 
	*/

	PetscErrorCode        ierr;
	PetscInt              iw, ik, ii, jj, N, k, Nw;
	Mat                   F_permuted;
	Vec                   V1, V2;

	PetscFunctionBeginUser;

	ierr = MatGetSize(F,&N,&Nw);CHKERRQ(ierr);
	k    = RSVDt->RSVD.k;
	Nw  /= k;
	ierr = MatDuplicate(F,MAT_DO_NOT_COPY_VALUES,&F_permuted);CHKERRQ(ierr);

	for (ik=0; ik<k; ik++) {
		for (iw=0; iw<Nw; iw++) {
			ii   = ik    + iw*k;
			jj   = ik*Nw + iw;
			ierr = MatDenseGetColumnVecRead(F,jj,&V1);CHKERRQ(ierr);
			ierr = MatDenseGetColumnVecWrite(F_permuted,ii,&V2);CHKERRQ(ierr);
			ierr = VecCopy(V1,V2);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(F_permuted,ii,&V2);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecRead(F,jj,&V1);CHKERRQ(ierr);
		}
	}
	
	ierr = MatAssemblyBegin(F_permuted,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(F_permuted,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatCopy(F_permuted,F,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	
	ierr = MatDestroy(&F_permuted);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}


