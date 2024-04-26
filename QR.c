
#include <QR.h>

PetscErrorCode QR(RSVD_matrices *RSVD_mat, PetscInt Nw)
{
	
	PetscErrorCode        ierr;
	Mat                   Y_temp, Y_hat_tall, Y_hat_re;
	BV                    Q;

	PetscFunctionBeginUser;

	ierr = MatCreate(PETSC_COMM_WORLD,&Y_hat_re);CHKERRQ(ierr);
	ierr = MatSetType(Y_hat_re,MATDENSE);CHKERRQ(ierr);
	ierr = PermuteForcing(RSVD_mat->Y_hat, Nw, Y_hat_re);CHKERRQ(ierr);

	ierr = MergeFrequencies(Y_hat_re,Nw,&Y_hat_tall);CHKERRQ(ierr);

	ierr = BVCreateFromMat(Y_hat_tall,&Q);CHKERRQ(ierr);
	ierr = BVSetType(Q,BVVECS);CHKERRQ(ierr);
	ierr = BVOrthogonalize(Q,NULL);CHKERRQ(ierr);
	ierr = BVCreateMat(Q,&Y_temp);CHKERRQ(ierr);
	ierr = BVDestroy(&Q);CHKERRQ(ierr);
	ierr = MatCopy(Y_temp,Y_hat_tall,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = MatDestroy(&Y_temp);CHKERRQ(ierr);

	ierr = SplitFrequencies(Y_hat_tall,Nw,&Y_hat_re);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->Y_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVD_mat->Y_hat,MATDENSE);CHKERRQ(ierr);
	ierr = ReversePermuteForcing(Y_hat_re, Nw, RSVD_mat->Y_hat);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

