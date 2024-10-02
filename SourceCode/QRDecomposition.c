
#include <slepcbv.h>

PetscErrorCode QRDecomposition(Mat V)
{
	/*
		Performs QR decomposition on a given matrix 
	*/  

	PetscErrorCode        ierr;
	Mat                   Q_temp;
	BV                    Q;

	PetscFunctionBeginUser;

	ierr = BVCreateFromMat(V,&Q);CHKERRQ(ierr);
	ierr = BVSetType(Q,BVVECS);CHKERRQ(ierr);
	ierr = BVOrthogonalize(Q,NULL);CHKERRQ(ierr);
	ierr = BVCreateMat(Q,&Q_temp);CHKERRQ(ierr);
	ierr = MatCopy(Q_temp,V,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = BVDestroy(&Q);CHKERRQ(ierr);
	ierr = MatDestroy(&Q_temp);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}



