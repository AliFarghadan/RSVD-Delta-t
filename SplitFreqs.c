
#include <SplitFreqs.h>

PetscErrorCode SplitFreqs(Mat Y_hat_tall, PetscInt Nw, Mat *Y_hat)
{

	PetscErrorCode      ierr;
	PetscInt            N,k,loc_rows;
	PetscScalar        *array_mat;
	Mat                 Y_temp;
	
	PetscFunctionBeginUser;

	ierr = MatGetSize(Y_hat_tall,&N,&k);CHKERRQ(ierr);
	N   /= Nw;
	ierr = MatGetLocalSize(Y_hat_tall, &loc_rows, NULL);CHKERRQ(ierr);
	loc_rows /= Nw;
	ierr = MatDenseGetArray(Y_hat_tall, &array_mat);CHKERRQ(ierr);
	ierr = MatCreateDense(PETSC_COMM_WORLD, loc_rows, PETSC_DECIDE, N, k*Nw, array_mat, &Y_temp);CHKERRQ(ierr);
	ierr = MatDuplicate(Y_temp,MAT_COPY_VALUES,Y_hat);CHKERRQ(ierr);
	ierr = MatDenseRestoreArray(Y_hat_tall, &array_mat);CHKERRQ(ierr);

	ierr = MatDestroy(&Y_hat_tall);CHKERRQ(ierr);
	ierr = MatDestroy(&Y_temp);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}
