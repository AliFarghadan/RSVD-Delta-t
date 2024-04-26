
#include <MergeFreqs.h>

PetscErrorCode MergeFreqs(Mat Y_hat, PetscInt Nw, Mat *Y_hat_tall)
{

	PetscErrorCode      ierr;
	PetscInt            N,k,loc_rows;
	PetscScalar        *array_mat;
	Mat                 Y_temp;

	PetscFunctionBeginUser;

	ierr = MatGetSize(Y_hat,&N,&k);CHKERRQ(ierr);
	k   /= Nw;
	ierr = MatGetLocalSize(Y_hat, &loc_rows, NULL);CHKERRQ(ierr);
	loc_rows *= Nw;
	ierr = MatDenseGetArray(Y_hat, &array_mat);CHKERRQ(ierr);
	ierr = MatCreateDense(PETSC_COMM_WORLD, loc_rows, PETSC_DECIDE, N*Nw, k, array_mat, &Y_temp);CHKERRQ(ierr);
	ierr = MatDuplicate(Y_temp,MAT_COPY_VALUES,Y_hat_tall);CHKERRQ(ierr);
	ierr = MatDenseRestoreArray(Y_hat, &array_mat);CHKERRQ(ierr);

	ierr = MatDestroy(&Y_hat);CHKERRQ(ierr);
	ierr = MatDestroy(&Y_temp);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
	
}

