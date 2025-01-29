
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode ApplyWeightMats(RSVD_matrices *RSVD, RSVDt_vars *RSVDt, Weight_matrices *Weight, PetscBool before)
{
	/*
		Applies weight, input and output matrices if defined
	*/

	PetscErrorCode        ierr=0;
	Mat                   Y;

	PetscFunctionBeginUser;

	if (DirAdj) { // direct 
		if (before) { // forcing 
			if (Weight_mat->InvInputWeightFlg)  {
				ierr = MatMatMult(Weight_mat->W_f_sqrt_inv,RSVD_mat->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatCopy(Y,RSVD_mat->Y_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}
			if (Weight_mat->InputMatrixFlg)  {
				ierr = MatMatMult(Weight_mat->B,RSVD_mat->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatDuplicate(Y,MAT_COPY_VALUES,&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}
		} else { // response 
			if (Weight_mat->OutputMatrixFlg) {
				ierr = MatMatMult(Weight_mat->C,RSVD_mat->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatDestroy(&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatDuplicate(Y,MAT_COPY_VALUES,&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}
			if (Weight_mat->OutputWeightFlg) {
				ierr = MatMatMult(Weight_mat->W_q_sqrt,RSVD_mat->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatCopy(Y,RSVD_mat->Y_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}		
		}
	} else { // adjoint
		if (before) { // forcing
			if (Weight_mat->OutputWeightFlg) {
				ierr = MatHermitianTranspose(Weight_mat->W_q_sqrt, MAT_INPLACE_MATRIX, &Weight_mat->W_q_sqrt);CHKERRQ(ierr);
				ierr = MatMatMult(Weight_mat->W_q_sqrt,RSVD_mat->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatCopy(Y,RSVD_mat->Y_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
				ierr = MatHermitianTranspose(Weight_mat->W_q_sqrt, MAT_INPLACE_MATRIX, &Weight_mat->W_q_sqrt);CHKERRQ(ierr);
			}
			if (Weight_mat->OutputMatrixFlg) {
				ierr = MatHermitianTranspose(Weight_mat->C, MAT_INPLACE_MATRIX, &Weight_mat->C);CHKERRQ(ierr);
				ierr = MatMatMult(Weight_mat->C,RSVD_mat->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatHermitianTranspose(Weight_mat->C, MAT_INPLACE_MATRIX, &Weight_mat->C);CHKERRQ(ierr);
				ierr = MatDestroy(&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatDuplicate(Y,MAT_COPY_VALUES,&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}
		} else { // response
			if (Weight_mat->InputMatrixFlg)  {
				ierr = MatHermitianTranspose(Weight_mat->B, MAT_INPLACE_MATRIX, &Weight_mat->B);CHKERRQ(ierr);
				ierr = MatMatMult(Weight_mat->B,RSVD_mat->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatHermitianTranspose(Weight_mat->B, MAT_INPLACE_MATRIX, &Weight_mat->B);CHKERRQ(ierr);
				ierr = MatDestroy(&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatDuplicate(Y,MAT_COPY_VALUES,&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}
			if (Weight_mat->InvInputWeightFlg)  {
				ierr = MatHermitianTranspose(Weight_mat->W_f_sqrt_inv, MAT_INPLACE_MATRIX, &Weight_mat->W_f_sqrt_inv);CHKERRQ(ierr);
				ierr = MatMatMult(Weight_mat->W_f_sqrt_inv,RSVD_mat->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatCopy(Y,RSVD_mat->Y_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);	
				ierr = MatHermitianTranspose(Weight_mat->W_f_sqrt_inv, MAT_INPLACE_MATRIX, &Weight_mat->W_f_sqrt_inv);CHKERRQ(ierr);
			}
		}
	}

	PetscFunctionReturn(0);
	
}




