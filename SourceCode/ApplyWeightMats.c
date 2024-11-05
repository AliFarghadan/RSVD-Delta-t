
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode ApplyWeightMats(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, Weight_matrices *Weight_mat, PetscBool before)
{
	/*
		Applies weight, input and output matrices if defined
	*/

	PetscErrorCode        ierr=0;
	Mat                   Y;

	PetscFunctionBeginUser;

	if (RSVDt->TS.DirAdj) { // direct 
		if (before) { // forcing 
			if (Weight_mat->InvInputWeightFlg)  ierr = MatMatMult(Weight_mat->W_f_sqrt_inv,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);
			if (Weight_mat->InputMatrixFlg)  {
				ierr = MatMatMult(Weight_mat->B,RSVD_mat->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatDestroy(&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatSetType(RSVD_mat->Y_hat,MATDENSE);CHKERRQ(ierr);
				ierr = MatSetSizes(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,RSVDt->RSVD.Nw_eff*RSVDt->RSVD.k);CHKERRQ(ierr);
				ierr = MatSetUp(RSVD_mat->Y_hat);CHKERRQ(ierr);				
				ierr = MatCopy(Y,RSVD_mat->Y_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}
		} else { // response 
			if (Weight_mat->OutputMatrixFlg) {
				ierr = MatMatMult(Weight_mat->C,RSVD_mat->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatDestroy(&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatSetType(RSVD_mat->Y_hat,MATDENSE);CHKERRQ(ierr);
				ierr = MatSetSizes(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.Nc,RSVDt->RSVD.Nw_eff*RSVDt->RSVD.k);CHKERRQ(ierr);
				ierr = MatSetUp(RSVD_mat->Y_hat);CHKERRQ(ierr);			
				ierr = MatCopy(Y,RSVD_mat->Y_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}
			if (Weight_mat->OutputWeightFlg) ierr = MatMatMult(Weight_mat->W_q_sqrt,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);
		}
	} else { // adjoint
		if (before) { // forcing
			if (Weight_mat->OutputWeightFlg) {
				ierr = MatHermitianTranspose(Weight_mat->W_q_sqrt, MAT_INPLACE_MATRIX, &Weight_mat->W_q_sqrt);CHKERRQ(ierr);
				ierr = MatMatMult(Weight_mat->W_q_sqrt,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatHermitianTranspose(Weight_mat->W_q_sqrt, MAT_INPLACE_MATRIX, &Weight_mat->W_q_sqrt);CHKERRQ(ierr);
			}
			if (Weight_mat->OutputMatrixFlg) {
				ierr = MatHermitianTranspose(Weight_mat->C, MAT_INPLACE_MATRIX, &Weight_mat->C);CHKERRQ(ierr);
				ierr = MatMatMult(Weight_mat->C,RSVD_mat->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatHermitianTranspose(Weight_mat->C, MAT_INPLACE_MATRIX, &Weight_mat->C);CHKERRQ(ierr);
				ierr = MatDestroy(&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatSetType(RSVD_mat->Y_hat,MATDENSE);CHKERRQ(ierr);
				ierr = MatSetSizes(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,RSVDt->RSVD.Nw_eff*RSVDt->RSVD.k);CHKERRQ(ierr);
				ierr = MatSetUp(RSVD_mat->Y_hat);CHKERRQ(ierr);				
				ierr = MatCopy(Y,RSVD_mat->Y_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}
		} else { // response
			if (Weight_mat->InputMatrixFlg)  {
				ierr = MatHermitianTranspose(Weight_mat->B, MAT_INPLACE_MATRIX, &Weight_mat->B);CHKERRQ(ierr);
				ierr = MatMatMult(Weight_mat->B,RSVD_mat->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatHermitianTranspose(Weight_mat->B, MAT_INPLACE_MATRIX, &Weight_mat->B);CHKERRQ(ierr);
				ierr = MatDestroy(&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatSetType(RSVD_mat->Y_hat,MATDENSE);CHKERRQ(ierr);
				ierr = MatSetSizes(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.Nb,RSVDt->RSVD.Nw_eff*RSVDt->RSVD.k);CHKERRQ(ierr);
				ierr = MatSetUp(RSVD_mat->Y_hat);CHKERRQ(ierr);				
				ierr = MatCopy(Y,RSVD_mat->Y_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}
			if (Weight_mat->InvInputWeightFlg)  {
				ierr = MatHermitianTranspose(Weight_mat->W_f_sqrt_inv, MAT_INPLACE_MATRIX, &Weight_mat->W_f_sqrt_inv);CHKERRQ(ierr);
				ierr = MatMatMult(Weight_mat->W_f_sqrt_inv,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);
				ierr = MatHermitianTranspose(Weight_mat->W_f_sqrt_inv, MAT_INPLACE_MATRIX, &Weight_mat->W_f_sqrt_inv);CHKERRQ(ierr);
			}
		}
	}

	PetscFunctionReturn(0);
	
}



