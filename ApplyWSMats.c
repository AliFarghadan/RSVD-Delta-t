
#include <ApplyWSMats.h>

PetscErrorCode ApplyWSMats(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, WS_matrices *WS_mat, PetscBool before)
{

	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	if (RSVDt->TS.flg_dir_adj) {
		if (before) {
			if (WS_mat->flg_Wout) {
				ierr = WS_mat->flg_diag ? MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->W_out_vec,NULL) : \
					MatMatMult(WS_mat->W_out,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);
			}
			if (WS_mat->flg_B) ierr = MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->B_vec,NULL);CHKERRQ(ierr);
		} else {
			if (WS_mat->flg_C) ierr = MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->C_vec,NULL);CHKERRQ(ierr);
			if (WS_mat->flg_Win) {
				ierr = WS_mat->flg_diag ? MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->W_in_vec,NULL) : \
					MatMatMult(WS_mat->W_in,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);
		}
	} else {
		if (before) {
			if (WS_mat->flg_Win) {
				ierr = WS_mat->flg_diag ? MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->W_in_vec,NULL) : \
					MatMatMult(WS_mat->W_in,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);
			}
			if (WS_mat->flg_C) ierr = MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->C_vec,NULL);CHKERRQ(ierr);
		} else {
			if (WS_mat->flg_B) ierr = MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->B_vec,NULL);CHKERRQ(ierr);
			if (WS_mat->flg_Wout) {
				ierr = WS_mat->flg_diag ? MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->W_out_vec,NULL) : \
					MatMatMult(WS_mat->W_out,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);
			}
	}

	PetscFunctionReturn(0);
	
}

