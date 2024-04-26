
#include <TSTimeDomainRunRK4.h>
#include <TSActionRK4.h>

PetscErrorCode TSTimeDomainRunRK4(RSVDt_vars *RSVDt, LNS_params *LNS_mat, WS_matrices *WS_mat, TS_matrices *TS_mat, Directories *dirs)
{
	
	PetscErrorCode        ierr;
	Mat                   Y_all, trans_norm;
	PetscBool             TSTimeDomainRun_flg, TSTimeDomainRun_dir;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetBool(NULL,NULL,"-TSTimeDomainRun_dir",&TSTimeDomainRun_dir,&TSTimeDomainRun_dir);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-TSTimeDomainRun_flg",&TSTimeDomainRun_flg,&TSTimeDomainRun_flg);CHKERRQ(ierr);

	if (TSTimeDomainRun_flg) {

		RSVDt->TS.flg_dir_adj = TSTimeDomainRun_dir;

		if (RSVDt->TS.flg_dir_adj) {	
			if (WS_mat->flg_Wout) {
				if (WS_mat->flg_diag) {
					ierr = MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->W_out_vec,NULL);CHKERRQ(ierr);
				}
				else {
					ierr = MatMatMult(WS_mat->W_out,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);
				}
			}
			if (WS_mat->flg_B)  ierr = MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->B_vec,NULL);CHKERRQ(ierr);
			if (Singular_flg) ierr = ProjectOut(RSVD_mat, RSVDt, Sing_mat->SingU);CHKERRQ(ierr);
		} else {
			if (WS_mat->flg_Win) {
				if (WS_mat->flg_diag) {
					ierr = MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->W_in_vec,NULL);CHKERRQ(ierr);
				}
				else {
					ierr = MatMatMult(WS_mat->W_in,RSVD_mat->Y_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&RSVD_mat->Y_hat);CHKERRQ(ierr);
				}
			}
			if (WS_mat->flg_C)  ierr = MatDiagonalScale(RSVD_mat->Y_hat,WS_mat->C_vec,NULL);CHKERRQ(ierr);
			if (Singular_flg) ierr = ProjectOut(RSVD_mat, RSVDt, &Sing_mat->SingV);CHKERRQ(ierr);
		}

		ierr = TSActionRK4(RSVD_mat,DFT_mat,LNS_mat,TS_mat,RSVDt,dirs);CHKERRQ(ierr);	

		ierr = SlepcFinalize();
		return ierr;
	}

	PetscFunctionReturn(0);
	
}
