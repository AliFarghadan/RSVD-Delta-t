
#include <TransientRun.h>
#include <TransientRK4.h>

PetscErrorCode TransientRunRK4(RSVDt_vars *RSVDt, LNS_params *LNS_mat, WS_matrices *WS_mat, TS_matrices *TS_mat, Directories *dirs)
{
	
	PetscErrorCode        ierr;
	Mat                   Y_all, trans_norm;
	PetscBool             TransientRun_flg, TransientRun_dir;
	PetscInt              k_trans, rseed, mod = 10;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetBool(NULL,NULL,"-TransientSave_flg",&TransientSave_flg,&TransientSave_flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-TransientRun_flg",&TransientRun_flg,&TransientRun_flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-TransientRun_dir",&TransientRun_dir,&TransientRun_dir);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-k_trans",&k_trans,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-rseed",&rseed,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-mod",&mod,NULL);CHKERRQ(ierr);

	if (TransientRun_flg) {

		RSVDt->TS.flg_dir_adj = TransientRun_dir;

		ierr = MatCreate(PETSC_COMM_WORLD,&Y_all);CHKERRQ(ierr);
		ierr = MatCreate(PETSC_COMM_WORLD,&trans_norm);CHKERRQ(ierr);

		ierr = TransientRK4(Y_all,trans_norm,k_trans,mod,rseed,LNS_mat,WS_mat,TS_mat,RSVDt,TransientSave_flg,dirs);CHKERRQ(ierr);

		ierr = SlepcFinalize();
		return ierr;

	}

	PetscFunctionReturn(0);
	
}
