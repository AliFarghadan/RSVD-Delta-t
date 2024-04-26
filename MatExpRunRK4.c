
#include <MatExpRunRK4.h>
#include <MatExpActionRK4.h>

PetscErrorCode MatExpRunRK4(RSVDt_vars *RSVDt, LNS_params *LNS_mat, WS_matrices *WS_mat, TS_matrices *TS_mat, Directories *dirs)
{
	
	PetscErrorCode        ierr;
	Mat                   Y_in, Y_out;
	PetscBool             period_flg;
	PetscInt              Ngap, st_ind = 1;
	PetscViewer           fd;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetString(NULL,NULL,"-input_subspace",(char*)&dirs.filename,PETSC_MAX_PATH_LEN,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-MatExpRun_flg",&MatExpRun_flg,&MatExpRun_flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-MatExpRun_dir",&MatExpRun_dir,&MatExpRun_dir);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-Ngap_period_flg",&period_flg,&period_flg);CHKERRQ(ierr);

	if (MatExpRun_flg) {

		RSVDt.TS.flg_dir_adj = MatExpAction_dir;
		
		ierr = MatCreate(PETSC_COMM_WORLD,&Y_out);CHKERRQ(ierr);
		Ngap = (period_flg) ? RSVDt->TS.Nt : RSVDt->TS.time_ratio_out;
		ierr = MatCreate(PETSC_COMM_WORLD,&Y_in);CHKERRQ(ierr);
		ierr = MatSetType(Y_in,MATDENSE);CHKERRQ(ierr);
		ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = MatLoad(Y_in,fd);CHKERRQ(ierr);
		ierr = MatExpActionRK4(Y_in,Y_out,Ngap,st_ind,LNS_mat,WS_mat,TS_mat,RSVDt,dirs);CHKERRQ(ierr);
		
		ierr = SlepcFinalize();
		return ierr;

	}

	PetscFunctionReturn(0);
	
}
