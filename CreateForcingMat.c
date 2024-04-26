
#include <CreateForcingMat.h>

PetscErrorCode CreateForcingMat(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, Directories *dirs)
{

	PetscErrorCode        ierr;
	PetscBool             flg_Fin;
	PetscViewer           fd;

	PetscFunctionBeginUser;

	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->Y_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVD_mat->Y_hat,MATDENSE);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL,NULL,"-Fin",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,&flg_Fin);CHKERRQ(ierr);
	if (flg_Fin) {
		ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the input forcing matrix!\n");
		ierr = MatLoad(RSVD_mat->Y_hat,fd);CHKERRQ(ierr);
		ierr = MatGetSize(RSVD_mat->Y_hat,NULL,&RSVDt->RSVD.k);CHKERRQ(ierr);
		RSVDt->RSVD.k   /= (RSVDt->TS.flg_real_A) ? RSVDt->RSVD.Nw/2 : RSVDt->RSVD.Nw;
		if (RSVDt->TS.flg_real_A) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"@ N = %d, @ k = %d, @ Nw = %d\n",(int)RSVDt->RSVD.N,(int)RSVDt->RSVD.k,(int)RSVDt->RSVD.Nw/2);CHKERRQ(ierr);
		} else {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"@ N = %d, @ k = %d, @ Nw = %d\n",(int)RSVDt->RSVD.N,(int)RSVDt->RSVD.k,(int)RSVDt->RSVD.Nw);CHKERRQ(ierr);
		}
	} else {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Forcing input is generated randomly (uniform distribution)!\n");CHKERRQ(ierr);
		if (RSVDt->TS.flg_real_A) {
			ierr = MatSetSizes(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,RSVDt->RSVD.Nw/2*RSVDt->RSVD.k);CHKERRQ(ierr);
		} else {
			ierr = MatSetSizes(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,RSVDt->RSVD.Nw*RSVDt->RSVD.k);CHKERRQ(ierr);
		}
		ierr = MatSetUp(RSVD_mat->Y_hat);CHKERRQ(ierr);
		ierr = MatSetRandom(RSVD_mat->Y_hat,NULL);CHKERRQ(ierr);
		if (RSVDt->TS.flg_real_A) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"@ N = %d, @ k = %d, @ Nw = %d\n",(int)RSVDt->RSVD.N,(int)RSVDt->RSVD.k,(int)RSVDt->RSVD.Nw/2);CHKERRQ(ierr);
		} else {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"@ N = %d, @ k = %d, @ Nw = %d\n",(int)RSVDt->RSVD.N,(int)RSVDt->RSVD.k,(int)RSVDt->RSVD.Nw);CHKERRQ(ierr);
		}
	}		

	PetscFunctionReturn(0);
	
}

