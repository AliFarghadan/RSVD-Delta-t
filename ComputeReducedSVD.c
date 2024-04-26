
#include <ComputeReducedSVD.h>
#include <SVD.h>

PetscErrorCode ComputeReducedSVD(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, resolvent_matrices *res_mat, Directories *dirs)
{
	
	PetscErrorCode        ierr;
	PetscViewer           fd;

	PetscFunctionBeginUser;

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"Y_hat_BeforeSVD_RK4");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = MatView(RSVD_mat->Y_hat,fd);CHKERRQ(ierr);	

	ierr = SVD(RSVD_mat,res_mat,RSVDt);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Saving onto the disk in progress...\n");

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"S_hat_clean_V1_RK4");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = VecView(res_mat->S_hat,fd);CHKERRQ(ierr);

	// if (WS_mat->flg_Wout) {
	// 	if (WS_mat->flg_diag) {
	// 		ierr = MatDiagonalScale(res_mat->U_hat,WS_mat->W_out_vec,NULL);CHKERRQ(ierr);
	// 		ierr = MatDiagonalScale(res_mat->V_hat,WS_mat->W_out_vec,NULL);CHKERRQ(ierr);
	// 	}
	// 	else {
	// 		ierr = MatMatMult(WS_mat->W_out,res_mat->U_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&res_mat->U_hat);CHKERRQ(ierr);
	// 		ierr = MatMatMult(WS_mat->W_out,res_mat->V_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&res_mat->V_hat);CHKERRQ(ierr);
	// 	}
	// }

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"V_hat_clean_V1_RK4");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = MatView(res_mat->V_hat,fd);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"U_hat_clean_V1_RK4");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = MatView(res_mat->U_hat,fd);CHKERRQ(ierr);	

	ierr = PetscPrintf(PETSC_COMM_WORLD,"DONE :))\n");CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

