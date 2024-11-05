
#include <petscmat.h>
#include <Variables.h>
#include <PermuteMat.h>
#include <SaveInputVarsCopy.h>

PetscErrorCode SaveResults(Resolvent_matrices *Res_mat, RSVDt_vars *RSVDt, Directories *dirs)
{
	/*
		Saves the input/output modes and gains in the output directory
	*/
	
	PetscErrorCode        ierr;
	Mat                   M;
	PetscInt              ik, iw;
	PetscViewer           fd;

	PetscFunctionBeginUser;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n********************************************\n"
			"************ Saving the output *************\n********************************************\n");CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"\nThe results directory: %s\n\n",dirs->FolderDir);CHKERRQ(ierr);

	/*
		Saves resolvent modes for each mode separately (accross all frequencies) 
	*/

	ierr = PetscPrintf(PETSC_COMM_WORLD,"One matrix of size %d (modes) x %d (frequencies) for gains\n", \
			(int) RSVDt->RSVD.k, (int) RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->FolderDir,"S_hat");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = MatView(Res_mat->S_hat,fd);CHKERRQ(ierr);

	/*
		Nw matrix of size N x k for response modes and the same for forcing modes (SaveResults == 1)
		OR 
		k matrix of size N x Nw for response modes and the same for forcing modes (SaveResults == 2)
	*/

	if (RSVDt->SaveResultsOpt == 1) {

	ierr = PetscPrintf(PETSC_COMM_WORLD,"%d matrices of size %d x %d for response modes\n%d matrices of size %d x %d for forcing modes\n\n", \
			(int) RSVDt->RSVD.k, (int) RSVDt->RSVD.Nc, (int) RSVDt->RSVD.Nw_eff, \
			(int) RSVDt->RSVD.k, (int) RSVDt->RSVD.Nb, (int) RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);

		ierr = PermuteMat(Res_mat->U_hat, RSVDt);CHKERRQ(ierr);
		ierr = PermuteMat(Res_mat->V_hat, RSVDt);CHKERRQ(ierr);

		for (ik=0; ik<RSVDt->RSVD.k; ik++) {
			
			ierr = MatDenseGetSubMatrix(Res_mat->U_hat,PETSC_DECIDE,PETSC_DECIDE,ik*RSVDt->RSVD.Nw_eff,(ik+1)*RSVDt->RSVD.Nw_eff,&M);CHKERRQ(ierr);
			ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d%s",dirs->FolderDir,"U_hat_k",(int) ik+1,"_allFreqs");CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
			ierr = MatView(M,fd);CHKERRQ(ierr);
			ierr = MatDenseRestoreSubMatrix(Res_mat->U_hat,&M);CHKERRQ(ierr);

			ierr = MatDenseGetSubMatrix(Res_mat->V_hat,PETSC_DECIDE,PETSC_DECIDE,ik*RSVDt->RSVD.Nw_eff,(ik+1)*RSVDt->RSVD.Nw_eff,&M);CHKERRQ(ierr);
			ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d%s",dirs->FolderDir,"V_hat_k",(int) ik+1,"_allFreqs");CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
			ierr = MatView(M,fd);CHKERRQ(ierr);
			ierr = MatDenseRestoreSubMatrix(Res_mat->V_hat,&M);CHKERRQ(ierr);
		}
	} else { // RSVDt->SaveResultsOpt == 2

		ierr = PetscPrintf(PETSC_COMM_WORLD,"%d matrices of size %d x %d for response modes\n%d matrices of size %d x %d for forcing modes\n\n", \
			(int) RSVDt->RSVD.Nw_eff, (int) RSVDt->RSVD.Nc, (int) RSVDt->RSVD.k, \
			(int) RSVDt->RSVD.Nw_eff, (int) RSVDt->RSVD.Nb, (int) RSVDt->RSVD.k);CHKERRQ(ierr);

		for (iw=0; iw<RSVDt->RSVD.Nw_eff; iw++) {
			
			ierr = MatDenseGetSubMatrix(Res_mat->U_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&M);CHKERRQ(ierr);
			ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d%s",dirs->FolderDir,"U_hat_Freq",(int) iw,"_allK");CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
			ierr = MatView(M,fd);CHKERRQ(ierr);
			ierr = MatDenseRestoreSubMatrix(Res_mat->U_hat,&M);CHKERRQ(ierr);

			ierr = MatDenseGetSubMatrix(Res_mat->V_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&M);CHKERRQ(ierr);
			ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d%s",dirs->FolderDir,"V_hat_Freq",(int) iw,"_allK");CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
			ierr = MatView(M,fd);CHKERRQ(ierr);
			ierr = MatDenseRestoreSubMatrix(Res_mat->V_hat,&M);CHKERRQ(ierr);
		}
	}


	/*
		Saves a copy of the input variables in the results folder and exit
	*/

	ierr = SaveInputVarsCopy(dirs);CHKERRQ(ierr); 

	ierr = PetscPrintf(PETSC_COMM_WORLD,"DONE :))\n\n");CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}



