
#include <petscmat.h>
#include <petscviewerhdf5.h>
#include "Variables.h"
#include "CreateResultsDir.h"
#include "PermuteMat.h"

PetscErrorCode SaveResults(Resolvent_matrices *Res_mat, RSVDt_vars *RSVDt, Directories *dirs)
{
	/*
		Saves the input/output modes and gains in the output directory
	*/
	
	PetscErrorCode        ierr;
	Mat                   M;
	PetscViewer           fd;

	PetscFunctionBeginUser;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n********************************************\n"
			"************ Saving the output *************\n********************************************\n");CHKERRQ(ierr);

	/*
		Create a folder for the results
	*/

	ierr = CreateResultsDir(dirs, "ResolventModes_");CHKERRQ(ierr); 

	ierr = PetscPrintf(PETSC_COMM_WORLD,"\nThe output directory is: %s\n\n",dirs->FolderDir);CHKERRQ(ierr);

	/*
		Saving the results
	*/

	ierr = PetscPrintf(PETSC_COMM_WORLD,"One matrix of size %d (modes) x %d (frequencies) for gains\n", \
							(int) RSVDt->RSVD.k, (int) RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);
	
	if (Res_mat->all_in_one_flg) {

		/*
			Saves resolvent modes (all test vectors across all frequencies) all together
			Two matrices, each of size N x k x Nw, for response and forcing modes
		*/

		ierr = PetscPrintf(PETSC_COMM_WORLD,"Two matrices, each of size %d x %d x %d, for response and forcing modes\n", \
						(int) RSVDt->RSVD.N, (int) RSVDt->RSVD.k, (int) RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);

		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->FolderDir,"S_hat_clean_V1_RK4");CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
		ierr = MatView(Res_mat->S_hat,fd);CHKERRQ(ierr);

		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->FolderDir,"U_hat_clean_V1_RK4");CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
		ierr = MatView(Res_mat->U_hat,fd);CHKERRQ(ierr);

		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->FolderDir,"V_hat_clean_V1_RK4");CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
		ierr = MatView(Res_mat->V_hat,fd);CHKERRQ(ierr);

	} else {

		/*
			Saves resolvent modes for each mode separately (accross all frequencies) 
			Nw matrix of size N x k for response modes and the same for forcing modes
		*/

		ierr = PetscPrintf(PETSC_COMM_WORLD,"%d matrices of size %d x %d for response modes\n%d matrices of size %d x %d for forcing modes\n\n", \
				(int) RSVDt->RSVD.k, (int) RSVDt->RSVD.N, (int) RSVDt->RSVD.Nw_eff, (int) RSVDt->RSVD.k, (int) RSVDt->RSVD.N, (int) RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);

		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->FolderDir,"S_hat_clean_V1_RK4");CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
		ierr = MatView(Res_mat->S_hat,fd);CHKERRQ(ierr);

		ierr = PermuteMat(Res_mat->U_hat, RSVDt);CHKERRQ(ierr);
		ierr = PermuteMat(Res_mat->V_hat, RSVDt);CHKERRQ(ierr);

		for (PetscInt ik=0; ik<RSVDt->RSVD.k; ik++) {
			
			ierr = MatDenseGetSubMatrix(Res_mat->U_hat,PETSC_DECIDE,PETSC_DECIDE,ik*RSVDt->RSVD.Nw_eff,(ik+1)*RSVDt->RSVD.Nw_eff,&M);CHKERRQ(ierr);
			
			ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d%s",dirs->FolderDir,"U_hat_clean_V1_RK4_k",(int) ik+1,"_allFreqs.h5");CHKERRQ(ierr);
			ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, dirs->IO_dir, FILE_MODE_WRITE, &fd); CHKERRQ(ierr);
		    // ierr = PetscViewerSetFormat(fd, PETSC_VIEWER_HDF5);CHKERRQ(ierr);
		    // ierr = PetscViewerHDF5SetBaseDimension2(fd, PETSC_FALSE);CHKERRQ(ierr);
			ierr = PetscObjectSetName((PetscObject)M, "response");CHKERRQ(ierr);
			ierr = MatView(M,fd);CHKERRQ(ierr);

			// ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d%s",dirs->FolderDir,"U_hat_clean_V1_RK4_k",(int) ik+1,"_allFreqs");CHKERRQ(ierr);
			// ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
			// ierr = MatView(M,fd);CHKERRQ(ierr);
			ierr = MatDenseRestoreSubMatrix(Res_mat->U_hat,&M);CHKERRQ(ierr);

			ierr = MatDenseGetSubMatrix(Res_mat->V_hat,PETSC_DECIDE,PETSC_DECIDE,ik*RSVDt->RSVD.Nw_eff,(ik+1)*RSVDt->RSVD.Nw_eff,&M);CHKERRQ(ierr);
			ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d%s",dirs->FolderDir,"V_hat_clean_V1_RK4_k",(int) ik+1,"_allFreqs");CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
			ierr = MatView(M,fd);CHKERRQ(ierr);
			ierr = MatDenseRestoreSubMatrix(Res_mat->V_hat,&M);CHKERRQ(ierr);

		}
	}

	ierr = PetscPrintf(PETSC_COMM_WORLD,"DONE :))\n\n");CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}



