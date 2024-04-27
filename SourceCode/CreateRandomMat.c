
#include <petscmat.h>
#include "Variables.h"
#include "ReversePermuteMat.h"

PetscErrorCode CreateRandomMat(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, Directories *dirs)
{
	/*
		Generates a random matrix of size N \times kNw (default)
		The input matrix can be read in using -Fin
		The latter option is useful for testing or resuming from a previous power iteration
	*/  

	PetscErrorCode        ierr;
	PetscBool             flg_Fin;
	Mat                   M;
	PetscViewer           fd;

	PetscFunctionBeginUser;

	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->Y_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVD_mat->Y_hat,MATDENSE);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL,NULL,"-Fin",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,&flg_Fin);CHKERRQ(ierr);
	if (flg_Fin) {
		// ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);
		// ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		// ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the input forcing matrix!\n");CHKERRQ(ierr);
		// ierr = MatLoad(RSVD_mat->Y_hat,fd);CHKERRQ(ierr);
		// ierr = MatGetSize(RSVD_mat->Y_hat,NULL,&RSVDt->RSVD.k);CHKERRQ(ierr);
		// RSVDt->RSVD.k /= RSVDt->RSVD.Nw_eff;

		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the input forcing matrix from previously obtained RSVD-delta-t!\n");CHKERRQ(ierr);
		ierr = MatSetSizes(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,RSVDt->RSVD.Nw_eff*RSVDt->RSVD.k);CHKERRQ(ierr);
		ierr = MatSetUp(RSVD_mat->Y_hat);CHKERRQ(ierr);
		for (PetscInt ik=0; ik<RSVDt->RSVD.k; ik++) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"k = %d ", (int) ik+1);CHKERRQ(ierr);
			ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%s%d%s",dirs->root_dir,dirs->output,"twin_resolvent_modes/V_hat_twin_k",(int) ik+1,"_allFreqs");CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
			ierr = MatDenseGetSubMatrix(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,ik*RSVDt->RSVD.Nw_eff,(ik+1)*RSVDt->RSVD.Nw_eff,&M);CHKERRQ(ierr);
			ierr = MatLoad(M,fd);CHKERRQ(ierr);
			ierr = MatDenseRestoreSubMatrix(RSVD_mat->Y_hat,&M);CHKERRQ(ierr);
		}
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
		ierr = ReversePermuteMat(RSVD_mat->Y_hat, RSVDt);CHKERRQ(ierr);
		
	} else {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Forcing input is generated randomly (uniform distribution)!\n");CHKERRQ(ierr);
		ierr = MatSetSizes(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,RSVDt->RSVD.Nw_eff*RSVDt->RSVD.k);CHKERRQ(ierr);
		ierr = MatSetUp(RSVD_mat->Y_hat);CHKERRQ(ierr);
		ierr = MatSetRandom(RSVD_mat->Y_hat,NULL);CHKERRQ(ierr);
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"@ N = %d, @ k = %d, @ q = %d, @ Nw = %d\n",\
				(int)RSVDt->RSVD.N,(int)RSVDt->RSVD.k,(int)RSVDt->RSVD.q,(int)RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

