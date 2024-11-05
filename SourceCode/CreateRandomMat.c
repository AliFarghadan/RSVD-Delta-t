
#include <petscmat.h>
#include <Variables.h>
#include <ReversePermuteMat.h>

PetscErrorCode CreateRandomMat(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, Directories *dirs)
{
	/*
		Generates a random matrix of size N x (k x Nw_eff) (default)
		The input matrix can be read in if desired
		The latter option is useful for testing or resuming from a previous power iteration
	*/  

	PetscErrorCode        ierr;
	PetscRandom           r;
	PetscViewer           fd;
	PetscInt              row;

	PetscFunctionBeginUser;

	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->Y_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVD_mat->Y_hat,MATDENSE);CHKERRQ(ierr);
	if (RSVDt->RSVD.InputForcingFlg) {
		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->RootDir,dirs->InputForcingDir);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the input forcing matrix         : %s\n", dirs->IO_dir);CHKERRQ(ierr);
		ierr = MatLoad(RSVD_mat->Y_hat,fd);CHKERRQ(ierr);
		ierr = MatGetSize(RSVD_mat->Y_hat,&row,&RSVDt->RSVD.k);CHKERRQ(ierr);
		if (row != RSVDt->RSVD.Nb) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Input forcing must be have %d rows, current size = %d x %d", (int)RSVDt->RSVD.Nb, (int)row, (int)RSVDt->RSVD.k);CHKERRQ(ierr);   
		RSVDt->RSVD.k /= RSVDt->RSVD.Nw_eff;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Based on the input forcing: k = %d\n", (int)RSVDt->RSVD.k);CHKERRQ(ierr);	
	} else {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Generating a random forcing matrix\n");CHKERRQ(ierr);
		ierr = MatSetSizes(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.Nb,RSVDt->RSVD.Nw_eff*RSVDt->RSVD.k);CHKERRQ(ierr);
		ierr = MatSetUp(RSVD_mat->Y_hat);CHKERRQ(ierr);
		ierr = PetscRandomCreate(PETSC_COMM_WORLD,&r);CHKERRQ(ierr);
		ierr = PetscRandomSetSeed(r, RSVDt->RSVD.RandSeed);CHKERRQ(ierr);
		ierr = PetscRandomSeed(r);CHKERRQ(ierr);
		ierr = MatSetRandom(RSVD_mat->Y_hat,NULL);CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
	
}

