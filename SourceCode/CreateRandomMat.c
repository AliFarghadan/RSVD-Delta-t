
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode CreateRandomMat(RSVD_matrices *RSVD, RSVDt_vars *RSVDt, Directories *dirs)
{
	/*
		Generates a random matrix of size N x (k x Nw_eff) (default)
		The input matrix can be read in if desired
		The latter option is useful for testing or resuming from a previous power iteration
	*/  

	PetscErrorCode        ierr;
	PetscRandom           r;
	PetscViewer           fd;
	PetscInt              row, k;

	PetscFunctionBeginUser;

	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD->Y_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVD->Y_hat,MATDENSE);CHKERRQ(ierr);
	if (RSVDt->RSVD.InputForcingFlg) {
		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->RootDir,dirs->InputForcingDir);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the input forcing matrix         : %s\n", dirs->IO_dir);CHKERRQ(ierr);
		ierr = MatLoad(RSVD->Y_hat,fd);CHKERRQ(ierr);
		ierr = MatGetSize(RSVD->Y_hat,&row,&k);CHKERRQ(ierr);
		if (row != RSVDt->RSVD.Nb) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Input forcing must be have %d rows, current size = %d x %d", (int)RSVDt->RSVD.Nb, (int)row, (int)RSVDt->RSVD.k);CHKERRQ(ierr); 
		if (PetscFmodReal(k,RSVDt->RSVD.Nw_eff) != 0) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER, \
					"Input forcing number of columns must be factor of %d, current number of columns = %d", \
								(int)RSVDt->RSVD.Nw_eff, (int)k);CHKERRQ(ierr);   
		k /= RSVDt->RSVD.Nw_eff;
		if (k != RSVDt->RSVD.k) {
			RSVDt->RSVD.k = k;
			if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"Based on the input forcing: k = %d\n", (int)RSVDt->RSVD.k);CHKERRQ(ierr);
		}
	} else {
		if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"Generating a random forcing matrix\n");CHKERRQ(ierr);
		ierr = MatSetSizes(RSVD->Y_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.Nb,RSVDt->RSVD.Nw_eff*RSVDt->RSVD.k);CHKERRQ(ierr);
		ierr = MatSetUp(RSVD->Y_hat);CHKERRQ(ierr);
		ierr = PetscRandomCreate(PETSC_COMM_WORLD,&r);CHKERRQ(ierr);
		ierr = PetscRandomSetSeed(r, RSVDt->RSVD.RandSeed);CHKERRQ(ierr);
		ierr = PetscRandomSeed(r);CHKERRQ(ierr);
		ierr = MatSetRandom(RSVD->Y_hat,r);CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
	
}

