
#include <petscmat.h>
#include "Variables.h"

PetscErrorCode LNSStructure(LNS_vars *LNS_mat, RSVDt_vars *RSVDt, Directories *dirs)
{
	/*
		Reads the LNS operator and applies discounting if desired
	*/

	PetscErrorCode        ierr;
	PetscLogDouble        t1, t2;
	PetscViewer           fd;

	PetscFunctionBeginUser;

	ierr = PetscTime(&t1);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->Mat_filename);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&LNS_mat->A);CHKERRQ(ierr);
	ierr = MatSetType(LNS_mat->A,MATMPIAIJ);CHKERRQ(ierr);
	ierr = MatLoad(LNS_mat->A,fd);CHKERRQ(ierr);
	ierr = MatGetSize(LNS_mat->A,&RSVDt->RSVD.N,NULL);CHKERRQ(ierr);

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"The matrix is loaded in %f seconds\n", t2-t1);CHKERRQ(ierr);

	/*
		Discounting
	*/

	if (LNS_mat->RSVDt.disc.flg_disc) {
		ierr = MatShift(LNS_mat->A,-LNS_mat->RSVDt.disc.Beta);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n---- Discounting with Beta = %g ----\n\n", LNS_mat->RSVDt.disc.Beta);CHKERRQ(ierr);
	}

	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}

