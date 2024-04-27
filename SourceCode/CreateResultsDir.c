
#include <unistd.h>
#include <petscmat.h>
#include "Variables.h"

PetscErrorCode CreateResultsDir(Directories *dirs, const char *FileName)
{
	/*
		Makes a folder in the root directory for the transient results
	*/

	PetscErrorCode        ierr;
	PetscInt              rank, FolderInd=0;
	char                  FolderName[PETSC_MAX_PATH_LEN];

	PetscFunctionBeginUser;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,(int*)&rank);CHKERRMPI(ierr);
	ierr = PetscSNPrintf((char*)&FolderName,PETSC_MAX_PATH_LEN,"%s%d/",FileName,(int)FolderInd);CHKERRQ(ierr);
	ierr = PetscSNPrintf((char*)&dirs->FolderDir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,FolderName);CHKERRQ(ierr);
	while (access(dirs->FolderDir, F_OK) == 0) {
		FolderInd++;
		ierr = PetscSNPrintf((char*)&FolderName,PETSC_MAX_PATH_LEN,"%s%d/",FileName,(int)FolderInd);CHKERRQ(ierr);
		ierr = PetscSNPrintf((char*)&dirs->FolderDir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,FolderName);CHKERRQ(ierr);
	}
	ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"mkdir %s",dirs->FolderDir);CHKERRQ(ierr);
	if ((int)rank == 0) ierr = system(dirs->IO_dir);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
	
}
