
#include <petscmat.h>
#include "Variables.h"

PetscErrorCode SaveInputVarsCopy(Directories *dirs)
{
	/*
		Saves a copy of the input variables in the results folder
	*/

	PetscErrorCode        ierr;
	char                  filename[PETSC_MAX_PATH_LEN];
	char                  pwd[PETSC_MAX_PATH_LEN];
	char                  src[PETSC_MAX_PATH_LEN];
	char                  dst[PETSC_MAX_PATH_LEN];
	char                  buffer[1024];
	FILE                 *src_file;
	FILE                 *dst_file;
	size_t                bytes_read;
	PetscBool             flg_set;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetString(NULL, NULL,"-inputs",(char*)&filename,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	if (flg_set) {
		ierr = PetscGetWorkingDirectory(pwd, PETSC_MAX_PATH_LEN);CHKERRQ(ierr);
		ierr = PetscSNPrintf((char*)&src,PETSC_MAX_PATH_LEN,"%s/%s",pwd,filename);CHKERRQ(ierr);
		ierr = PetscSNPrintf((char*)&dst,PETSC_MAX_PATH_LEN,"%s%s",dirs->FolderDir,filename);CHKERRQ(ierr);

		ierr = PetscFOpen(MPI_COMM_SELF, src, "r", &src_file);CHKERRQ(ierr);
		ierr = PetscFOpen(MPI_COMM_SELF, dst, "w", &dst_file);CHKERRQ(ierr);
		while ((bytes_read = fread(buffer, 1, 1024, src_file)) > 0) {
		    fwrite(buffer, 1, bytes_read, dst_file);
		}
		ierr = PetscFClose(MPI_COMM_SELF, src_file);CHKERRQ(ierr);
		ierr = PetscFClose(MPI_COMM_SELF, dst_file);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
	
}
