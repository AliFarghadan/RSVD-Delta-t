
#include <ReadSingularModes.h>

PetscErrorCode ReadSingularModes(Sing_matrices *Sing_mat, RSVDt_vars *RSVDt, Directories *dirs)
{

	PetscErrorCode        ierr;
	Mat                   V, U;

	PetscFunctionBeginUser;

	if (Singular_flg) {

		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the singular modes!\n");

		ierr = MatCreate(PETSC_COMM_WORLD,&V);CHKERRQ(ierr);
		ierr = MatSetType(V,MATDENSE);CHKERRQ(ierr);
		ierr = PetscOptionsGetString(NULL,NULL,"-SingV",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,NULL);CHKERRQ(ierr);
		ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs.filename);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = MatLoad(V,fd);CHKERRQ(ierr);	
		ierr = MergeFreqs(V,RSVDt->RSVD.Nw,&Sing_mat->SingV);CHKERRQ(ierr);

		ierr = MatCreate(PETSC_COMM_WORLD,&U);CHKERRQ(ierr);
		ierr = MatSetType(U,MATDENSE);CHKERRQ(ierr);
		ierr = PetscOptionsGetString(NULL,NULL,"-SingU",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,NULL);CHKERRQ(ierr);
		ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs.filename);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = MatLoad(U,fd);CHKERRQ(ierr);	
		ierr = MergeFreqs(U,RSVDt->RSVD.Nw,&Sing_mat->SingU);CHKERRQ(ierr);

	}

	PetscFunctionReturn(0);
	
}

