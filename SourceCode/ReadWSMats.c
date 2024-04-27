
#include <petscmat.h>
#include "Variables.h"

PetscErrorCode ReadWSMats(WS_matrices *WS_mat, Directories *dirs)
{
	/*
		Read the weight, input and outputs matrices if given in the jobscript
	*/

	PetscErrorCode        ierr;
	PetscInt              N;
	Mat                   B,C;
	PetscViewer           fd;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetString(NULL,NULL,"-Win",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,&WS_mat->flg_Win);CHKERRQ(ierr);
	if (WS_mat->flg_Win) {
		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);
		ierr = MatCreate(PETSC_COMM_WORLD,&WS_mat->W_in);CHKERRQ(ierr);
		ierr = MatSetType(WS_mat->W_in,MATMPIAIJ);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the input weight matrix!\n");
		ierr = MatLoad(WS_mat->W_in,fd);CHKERRQ(ierr);
		if (WS_mat->flg_diag) {
			ierr = MatGetSize(WS_mat->W_in,&N,NULL);CHKERRQ(ierr);
			ierr = VecCreate(PETSC_COMM_WORLD,&WS_mat->W_in_vec);CHKERRQ(ierr);
			ierr = VecSetSizes(WS_mat->W_in_vec,PETSC_DECIDE,N);CHKERRQ(ierr);
			ierr = VecSetUp(WS_mat->W_in_vec);CHKERRQ(ierr);
			ierr = MatGetDiagonal(WS_mat->W_in,WS_mat->W_in_vec);CHKERRQ(ierr);
			ierr = MatDestroy(&WS_mat->W_in);CHKERRQ(ierr);
		}       
	}

	ierr = PetscOptionsGetString(NULL,NULL,"-Wout",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,&WS_mat->flg_Wout);CHKERRQ(ierr);
	if (WS_mat->flg_Wout) {
		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);
		ierr = MatCreate(PETSC_COMM_WORLD,&WS_mat->W_out);CHKERRQ(ierr);
		ierr = MatSetType(WS_mat->W_out,MATMPIAIJ);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the output weight matrix!\n");
		ierr = MatLoad(WS_mat->W_out,fd);CHKERRQ(ierr);
		if (WS_mat->flg_diag) {
			ierr = MatGetSize(WS_mat->W_out,&N,NULL);CHKERRQ(ierr);
			ierr = VecCreate(PETSC_COMM_WORLD,&WS_mat->W_out_vec);CHKERRQ(ierr);
			ierr = VecSetSizes(WS_mat->W_out_vec,PETSC_DECIDE,N);CHKERRQ(ierr);
			ierr = VecSetUp(WS_mat->W_out_vec);CHKERRQ(ierr);
			ierr = MatGetDiagonal(WS_mat->W_out,WS_mat->W_out_vec);CHKERRQ(ierr);
			ierr = MatDestroy(&WS_mat->W_out);CHKERRQ(ierr);
		}       
	}

	ierr = PetscOptionsGetString(NULL,NULL,"-B",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,&WS_mat->flg_B);CHKERRQ(ierr);
	if (WS_mat->flg_B) {
		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);
		ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
		ierr = MatSetType(B,MATMPIAIJ);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the forcing spatial filter matrix!\n");
		ierr = MatLoad(B,fd);CHKERRQ(ierr);
		ierr = MatGetSize(B,&N,NULL);CHKERRQ(ierr);
		ierr = VecCreate(PETSC_COMM_WORLD,&WS_mat->B_vec);CHKERRQ(ierr);
		ierr = VecSetSizes(WS_mat->B_vec,PETSC_DECIDE,N);CHKERRQ(ierr);
		ierr = VecSetUp(WS_mat->B_vec);CHKERRQ(ierr);
		ierr = MatGetDiagonal(B,WS_mat->B_vec);CHKERRQ(ierr);
		ierr = MatDestroy(&B);CHKERRQ(ierr);
	}

	ierr = PetscOptionsGetString(NULL,NULL,"-C",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,&WS_mat->flg_C);CHKERRQ(ierr);
	if (WS_mat->flg_C) {
		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);
		ierr = MatCreate(PETSC_COMM_WORLD,&C);CHKERRQ(ierr);
		ierr = MatSetType(C,MATMPIAIJ);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the response spatial filter matrix!\n");
		ierr = MatLoad(C,fd);CHKERRQ(ierr);
		ierr = MatGetSize(C,&N,NULL);CHKERRQ(ierr);
		ierr = VecCreate(PETSC_COMM_WORLD,&WS_mat->C_vec);CHKERRQ(ierr);
		ierr = VecSetSizes(WS_mat->C_vec,PETSC_DECIDE,N);CHKERRQ(ierr);
		ierr = VecSetUp(WS_mat->C_vec);CHKERRQ(ierr);
		ierr = MatGetDiagonal(C,WS_mat->C_vec);CHKERRQ(ierr);
		ierr = MatDestroy(&C);CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
	
}

