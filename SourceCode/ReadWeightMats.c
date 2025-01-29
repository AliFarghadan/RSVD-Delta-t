
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode ReadWeightMats(RSVDt_vars *RSVDt, Weight_matrices *Weight, Directories *dirs)
{
	/*
		Reads the weight, input and outputs matrices if defined
		
		Assumptions:
		- B is a matrix of size N x Nb (Nb can equal N).
		- W_f_sqrt_inv is a matrix of size Nb x Nb.
		- C is a matrix of size Nc x N (Nc can equal N).
		- W_q_sqrt_inv and W_q_sqrt are matrices of size Nc x Nc.

		If any of these matrix dimensions mismatch the expected sizes, an error will occur.
		If the flag is off, an identity matrix is assumed.
	*/

	PetscErrorCode        ierr;
	PetscViewer           fd;
	PetscInt              row1, col1, row2, col2;

	PetscFunctionBeginUser;

	if (Weight->InvInputWeightFlg) {
		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->RootDir,dirs->InvInputWeightDir);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the inverse input weight matrix  : %s\n", dirs->IO_dir);
		ierr = MatCreate(PETSC_COMM_WORLD,&Weight->W_f_sqrt_inv);CHKERRQ(ierr);
		ierr = MatSetType(Weight->W_f_sqrt_inv,MATMPIAIJ);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = MatLoad(Weight->W_f_sqrt_inv,fd);CHKERRQ(ierr); 
		ierr = MatGetSize(Weight->W_f_sqrt_inv,&row1,&col1);CHKERRQ(ierr); 
		if (row1 != col1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Input weight matrix (W_f_sqrt_inv) must be square, current size = %d x %d", (int)row1, (int)col1);CHKERRQ(ierr);   
	} else {
		row1 = RSVDt->RSVD.N;
		col1 = RSVDt->RSVD.N;
	}

	if (Weight->InputMatrixFlg) {
		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->RootDir,dirs->InputMatrixDir);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the input matrix                 : %s\n", dirs->IO_dir);
		ierr = MatCreate(PETSC_COMM_WORLD,&Weight->B);CHKERRQ(ierr);
		ierr = MatSetType(Weight->B,MATMPIAIJ);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);	
		ierr = MatLoad(Weight->B,fd);CHKERRQ(ierr);
		ierr = MatGetSize(Weight->B,&row2,&col2);CHKERRQ(ierr); 
		if (col2 != col1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Size mismatch between input matrix (B) and input weight matrix (W_f_sqrt_inv), %d != %d", (int)col2, (int)col1);CHKERRQ(ierr);
		if (row2 != RSVDt->RSVD.N) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Input matrix (B) must have %d rows, current size = %d x %d", (int)RSVDt->RSVD.N, (int)row2, (int)col2);CHKERRQ(ierr);
		RSVDt->RSVD.Nb = col2;
	} else {
		RSVDt->RSVD.Nb = RSVDt->RSVD.N;
		if (Weight->InvInputWeightFlg && RSVDt->RSVD.Nb != col1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Size mismatch between input matrix (B) and input weight matrix (W_f_sqrt_inv)");CHKERRQ(ierr);
	}

	if (Weight->InvOutputWeightFlg) {
		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->RootDir,dirs->InvOutputWeightDir);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the inverse output weight matrix : %s\n", dirs->IO_dir);
		ierr = MatCreate(PETSC_COMM_WORLD,&Weight->W_q_sqrt_inv);CHKERRQ(ierr);
		ierr = MatSetType(Weight->W_q_sqrt_inv,MATMPIAIJ);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = MatLoad(Weight->W_q_sqrt_inv,fd);CHKERRQ(ierr);     
		ierr = MatGetSize(Weight->W_q_sqrt_inv,&row1,&col1);CHKERRQ(ierr); 
		if (row1 != col1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Inverse output weight matrix (W_q_sqrt_inv) must be square, current size = %d x %d", (int)row1, (int)col1);CHKERRQ(ierr); 
	}	

	if (Weight->OutputWeightFlg) {
		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->RootDir,dirs->OutputWeightDir);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the output weight matrix         : %s\n", dirs->IO_dir);
		ierr = MatCreate(PETSC_COMM_WORLD,&Weight->W_q_sqrt);CHKERRQ(ierr);
		ierr = MatSetType(Weight->W_q_sqrt,MATMPIAIJ);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = MatLoad(Weight->W_q_sqrt,fd);CHKERRQ(ierr);      
		ierr = MatGetSize(Weight->W_q_sqrt,&row1,&col1);CHKERRQ(ierr); 
		if (row1 != col1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Output weight matrix (W_q_sqrt) must be square, current size = %d x %d", (int)row1, (int)col1);CHKERRQ(ierr);
	} else {
		row1 = RSVDt->RSVD.N;
		col1 = RSVDt->RSVD.N;
	}

	if (Weight->OutputMatrixFlg) {
		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->RootDir,dirs->OutputMatrixDir);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the output matrix                : %s\n", dirs->IO_dir);
		ierr = MatCreate(PETSC_COMM_WORLD,&Weight->C);CHKERRQ(ierr);
		ierr = MatSetType(Weight->C,MATMPIAIJ);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = MatLoad(Weight->C,fd);CHKERRQ(ierr);
		ierr = MatGetSize(Weight->C,&row2,&col2);CHKERRQ(ierr); 
		if (col2 != RSVDt->RSVD.N) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Output matrix (C) must have %d columns, current size = %d x %d", (int)RSVDt->RSVD.N, (int)row2, (int)col2);CHKERRQ(ierr);
		RSVDt->RSVD.Nc = row2;
		if (row2 != row1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Size mismatch between output matrix (C) and output weight matrix (W_q_sqrt), %d != %d", (int)row2, (int)row1);CHKERRQ(ierr);
		if (Weight->InvOutputWeightFlg) {
			ierr = MatGetSize(Weight->W_q_sqrt_inv,&row1,NULL);CHKERRQ(ierr);
			if (row2 != row1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Size mismatch between output matrix (C) and inverse output weight matrix (W_q_sqrt_inv), %d != %d", (int)row2, (int)row1);CHKERRQ(ierr);
		}
	} else {
		RSVDt->RSVD.Nc = RSVDt->RSVD.N;
		if (Weight->OutputWeightFlg && RSVDt->RSVD.Nc != row1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Size mismatch between output matrix (C) and output weight matrix (W_q_sqrt)");CHKERRQ(ierr);
		if (Weight->InvOutputWeightFlg) ierr = MatGetSize(Weight->W_q_sqrt_inv,&row1,NULL);CHKERRQ(ierr);
		if (Weight->InvOutputWeightFlg && RSVDt->RSVD.Nc != row1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Size mismatch between output matrix (C) and inverse output weight matrix (W_q_sqrt_inv)");CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
	
}

