
#include <LNSStructure.h>

PetscErrorCode LNSStructure(LNS_params *LNS_mat, TS_matrices *TS_mat, Directories *dirs)
{

	PetscErrorCode      ierr;
	Mat                 A_org;
	MatInfo             info;
	PetscInt            N,loc_row;
	const PetscInt     *ia, *ja, *ia_mat, *ja_mat;
	PetscLogDouble      t1, t2;
	PetscViewer         fd;

	PetscFunctionBeginUser;

	ierr = PetscTime(&t1);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s%d",dirs->root_dir,dirs->input,dirs->prename,1);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&A_org);CHKERRQ(ierr);
	ierr = MatSetType(A_org,MATMPIAIJ);CHKERRQ(ierr);
	ierr = MatLoad(A_org,fd);CHKERRQ(ierr);
	ierr = MatGetSize(A_org,&N,NULL);CHKERRQ(ierr);

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"The matrix is loaded in %f seconds (N = %d)\n", t2-t1, (int)N);CHKERRQ(ierr);	

	ierr = MatGetInfo(A_org,MAT_LOCAL,&info);CHKERRQ(ierr);
	LNS_mat->nz_u = info.nz_used;

	ierr = MatGetLocalSize(A_org, &loc_row, NULL);CHKERRQ(ierr);

	ierr = PetscMalloc1(loc_row+1, &ia_mat);CHKERRQ(ierr);
	ierr = PetscMalloc1(LNS_mat->nz_u, &ja_mat);CHKERRQ(ierr);

	ierr = MatGetRowIJ(A_org, 0, PETSC_FALSE, PETSC_FALSE, NULL, &ia, &ja, NULL);CHKERRQ(ierr); 
	ierr = PetscArraycpy((PetscInt*)ia_mat, ia, loc_row+1);CHKERRQ(ierr);	
	ierr = PetscArraycpy((PetscInt*)ja_mat, ja, LNS_mat->nz_u);CHKERRQ(ierr);
	ierr = MatRestoreRowIJ(A_org, 0, PETSC_FALSE, PETSC_FALSE, NULL, &ia, &ja, NULL);CHKERRQ(ierr);

	ierr = MatDestroy(&A_org);CHKERRQ(ierr);

	ierr = PetscTime(&t1);CHKERRQ(ierr);	
	ierr = MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, loc_row, loc_row, N, N, ia_mat, ja_mat, NULL, &LNS_mat->A);CHKERRQ(ierr);
	ierr = MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, loc_row, loc_row, N, N, ia_mat, ja_mat, NULL, &TS_mat->A1);CHKERRQ(ierr);
	ierr = MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, loc_row, loc_row, N, N, ia_mat, ja_mat, NULL, &TS_mat->A2);CHKERRQ(ierr);
	ierr = MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, loc_row, loc_row, N, N, ia_mat, ja_mat, NULL, &TS_mat->A3);CHKERRQ(ierr);
	ierr = PetscTime(&t2);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"MatCreateMPIAIJWithArrays takes %f seconds\n", t2-t1);CHKERRQ(ierr);	

	ierr = VecCreate(PETSC_COMM_WORLD, &LNS_mat->A_nnz);CHKERRQ(ierr);
	ierr = VecSetSizes(LNS_mat->A_nnz, LNS_mat->nz_u, PETSC_DETERMINE);CHKERRQ(ierr);
	ierr = VecSetUp(LNS_mat->A_nnz);CHKERRQ(ierr);

	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

