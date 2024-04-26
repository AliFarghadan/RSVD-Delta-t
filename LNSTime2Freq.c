
#include <LNSTime2Freq.h>

PetscErrorCode LNSTime2Freq(LNS_params *LNS_mat, Directories *dirs)
{

	PetscErrorCode      ierr;
	const PetscScalar  *vals, *A_array;
	PetscScalar        *A_temp_arr, *A_mean_arr;
	const PetscInt     *ia, *ja, *ia_mat, *ja_mat;
	PetscReal           norm;
	Mat                 A_temp, A_snap, A_mean;
	Vec                 V1, V2, A_vec;
	MatInfo             info;
	PetscInt            N,i,ip,rstart,rend,nnz_pr,i_col,count,loc_row;
	PetscViewer         fd;
	PetscLogDouble      t1, t2;

	PetscFunctionBeginUser;

	/*
		Read in LNS operators in Fourier space if true
	*/

	if (LNS_mat->flg_read_A_hat) {
		ierr = PetscOptionsGetString(NULL,NULL,"-A_hat_filename",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,NULL);CHKERRQ(ierr);
		ierr = MatCreate(PETSC_COMM_WORLD,&LNS_mat->A_hat);CHKERRQ(ierr);
		ierr = MatSetType(LNS_mat->A_hat,MATDENSE);CHKERRQ(ierr);
		ierr = MatSetSizes(LNS_mat->A_hat,LNS_mat->nz_u,PETSC_DECIDE,PETSC_DETERMINE,LNS_mat->RSVDt.LNS.Nb);CHKERRQ(ierr);
		ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscTime(&t1);CHKERRQ(ierr);
		ierr = MatLoad(LNS_mat->A_hat,fd);CHKERRQ(ierr);
		ierr = PetscTime(&t2);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the dense A_hat with Nb = %d took %f seconds!\n",(int)LNS_mat->RSVDt.LNS.Nb,t2-t1);
		// ierr = MatGetSize(LNS_mat->A_hat,&N,&k);CHKERRQ(ierr);
		// ierr = PetscPrintf(PETSC_COMM_WORLD,"Size of A_hat = %d x %d\n", (int)N, (int)k);CHKERRQ(ierr);
	} else {

		/*
			Read in LNS operators in time and compute Nb A_hat coefficients
		*/

		ierr = PetscTime(&t1);CHKERRQ(ierr);

		ierr = MatCreate(PETSC_COMM_WORLD,&A_snap);CHKERRQ(ierr);
		ierr = MatSetType(A_snap,MATDENSE);CHKERRQ(ierr);

		for (ip = 1; ip <= LNS_mat->RSVDt.LNS.Np; ip++) {

			ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s%d",dirs->root_dir,dirs->input,dirs->prename,(int)ip);CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);

			ierr = MatCreate(PETSC_COMM_WORLD,&A_temp);CHKERRQ(ierr);
			ierr = MatSetType(A_temp,MATMPIAIJ);CHKERRQ(ierr);
			ierr = MatLoad(A_temp,fd);CHKERRQ(ierr);
			if (LNS_mat->RSVDt.disc.flg_disc) ierr = MatShift(A_temp,-LNS_mat->RSVDt.disc.beta);CHKERRQ(ierr);

			if (ip == 1) {
				ierr = MatGetSize(A_temp,&N,NULL);CHKERRQ(ierr);
				ierr = MatSetSizes(A_snap,LNS_mat->nz_u,PETSC_DECIDE,PETSC_DETERMINE,LNS_mat->RSVDt.LNS.Np);CHKERRQ(ierr);
				ierr = MatSetUp(A_snap);CHKERRQ(ierr);
			}

			ierr = PetscMalloc1(LNS_mat->nz_u, &A_temp_arr);CHKERRQ(ierr);

			ierr = MatGetOwnershipRange(A_temp, &rstart, &rend);CHKERRQ(ierr);
			count = 0;
			for (i = rstart; i < rend; i++) {
				ierr = MatGetRow(A_temp, i, &nnz_pr, NULL, &vals);CHKERRQ(ierr);
				for (i_col = 0; i_col < nnz_pr; i_col++){
					A_temp_arr[count + i_col] = vals[i_col];
				}
				count += nnz_pr;
				ierr = MatRestoreRow(A_temp, i, &nnz_pr, NULL, &vals);CHKERRQ(ierr);
			}

			ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, LNS_mat->nz_u, PETSC_DECIDE, A_temp_arr, &V1);CHKERRQ(ierr);
			ierr = MatDenseGetColumnVecWrite(A_snap,ip-1,&V2);CHKERRQ(ierr);
			ierr = VecCopy(V1,V2);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(A_snap,ip-1,&V2);CHKERRQ(ierr);

			ierr = MatDestroy(&A_temp);CHKERRQ(ierr);
			ierr = VecDestroy(&V1);CHKERRQ(ierr);
			ierr = PetscFree(A_temp_arr);CHKERRQ(ierr);
		}

		if (LNS_mat->RSVDt.disc.flg_disc) ierr = PetscPrintf(PETSC_COMM_WORLD, \
				"**** Discounting with beta = %g ****\n", LNS_mat->RSVDt.disc.beta);CHKERRQ(ierr);

		ierr = MatAssemblyBegin(A_snap,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(A_snap,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

		ierr = PetscTime(&t2);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"%d LNS operators of size N = %d are loaded in %f seconds\n", (int)LNS_mat->RSVDt.LNS.Np, (int)N, t2-t1);CHKERRQ(ierr);

		ierr = MatMatMult(A_snap,LNS_mat->LNS_DFT.dft_LNS,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&LNS_mat->A_hat);CHKERRQ(ierr);

		if (LNS_mat->flg_save_A_hat_den_matrix) {
			ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s%d",dirs->root_dir,dirs->output,"A_hat_dense_Nb_",(int)LNS_mat->RSVDt.LNS.Nb);CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
			ierr = MatView(LNS_mat->A_hat,fd);CHKERRQ(ierr);
		}

		if (LNS_mat->flg_save_A_hat_sp_matrix) {
			for (ip = 0; ip < LNS_mat->RSVDt.LNS.Nb; ip++) {
				ierr = MatDenseGetColumnVecRead(LNS_mat->A_hat,ip,&A_vec);CHKERRQ(ierr);
				ierr = VecGetArrayRead(A_vec, &A_array);CHKERRQ(ierr);
				ierr = MatUpdateMPIAIJWithArray(LNS_mat->A, A_array);CHKERRQ(ierr);
				ierr = VecRestoreArrayRead(A_vec, &A_array);CHKERRQ(ierr);	
				ierr = MatDenseRestoreColumnVecRead(LNS_mat->A_hat,ip,&A_vec);CHKERRQ(ierr);		
				ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s%d",dirs->root_dir,dirs->output,"A_hat",(int)ip);CHKERRQ(ierr);
				ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
				ierr = MatView(LNS_mat->A,fd);CHKERRQ(ierr);
			}
		}

		if (LNS_mat->flg_Frob_norm) {		
			for (ip = 0; ip < LNS_mat->RSVDt.LNS.Nb; ip++) {
				ierr = MatDenseGetColumnVecRead(LNS_mat->A_hat,ip,&A_vec);CHKERRQ(ierr);
				ierr = VecNorm(A_vec, NORM_2, &norm);CHKERRQ(ierr);
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Frobenius norm of iw = %d is %f\n", (int)ip, (double) norm);CHKERRQ(ierr);
				ierr = MatDenseRestoreColumnVecRead(LNS_mat->A_hat,ip,&A_vec);CHKERRQ(ierr);	
			}
		}
	}

	PetscFunctionReturn(0);
	
}
