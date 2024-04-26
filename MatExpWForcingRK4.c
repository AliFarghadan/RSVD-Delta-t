
#include <MatExpWForcingRK4.h>
#include <TSRK4.h>
#include <DestroyTSmat.h>

PetscErrorCode MatExpWForcingRK4(RSVD_matrices *RSVD_mat, DFT_matrices *DFT_mat, \
		LNS_params *LNS_mat, TS_matrices *TS_mat, RSVDt_vars *RSVDt, Directories *dirs)
{

	PetscErrorCode       ierr;
	Mat                  Y_all,F_hat_re,F_temp,Y_all_k,F_hat_k,Y_IC;
	Vec                  y,y_temp;
	PetscInt             N,k,Nt_saved,Nw_eff,i,ik,rend,prg_cnt,pos,mod,st_ind=1;
	PetscReal            norm;
	PetscBool            flg_IC;
	PetscViewer          fd;
	PetscLogDouble       t1,t2;

	PetscFunctionBeginUser;

	ierr = RSVDt->TS.flg_dir_adj ? PetscPrintf(PETSC_COMM_WORLD,"\n*** Direct action in process! ***\n") : \
				PetscPrintf(PETSC_COMM_WORLD,"\n*** Adjoint action in process! ***\n");CHKERRQ(ierr);

	/*
		Define TS parameters
	*/

	// rend  = RSVDt->TS.time_ratio;
	rend  = RSVDt->TS.Nt_transient;
	mod   = RSVDt->TS.time_ratio;
	Nt_saved = PetscFloorReal(rend/mod) + 1;
	rend += st_ind;

	/*
		Get the size of matrices
	*/

	ierr = MatGetSize(RSVD_mat->Y_hat,&N,&k);CHKERRQ(ierr);
	Nw_eff  = (RSVDt->TS.flg_real_A) ? RSVDt->RSVD.Nw/2 : RSVDt->RSVD.Nw;
	k   /= Nw_eff;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Total snapshots to be saved = %d for k = %d (The first snapshot is the I.C.)\n", (int)Nt_saved, (int)k);CHKERRQ(ierr);

	/*
		Apply solution from previous action as forcing
	*/

	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->F_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVD_mat->F_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(RSVD_mat->F_hat,PETSC_DECIDE,PETSC_DECIDE,N,Nw_eff*k);CHKERRQ(ierr);
	ierr = MatSetUp(RSVD_mat->F_hat);CHKERRQ(ierr);	
	ierr = MatCopy(RSVD_mat->Y_hat,RSVD_mat->F_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = MatDestroy(&RSVD_mat->Y_hat);CHKERRQ(ierr);		

	/*
		Create all required matrices (including solution, forcing, RHS, norm, etc.)
	*/

	ierr = VecCreate(PETSC_COMM_WORLD,&y);CHKERRQ(ierr);
	ierr = VecSetType(y,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(y,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = CreateRK4F(TS_mat,N);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&Y_all_k);CHKERRQ(ierr);
	ierr = MatSetType(Y_all_k,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Y_all_k,PETSC_DECIDE,PETSC_DECIDE,N,Nt_saved);CHKERRQ(ierr);
	ierr = MatSetUp(Y_all_k);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&Y_IC);CHKERRQ(ierr);
	ierr = MatSetType(Y_IC,MATDENSE);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL,NULL,"-IC",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,&flg_IC);CHKERRQ(ierr);
	if (flg_IC) {
		ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,dirs->filename);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);			
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading IC before forced run!\n");CHKERRQ(ierr);
		ierr = MatLoad(Y_IC,fd);CHKERRQ(ierr);
		ierr = MatGetSize(Y_IC,NULL,&k);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of test vectors = %d\n", (int)k);CHKERRQ(ierr);
	} else {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Starting from zero IC for all test vectors!\n");CHKERRQ(ierr);
		ierr = MatSetSizes(Y_IC,PETSC_DECIDE,PETSC_DECIDE,N,k);CHKERRQ(ierr);		
	}
	ierr = MatSetUp(Y_IC);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&Y_all);CHKERRQ(ierr);
	ierr = MatSetType(Y_all,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Y_all,PETSC_DECIDE,PETSC_DECIDE,N,Nt_saved*k);CHKERRQ(ierr);
	ierr = MatSetUp(Y_all);CHKERRQ(ierr);

	/*
		Reshape the forcing matrix
	*/

	ierr = MatCreate(PETSC_COMM_WORLD,&F_hat_re);CHKERRQ(ierr);
	ierr = MatSetType(F_hat_re,MATDENSE);CHKERRQ(ierr);
	ierr = PermuteForcing(RSVD_mat->F_hat, Nw_eff, F_hat_re);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&F_hat_k);CHKERRQ(ierr);
	ierr = MatSetType(F_hat_k,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(F_hat_k,PETSC_DECIDE,PETSC_DECIDE,N,Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(F_hat_k);CHKERRQ(ierr);

	/*
		Extract the forcing submatrix for each mode and loop-over all modes
	*/	

	for (ik=0; ik<k; ik++) {

		prg_cnt = 0;

		ierr = MatDenseGetColumnVecWrite(Y_IC,ik,&y_temp);CHKERRQ(ierr);
		ierr = VecCopy(y_temp,y);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(Y_IC,ik,&y_temp);CHKERRQ(ierr);

		ierr = MatDenseGetSubMatrix(F_hat_re,PETSC_DECIDE,PETSC_DECIDE,ik*Nw_eff,(ik+1)*Nw_eff,&F_temp);CHKERRQ(ierr);
		ierr = MatCopy(F_temp,F_hat_k,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		ierr = MatAssemblyBegin(F_hat_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(F_hat_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(F_hat_re,&F_temp);CHKERRQ(ierr);

		/*
			Create initial forcing and LNS operator
		*/

		ierr = RSVDt->TS.flg_dir_adj ? CreateStreamingForcing(F_hat_k,DFT_mat->inv_dft_dir,0,TS_mat->F3) : \
					CreateStreamingForcing(F_hat_k,DFT_mat->inv_dft_dir,0,TS_mat->F3);CHKERRQ(ierr);

		ierr = RSVDt->TS.flg_dir_adj ? CreateStreamingLNS(LNS_mat,0,TS_mat->A3) : \
				CreateStreamingLNS(LNS_mat,0,TS_mat->A3);CHKERRQ(ierr);

		/*
			Measure time-stepping wall-time
		*/

		ierr = PetscTime(&t1);CHKERRQ(ierr);

		/*
			Perform time-stepping (save "initial conditions" instead of "responses")
		*/

		for (i=st_ind; i<rend; i++) {

			if (PetscFmodReal(i-st_ind,mod) == 0) {
				pos  = (i-st_ind)/mod;
				ierr = MatDenseGetColumnVecWrite(Y_all_k,pos,&y_temp);CHKERRQ(ierr);
				ierr = VecCopy(y,y_temp);CHKERRQ(ierr);
				ierr = MatDenseRestoreColumnVecWrite(Y_all_k,pos,&y_temp);CHKERRQ(ierr);

				ierr = VecNorm(y,NORM_2,&norm);CHKERRQ(ierr);
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of the solution: %f @ position index = %d\n",(double)norm,(int)PetscFmodReal(pos,Nw_eff)+1);CHKERRQ(ierr);
			}

			ierr = TSRK4(F_hat_k,DFT_mat,LNS_mat,TS_mat,RSVDt,i,y);CHKERRQ(ierr);

			if ((double)(i-st_ind+1)/(rend-st_ind+1) >= 0.01*prg_cnt){
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Progress percentage: %d% \n",(int)prg_cnt);CHKERRQ(ierr);
				ierr = PetscTime(&t2);CHKERRQ(ierr);
				ierr = PetscPrintf(PETSC_COMM_WORLD,"TS took %f seconds (wall-time)\n", t2-t1);CHKERRQ(ierr);		
				prg_cnt += 10;
			}
		}

		ierr = PetscTime(&t2);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"TS took %f seconds (k = %d)\n", t2-t1, (int)ik);CHKERRQ(ierr);

		ierr = MatAssemblyBegin(Y_all_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(Y_all_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nt_saved,(ik+1)*Nt_saved,&F_temp);CHKERRQ(ierr);
		ierr = MatCopy(Y_all_k,F_temp,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Y_all,&F_temp);CHKERRQ(ierr);	

	}	

	ierr = RSVDt->TS.flg_dir_adj ? PetscPrintf(PETSC_COMM_WORLD,"Direct action DONE!\n") : \
				PetscPrintf(PETSC_COMM_WORLD,"Adjoint action DONE!\n");CHKERRQ(ierr);		

	/*
		Remove the forcing input as it is no longer needed
	*/

	ierr = VecDestroy(&y);CHKERRQ(ierr);
	ierr = MatDestroy(&F_hat_re);CHKERRQ(ierr);
	ierr = DestroyTSmat(TS_mat);CHKERRQ(ierr);

	/*
		Save the output
	*/

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"Y_forced_dir");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = MatView(Y_all,fd);CHKERRQ(ierr);			

	PetscFunctionReturn(0);
	
}

