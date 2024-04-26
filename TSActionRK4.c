
#include <TSActionRK4.h>
#include <TSRK4.h>
#include <SaveSnapshots.h>
#include <DestroyTSmat.h>

PetscErrorCode TSActionRK4(RSVD_matrices *RSVD_mat, DFT_matrices *DFT_mat, \
		LNS_params *LNS_mat, TS_matrices *TS_mat, RSVDt_vars *RSVDt, Directories *dirs)
{

	PetscErrorCode       ierr;
	KSP                  ksp;
	Mat                  Y_all,F_hat_re,F_temp,Y_all_k,F_hat_k,Q_basis,EQ,M_tilde;
	Vec                  y,y_temp,y1,y2,b;
	PetscInt             N,k,Nw_output,Nw_input,i,ik,rend,prg_cnt=0,period_ind,pos;
	PetscReal            norm;
	PetscLogDouble       t1,t2;
	PetscViewer          fd;

	PetscFunctionBeginUser;

	ierr = RSVDt->TS.flg_dir_adj ? PetscPrintf(PETSC_COMM_WORLD,"\n*** Direct action in process! ***\n") : \
				PetscPrintf(PETSC_COMM_WORLD,"\n*** Adjoint action in process! ***\n");CHKERRQ(ierr);

	/*
		Define TS parameters
	*/

	k     = RSVDt->RSVD.k;			
	rend  = RSVDt->TS.Nt_transient + RSVDt->TS.Nt;
	rend += (RSVDt->TS.flg_eff_trans) ? RSVDt->TS.time_ratio : 0; 

	// Nw_eff  = (RSVDt->TS.flg_real_A) ? RSVDt->RSVD.Nw/2 : RSVDt->RSVD.Nw;
	Nw_output = (RSVDt->TS.flg_dir_adj) ? RSVDt->RSVD.Nw_out : RSVDt->RSVD.Nw;

	/*
		Get the size of matrices
	*/

	ierr = MatGetSize(RSVD_mat->Y_hat,&N,&Nw_input);CHKERRQ(ierr);
	Nw_input /= k;
	
	/*
		Apply solution from previous action as forcing
	*/

	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->F_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVD_mat->F_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(RSVD_mat->F_hat,PETSC_DECIDE,PETSC_DECIDE,N,Nw_input*k);CHKERRQ(ierr);
	ierr = MatSetUp(RSVD_mat->F_hat);CHKERRQ(ierr);	
	ierr = MatCopy(RSVD_mat->Y_hat,RSVD_mat->F_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = MatDestroy(&RSVD_mat->Y_hat);CHKERRQ(ierr);

	/*
		Create all required matrices (including solution, forcing, RHS, norm, etc.)
	*/

	ierr = VecCreate(PETSC_COMM_WORLD,&y);CHKERRQ(ierr);
	ierr = VecSetType(y,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(y,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&F_hat_re);CHKERRQ(ierr);
	ierr = MatSetType(F_hat_re,MATDENSE);CHKERRQ(ierr);

	ierr = CreateRK4Mats(TS_mat,N);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&Y_all_k);CHKERRQ(ierr);
	ierr = MatSetType(Y_all_k,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Y_all_k,PETSC_DECIDE,PETSC_DECIDE,N,Nw_output);CHKERRQ(ierr);
	ierr = MatSetUp(Y_all_k);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&Y_all);CHKERRQ(ierr);
	ierr = MatSetType(Y_all,MATDENSE);CHKERRQ(ierr);
	ierr = !RSVDt->TS.flg_eff_trans ? MatSetSizes(Y_all,PETSC_DECIDE,PETSC_DECIDE,N,Nw_output*k) : \
						MatSetSizes(Y_all,PETSC_DECIDE,PETSC_DECIDE,N,(Nw_output+1)*k);CHKERRQ(ierr);
	ierr = MatSetUp(Y_all);CHKERRQ(ierr);

	/*
		Reshape the forcing matrix
	*/

	ierr = PermuteMat(RSVD_mat->F_hat, Nw_input, F_hat_re);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&F_hat_k);CHKERRQ(ierr);
	ierr = MatSetType(F_hat_k,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(F_hat_k,PETSC_DECIDE,PETSC_DECIDE,N,Nw_input);CHKERRQ(ierr);
	ierr = MatSetUp(F_hat_k);CHKERRQ(ierr);

	/*
		Extract the forcing submatrix for each mode and loop-over all modes
	*/	

	for (ik=0; ik<k; ik++) {

		ierr = VecZeroEntries(y);CHKERRQ(ierr);

		ierr = MatDenseGetSubMatrix(F_hat_re,PETSC_DECIDE,PETSC_DECIDE,ik*Nw_input,(ik+1)*Nw_input,&F_temp);CHKERRQ(ierr);
		ierr = MatCopy(F_temp,F_hat_k,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		ierr = MatAssemblyBegin(F_hat_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(F_hat_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(F_hat_re,&F_temp);CHKERRQ(ierr);

		/*
			Create initial forcing and LNS operator
		*/

		ierr = RSVDt->TS.flg_dir_adj ? CreateStreamingForcing(F_hat_k,DFT_mat->inv_dft_dir,0,TS_mat->F3) : \
				CreateStreamingForcing(F_hat_k,DFT_mat->inv_dft_adj,0,TS_mat->F3);CHKERRQ(ierr);

		ierr = RSVDt->TS.flg_dir_adj ? CreateStreamingLNS(LNS_mat,0,TS_mat->A3) : \
				CreateStreamingLNS(LNS_mat,0,TS_mat->A3);CHKERRQ(ierr);

		/*
			Measure time-stepping wall-time
		*/

		ierr = PetscTime(&t1);CHKERRQ(ierr);

		/*
			Perform time-stepping
		*/

		if (RSVDt->TS.flg_update_IC) { // Updating I.C.

			ierr = PetscPrintf(PETSC_COMM_WORLD,"Updating I.C. is in process!\n");

			/*
				Read-in basis and evolved basis (Q and EQ)
			*/

			if (ik == 0) {							

				ierr = MatCreate(PETSC_COMM_WORLD,&Q_basis);CHKERRQ(ierr);
				ierr = MatSetType(Q_basis,MATDENSE);CHKERRQ(ierr);
				ierr = RSVDt->TS.flg_dir_adj ? PetscOptionsGetString(NULL,NULL,"-Q_basis_dir",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,NULL) : \
							PetscOptionsGetString(NULL,NULL,"-Q_basis_adj",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,NULL);CHKERRQ(ierr);
				ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);
				ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
				ierr = MatLoad(Q_basis,fd);CHKERRQ(ierr);
				ierr = MatCreateVecs(Q_basis,&b,&y1);CHKERRQ(ierr);
				ierr = MatCreateVecs(Q_basis,NULL,&y2);CHKERRQ(ierr);
				ierr = MatHermitianTranspose(Q_basis,MAT_INPLACE_MATRIX,&Q_basis);CHKERRQ(ierr);

				ierr = MatCreate(PETSC_COMM_WORLD,&EQ);CHKERRQ(ierr);
				ierr = MatSetType(EQ,MATDENSE);CHKERRQ(ierr);
				ierr = RSVDt->TS.flg_dir_adj ? PetscOptionsGetString(NULL,NULL,"-EQ_dir",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,NULL) : \
							PetscOptionsGetString(NULL,NULL,"-EQ_adj",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,NULL);CHKERRQ(ierr);
				ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->input,dirs->filename);CHKERRQ(ierr);
				ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
				ierr = MatLoad(EQ,fd);CHKERRQ(ierr);

				ierr = MatMatMult(Q_basis,EQ,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&M_tilde);CHKERRQ(ierr);
				ierr = MatShift(M_tilde, -1.);CHKERRQ(ierr);
				ierr = MatHermitianTranspose(Q_basis,MAT_INPLACE_MATRIX,&Q_basis);CHKERRQ(ierr);

				ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
				ierr = KSPSetOperators(ksp,M_tilde,M_tilde);CHKERRQ(ierr);
				ierr = KSPSetType(ksp,KSPPGMRES);CHKERRQ(ierr);
				ierr = KSPSetTolerances(ksp,1.e-15,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
			
			}	

			/*
				Updating I.C.
			*/

			period_ind = PetscFloorReal(RSVDt->TS.Nt_transient/RSVDt->TS.Nt);

			for (i=1; i<=period_ind*RSVDt->TS.Nt+1; i++) {

				ierr = TS_RK4(F_hat_k,DFT_mat,LNS_mat,TS_mat,RSVDt,i,y);CHKERRQ(ierr);

				if ((double)i/(period_ind*RSVDt->TS.Nt) >= 0.01*prg_cnt){
					ierr = PetscPrintf(PETSC_COMM_WORLD,"Progress percentage: %d% \n",(int)prg_cnt);CHKERRQ(ierr);
					ierr = PetscTime(&t2);CHKERRQ(ierr);
					ierr = PetscPrintf(PETSC_COMM_WORLD,"TS took %f seconds (wall-time)\n", t2-t1);CHKERRQ(ierr);		
					prg_cnt += 10;
				}

				if (i == (period_ind-1)*RSVDt->TS.Nt) ierr = VecCopy(y,y1);CHKERRQ(ierr);
				if (i == period_ind*RSVDt->TS.Nt) ierr = VecCopy(y,y2);CHKERRQ(ierr);
			}

			ierr = VecNorm(y1, NORM_2, &norm);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of y1 is %f\n", (double)norm);CHKERRQ(ierr);
			ierr = VecNorm(y2, NORM_2, &norm);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of y2 is %f\n", (double)norm);CHKERRQ(ierr);
			
			ierr = VecAXPY(y2,-1.,y1);CHKERRQ(ierr);
			ierr = MatMultHermitianTranspose(Q_basis, y2, b);CHKERRQ(ierr);

			ierr = KSPSolve(ksp,b,b);CHKERRQ(ierr);

			ierr = MatMult(Q_basis,b,y2);CHKERRQ(ierr);		
			ierr = VecAXPY(y1,-1.,y2);CHKERRQ(ierr);

			ierr = VecNorm(y1, NORM_2, &norm);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of y_IC is %f\n", (double)norm);CHKERRQ(ierr);				

			ierr = VecCopy(y1,y);CHKERRQ(ierr);

			/*
				One period to obtain steady-state snapshots
			*/
			
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Final period before saving steady-state snapshots!\n");
			if (ik == 0) prg_cnt = 0;

			for (i=1; i<=RSVDt->TS.Nt; i++) {

				ierr = TS_RK4(F_hat_k,DFT_mat,LNS_mat,TS_mat,RSVDt,i,y);CHKERRQ(ierr);

				if (PetscFmodReal(i,RSVDt->TS.time_ratio) == 0) {
					pos  = i/RSVDt->TS.time_ratio;
					pos  = PetscFmodReal(pos, Nw_output);
					pos  = (RSVDt->TS.flg_dir_adj) ? pos : Nw_output-pos;
					pos  = PetscFmodReal(pos, Nw_output);
					ierr = MatDenseGetColumnVecWrite(Y_all_k,pos,&y_temp);CHKERRQ(ierr);
					ierr = VecCopy(y,y_temp);CHKERRQ(ierr);
					ierr = MatDenseRestoreColumnVecWrite(Y_all_k,pos,&y_temp);CHKERRQ(ierr);
				}

				if ((double)i/RSVDt->TS.Nt >= 0.01*prg_cnt){
					ierr = PetscPrintf(PETSC_COMM_WORLD,"Progress percentage: %d% \n",(int)prg_cnt);CHKERRQ(ierr);
					ierr = PetscTime(&t2);CHKERRQ(ierr);
					ierr = PetscPrintf(PETSC_COMM_WORLD,"TS took %f seconds (wall-time)\n", t2-t1);CHKERRQ(ierr);		
					prg_cnt += 10;
				}	

			}

			ierr = VecNorm(y, NORM_2, &norm);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of y_IC(+T) is %f\n", (double)norm);CHKERRQ(ierr);

			if (ik == k-1) {
				ierr = VecDestroy(&y1);CHKERRQ(ierr);
				ierr = VecDestroy(&y2);CHKERRQ(ierr);
				ierr = VecDestroy(&b);CHKERRQ(ierr);
				ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
				ierr = MatDestroy(&M_tilde);CHKERRQ(ierr);
			}	

		} else { // Traditional time-stepping

			for (i=1; i<=rend; i++) {

				ierr = TSRK4(F_hat_k,DFT_mat,LNS_mat,TS_mat,RSVDt,i,y);CHKERRQ(ierr);

				ierr = SaveSnapshots(y,i,RSVDt,Y_all_k);CHKERRQ(ierr);

				// ierr = DisplayProgress(y,i,RSVDt,Y_all_k);CHKERRQ(ierr);

				if ((double)i/rend >= 0.01*prg_cnt){
					ierr = PetscPrintf(PETSC_COMM_WORLD,"Progress percentage: %d% \n",(int)prg_cnt);CHKERRQ(ierr);
					ierr = PetscTime(&t2);CHKERRQ(ierr);
					ierr = PetscPrintf(PETSC_COMM_WORLD,"TS took %f seconds (wall-time)\n", t2-t1);CHKERRQ(ierr);		
					prg_cnt += 10;
				}
			}

			ierr = PetscTime(&t2);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,"TS took %f seconds (k = %d)\n", t2-t1, (int)ik);CHKERRQ(ierr);

		}

		ierr = MatAssemblyBegin(Y_all_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(Y_all_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nw_output,(ik+1)*Nw_output,&F_temp);CHKERRQ(ierr);
		ierr = MatCopy(Y_all_k,F_temp,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Y_all,&F_temp);CHKERRQ(ierr);

	}

	/*
		Remove the forcing input (no longer needed)
	*/

	ierr = VecDestroy(&y);CHKERRQ(ierr);
	ierr = MatDestroy(&F_hat_re);CHKERRQ(ierr);
	ierr = DestroyTSmat(TS_mat);CHKERRQ(ierr);

	/*
		Back to frequency domain
	*/

	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->Y_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVD_mat->Y_hat,MATDENSE);CHKERRQ(ierr);

	ierr = RSVDt->TS.flg_dir_adj ? DFT(Y_all,DFT_mat->dft_dir,RSVDt->TS.time_ratio_out,RSVD_mat->Y_hat,RSVDt->TS.flg_real_A) : \
						DFT(Y_all,DFT_mat->dft_adj,RSVDt->TS.time_ratio,RSVD_mat->Y_hat,RSVDt->TS.flg_real_A);CHKERRQ(ierr);
	
	ierr = RSVDt->TS.flg_dir_adj ? PetscPrintf(PETSC_COMM_WORLD,"Direct action DONE!\n") : \
				PetscPrintf(PETSC_COMM_WORLD,"Adjoint action DONE!\n");CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
