
#include <TransientRK4.h>
#include <TSTransRK4.h>
#include <DestroyTSMats.h>

PetscErrorCode TransientRK4(Mat Y_all, Mat trans_norm, PetscInt k_trans, PetscInt mod, PetscInt rseed, \
		LNS_params *LNS_mat, WS_matrices *WS_mat, TS_matrices *TS_mat, RSVDt_vars *RSVDt, PetscBool Y_all_flg, Directories *dirs)
{

	PetscErrorCode       ierr;
	Mat                  Y_IC_rand,Y_all_k,Y_temp;
	Vec                  y0,norm_vec,y_temp;
	PetscInt             Nt_saved,i,ik,rend,pos;
	PetscReal            norm_real,div_val=1e3,con_val=1e-16;
	PetscScalar          norm;
	PetscBool            flg_div=0,flg_IC;
	PetscRandom          r;
	PetscViewer          fd;
	PetscLogDouble       t1, t2;

	PetscFunctionBeginUser;

	ierr = RSVDt->TS.flg_dir_adj ? PetscPrintf(PETSC_COMM_WORLD,"\n*** Transient run (direct action) in process! ***\n") : \
				PetscPrintf(PETSC_COMM_WORLD,"\n*** Transient run (adjoint action) in process! ***\n");CHKERRQ(ierr);

	/*
		Define TS parameters
	*/

	rend  = RSVDt->TS.Nt_transient;
	mod   = (mod <= 10) ? mod : RSVDt->TS.time_ratio;
	// mod   = 129;
	// ierr = PetscPrintf(PETSC_COMM_WORLD,"mod = %d\n", (int)mod);CHKERRQ(ierr);
	Nt_saved = PetscFloorReal(rend/mod) + 1;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Total snapshots to be saved = %d\n", (int) Y_all_flg * (int) Nt_saved);CHKERRQ(ierr);

	/*
		Create and set the size of matrices/vecs
	*/

	if (Y_all_flg) {
		ierr = MatSetType(Y_all,MATDENSE);CHKERRQ(ierr);
		ierr = MatSetSizes(Y_all,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,Nt_saved*k_trans);CHKERRQ(ierr);
		ierr = MatSetUp(Y_all);CHKERRQ(ierr);

		ierr = MatCreate(PETSC_COMM_WORLD,&Y_all_k);CHKERRQ(ierr);
		ierr = MatSetType(Y_all_k,MATDENSE);CHKERRQ(ierr);
		ierr = MatSetSizes(Y_all_k,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,Nt_saved);CHKERRQ(ierr);
		ierr = MatSetUp(Y_all_k);CHKERRQ(ierr);
	}

	ierr = MatSetType(trans_norm,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(trans_norm,PETSC_DECIDE,PETSC_DECIDE,Nt_saved,k_trans);CHKERRQ(ierr);
	ierr = MatSetUp(trans_norm);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&y0);CHKERRQ(ierr);
	ierr = VecSetType(y0,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(y0,PETSC_DECIDE,RSVDt->RSVD.N);CHKERRQ(ierr);	

	if (k_trans > 1) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Random matrix of size %d x %d!\n", (int)RSVDt->RSVD.N, (int)k_trans);CHKERRQ(ierr);
		ierr = PetscRandomCreate(PETSC_COMM_WORLD,&r);CHKERRQ(ierr);
		ierr = PetscRandomSetSeed(r, rseed);CHKERRQ(ierr);
		ierr = PetscRandomSeed(r);CHKERRQ(ierr);		
		ierr = MatCreate(PETSC_COMM_WORLD,&Y_IC_rand);CHKERRQ(ierr);
		ierr = MatSetType(Y_IC_rand,MATDENSE);CHKERRQ(ierr);
		ierr = MatSetSizes(Y_IC_rand,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,k_trans);CHKERRQ(ierr);
		ierr = MatSetUp(Y_IC_rand);CHKERRQ(ierr);
		ierr = MatSetRandom(Y_IC_rand,r);CHKERRQ(ierr);
	} else {	
		ierr = PetscOptionsGetString(NULL,NULL,"-IC",(char*)&dirs->filename,PETSC_MAX_PATH_LEN,&flg_IC);CHKERRQ(ierr);
		if (flg_IC){
			ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,dirs->filename);CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);			
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading a vector!\n");CHKERRQ(ierr);
			ierr = VecLoad(y0,fd);CHKERRQ(ierr);
			if (WS_mat->flg_Win) {
				ierr = VecPointwiseMult(y_temp,WS_mat->W_in_vec,y0);CHKERRQ(ierr);
				ierr = VecDot(y_temp,y0,&norm);CHKERRQ(ierr);
				norm = PetscSqrtReal(PetscRealPart(norm));
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of test vector: %f\n",(double)norm);CHKERRQ(ierr);  
				// ierr = VecScale(y0,1/norm);CHKERRQ(ierr);
			} else {
				// ierr = VecNormalize(y0,NULL);CHKERRQ(ierr);
			}		
		} else {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Random vector!\n");CHKERRQ(ierr);
			ierr = PetscRandomCreate(PETSC_COMM_WORLD,&r);CHKERRQ(ierr);
			ierr = PetscRandomSetSeed(r, rseed);CHKERRQ(ierr);
			ierr = PetscRandomSeed(r);CHKERRQ(ierr);
			ierr = VecSetRandom(y0,r);CHKERRQ(ierr);
			if (WS_mat->flg_Win && WS_mat->flg_diag) {
				ierr = VecPointwiseMult(y_temp,WS_mat->W_in_vec,y0);CHKERRQ(ierr);
				ierr = VecDot(y_temp,y0,&norm);CHKERRQ(ierr);
				norm = PetscSqrtReal(PetscRealPart(norm));
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of test vector: %f\n",(double)norm);CHKERRQ(ierr);  
				ierr = VecScale(y0,1/norm);CHKERRQ(ierr);
			} else if (WS_mat->flg_Win) {
				ierr = MatMult(WS_mat->W_in,y0,y_temp);CHKERRQ(ierr);
				ierr = VecDot(y_temp,y0,&norm);CHKERRQ(ierr);
				norm = PetscSqrtReal(PetscRealPart(norm));
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of test vector: %f\n",(double)norm);CHKERRQ(ierr);  
				ierr = VecScale(y0,1/norm);CHKERRQ(ierr); 
			} else {
				ierr = VecNormalize(y0,NULL);CHKERRQ(ierr);
			}
		}
	}

	ierr = CreateRK4F(TS_mat,RSVDt->RSVD.N);CHKERRQ(ierr);

	/*
		Loop over all test vectors
	*/	

	// ierr = PetscTime(&t1);CHKERRQ(ierr);

	for (ik=0; ik<k_trans; ik++) {

		if (flg_div) break;
		pos = 0;

		/*
			Extract random IC's for each ik
		*/

		if (k_trans > 1) {
			ierr = MatDenseGetColumnVecRead(Y_IC_rand,ik,&y_temp);CHKERRQ(ierr);
			ierr = VecCopy(y_temp,y0);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecRead(Y_IC_rand,ik,&y_temp);CHKERRQ(ierr);
		}

		if (WS_mat->flg_Win && WS_mat->flg_diag) {
			ierr = VecPointwiseMult(y_temp,WS_mat->W_in_vec,y0);CHKERRQ(ierr);
			ierr = VecDot(y_temp,y0,&norm);CHKERRQ(ierr);
			norm = PetscSqrtReal(PetscRealPart(norm));
			ierr = VecScale(y0,1/norm);CHKERRQ(ierr);
		} else if (WS_mat->flg_Win) {
			ierr = MatMult(WS_mat->W_in,y0,y_temp);CHKERRQ(ierr);
			ierr = VecDot(y_temp,y0,&norm);CHKERRQ(ierr);
			norm_real = PetscSqrtReal(PetscRealPart(norm));
			ierr = VecScale(y0,1/norm);CHKERRQ(ierr);
		} else{
			ierr = VecNorm(y0,NORM_2,&norm_real);CHKERRQ(ierr);
			norm = norm_real;
			// ierr = VecScale(y0,1/norm);CHKERRQ(ierr);
			ierr = VecNorm(y0,NORM_2,&norm_real);CHKERRQ(ierr);
			norm = norm_real;			
		}

		ierr = VecCreate(PETSC_COMM_WORLD,&norm_vec);CHKERRQ(ierr);
		ierr = VecSetType(norm_vec,VECMPI);CHKERRQ(ierr);
		ierr = VecSetSizes(norm_vec,PETSC_DECIDE,Nt_saved);CHKERRQ(ierr);	
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of test vector after normalization: %f @ t = %d\n",(double)norm,(int)0);CHKERRQ(ierr);           
		ierr = VecSetValues(norm_vec,1,&pos,&norm,INSERT_VALUES);CHKERRQ(ierr);			

		if (Y_all_flg) {
			ierr = MatDenseGetColumnVecWrite(Y_all_k,ik,&y_temp);CHKERRQ(ierr);
			ierr = VecCopy(y0,y_temp);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(Y_all_k,ik,&y_temp);CHKERRQ(ierr);		
		}

		/*
			Create initial LNS operator
		*/

		ierr = CreateStreamingLNS(LNS_mat,0,TS_mat->A3);CHKERRQ(ierr);

		/*
			Perform time-stepping
		*/

		ierr = PetscTime(&t1);CHKERRQ(ierr);

		for (i=1; i<=rend; i++) {

			ierr = TSTransRK4(LNS_mat,TS_mat,RSVDt,i,y0);CHKERRQ(ierr);

			if (PetscFmodReal(i,mod) == 0) {
				pos  = i/mod;
				// pos  += 1;
				if (Y_all_flg) {
					ierr = MatDenseGetColumnVecWrite(Y_all_k,pos,&y_temp);CHKERRQ(ierr);
					ierr = VecCopy(y0,y_temp);CHKERRQ(ierr);
					ierr = MatDenseRestoreColumnVecWrite(Y_all_k,pos,&y_temp);CHKERRQ(ierr);
				}

				if (WS_mat->flg_Win && WS_mat->flg_diag) {
					ierr = VecPointwiseMult(y_temp,WS_mat->W_in_vec,y0);CHKERRQ(ierr);
					ierr = VecDot(y_temp,y0,&norm);CHKERRQ(ierr);
					norm_real = PetscSqrtReal(PetscRealPart(norm));
				} else if (WS_mat->flg_Win) {
					ierr = MatMult(WS_mat->W_in,y0,y_temp);CHKERRQ(ierr);
					ierr = VecDot(y_temp,y0,&norm);CHKERRQ(ierr);
					norm_real = PetscSqrtReal(PetscRealPart(norm));
				} else{
					ierr = VecNorm(y0,NORM_2,&norm_real);CHKERRQ(ierr);
					norm = norm_real;
				}
				
				if (norm_real > div_val) {
					PetscPrintf(PETSC_COMM_WORLD,"Diverged!\n");CHKERRQ(ierr); 
					break;
					flg_div = 1;
				}

				if (norm_real < con_val) {
					PetscPrintf(PETSC_COMM_WORLD,"Converged!\n");CHKERRQ(ierr); 
					break;
				}
				
				norm = norm_real;
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of test vector: %f @ t = %f\n",(double)norm,i*RSVDt->TS.dt);CHKERRQ(ierr);
				ierr = VecSetValues(norm_vec,1,&pos,&norm,INSERT_VALUES);CHKERRQ(ierr);
			}
		}

		// ierr = PetscTime(&t2);CHKERRQ(ierr);
		// ierr = PetscPrintf(PETSC_COMM_WORLD,"TS took %f seconds (k = %d)\n", t2-t1, (int)ik);CHKERRQ(ierr);

		ierr = VecAssemblyBegin(norm_vec);CHKERRQ(ierr);
		ierr = VecAssemblyEnd(norm_vec);CHKERRQ(ierr);
		ierr = MatDenseGetColumnVecWrite(trans_norm,ik,&y_temp);CHKERRQ(ierr);
		ierr = VecCopy(norm_vec,y_temp);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(trans_norm,ik,&y_temp);CHKERRQ(ierr);

		if (Y_all_flg) {
			ierr = MatAssemblyBegin(Y_all_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
			ierr = MatAssemblyEnd(Y_all_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
			ierr = MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nt_saved,(ik+1)*Nt_saved,&Y_temp);CHKERRQ(ierr);
			ierr = MatCopy(Y_all_k,Y_temp,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
			ierr = MatDenseRestoreSubMatrix(Y_all,&Y_temp);CHKERRQ(ierr);
		}

		ierr = PetscTime(&t2);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Transient took %f seconds!\n", t2-t1);CHKERRQ(ierr);
	
	}

	/*
		Save the outputs
	*/	

	if (Y_all_flg && !flg_div) {
		ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"H_Y_all");CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
		ierr = MatView(Y_all,fd);CHKERRQ(ierr);
	}

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"H_trans_norm");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = MatView(trans_norm,fd);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"H_last_snapshot");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = VecView(y0,fd);CHKERRQ(ierr);	

	/*
		Remove local variables
	*/

	ierr = VecDestroy(&y0);CHKERRQ(ierr);
	ierr = DestroyTSmat(TS_mat);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}
