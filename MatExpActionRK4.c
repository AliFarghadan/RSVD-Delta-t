
#include <MatExpActionRK4.h>
#include <TSTransRK4.h>
#include <DestroyTSmat.h>

PetscErrorCode MatExpActionRK4(Mat Y_in, Mat Y_out, PetscInt Ngap, PetscInt st_ind, LNS_params *LNS_mat, \
					WS_matrices *WS_mat, TS_matrices *TS_mat, RSVDt_vars *RSVDt, Directories *dirs)
{

	PetscErrorCode       ierr;
	Vec                  y0,y_temp;
	PetscInt             N,k,i,ik,rend=st_ind+Ngap,prg_cnt=0;
	PetscViewer          fd;
	PetscLogDouble       t1,t2;

	PetscFunctionBeginUser;

	ierr = RSVDt->TS.flg_dir_adj ? PetscPrintf(PETSC_COMM_WORLD,"\n*** Computing the subspace evolution (direct action)! ***\n") : \
				PetscPrintf(PETSC_COMM_WORLD,"\n*** Computing the subspace evolution (adjoint action)! ***\n");CHKERRQ(ierr);

	/*
		Create and set the size of matrices/vecs
	*/

	ierr = MatGetSize(Y_in,&N,&k);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** time interval = %f , subspace dim = %d ***\n", (double)Ngap*RSVDt->TS.dt, (int)k);CHKERRQ(ierr);

	ierr = MatSetType(Y_out,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Y_out,PETSC_DECIDE,PETSC_DECIDE,N,k);CHKERRQ(ierr);
	ierr = MatSetUp(Y_out);CHKERRQ(ierr);	

	ierr = VecCreate(PETSC_COMM_WORLD,&y0);CHKERRQ(ierr);
	ierr = VecSetType(y0,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(y0,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = CreateRK4F(TS_mat,N);CHKERRQ(ierr);

	/*
		Loop over all test vectors
	*/

	for (ik=0; ik<k; ik++) {

		/*
			Extract IC's for each ik
		*/

		ierr = MatDenseGetColumnVecRead(Y_in,ik,&y_temp);CHKERRQ(ierr);
		ierr = VecCopy(y_temp,y0);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecRead(Y_in,ik,&y_temp);CHKERRQ(ierr);	

		/*
			Create initial LNS operator
		*/

		ierr = CreateStreamingLNS(LNS_mat,PetscFmodReal(2*(st_ind-1),2*RSVDt->TS.Nt),TS_mat->A3);CHKERRQ(ierr);

		/*
			Perform time-stepping
		*/

		ierr = PetscTime(&t1);CHKERRQ(ierr);

		for (i=st_ind; i<rend; i++) {

			ierr = TSTransRK4(LNS_mat,TS_mat,RSVDt,i,y0);CHKERRQ(ierr);

			if ((double)(i-st_ind)/(rend-st_ind) >= 0.01*prg_cnt){
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Progress percentage: %d% \n",(int)prg_cnt);CHKERRQ(ierr);
				ierr = PetscTime(&t2);CHKERRQ(ierr);
				ierr = PetscPrintf(PETSC_COMM_WORLD,"TS took %f seconds (wall-time)\n", t2-t1);CHKERRQ(ierr);		
				prg_cnt += 10;
			}
		}

		ierr = PetscTime(&t2);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Transient took %f seconds! (i = %d)\n", t2-t1, (int)ik+1);CHKERRQ(ierr);

		ierr = MatDenseGetColumnVecWrite(Y_out,ik,&y_temp);CHKERRQ(ierr);
		ierr = VecCopy(y0,y_temp);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(Y_out,ik,&y_temp);CHKERRQ(ierr);	
	}

	/*
		Save the output
	*/	

	ierr = PetscSNPrintf((char*)&dirs->file_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->root_dir,dirs->output,"EQ_test");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->file_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = MatView(Y_out,fd);CHKERRQ(ierr);

	/*
		Remove local variables
	*/

	ierr = VecDestroy(&y0);CHKERRQ(ierr);
	ierr = DestroyTSmat(TS_mat);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}
