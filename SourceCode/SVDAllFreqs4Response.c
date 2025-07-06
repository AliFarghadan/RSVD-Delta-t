
#include <slepcsvd.h>
#include <Variables.h>
#include <PermuteMat.h>

PetscErrorCode SVDAllFreqs4Response(RSVD_matrices *RSVD, RSVDt_vars *RSVDt, Weight_matrices *Weight, Resolvent_matrices *Res, Directories *dirs)
{
	/*
		Performs the economy SVD of matrices of size N \times k for Nw frequencies
		We perform SVD instead of QR to obtain U 
		This is generally more accurate than performing QR and recovering it later
		We save the response modes onto disk in this function to avoid keeping a large N x k x Nw matrix in memory
	*/
	
	PetscErrorCode        ierr;
	PetscInt              iw, ik, hh, mm, ss;
	Mat                   Y, Y1, Y2;
	Vec                   U;
	SVD                   svd;
	PetscViewer           fd;
	PetscLogDouble        t1, t2;


	PetscFunctionBeginUser;

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** Reduced SVD begins! ***\n");CHKERRQ(ierr);

	ierr = MatDuplicate(RSVD->Y_hat,MAT_DO_NOT_COPY_VALUES,&Y);CHKERRQ(ierr);

	for (iw=0; iw<RSVDt->RSVD.Nw_eff; iw++) {

		ierr = MatDenseGetSubMatrix(Y,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&Y1);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(RSVD->Y_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&Y2);CHKERRQ(ierr);
		ierr = SVDCreate(PETSC_COMM_WORLD,&svd);CHKERRQ(ierr);
		ierr = SVDSetOperators(svd,Y2,NULL);CHKERRQ(ierr);
		ierr = SVDSetDimensions(svd,RSVDt->RSVD.k,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
		ierr = SVDSolve(svd);CHKERRQ(ierr);

		for (ik=0; ik<RSVDt->RSVD.k; ik++) {
			ierr = MatDenseGetColumnVecWrite(Y1,ik,&U);CHKERRQ(ierr);
			ierr = SVDGetSingularTriplet(svd,ik,NULL,U,NULL);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(Y1,ik,&U);CHKERRQ(ierr);
		}
		ierr = MatDenseRestoreSubMatrix(RSVD->Y_hat,&Y2);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Y,&Y1);CHKERRQ(ierr);
		
		ierr = SVDDestroy(&svd);CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(Y,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Y,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatCopy(Y,RSVD->Y_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = MatDestroy(&Y);CHKERRQ(ierr);

	/*
		Prints out the elapsed time
	*/

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Reduced SVD elapsed time = %02d:%02d:%02d ***\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);	

	/*
		Saves the response modes
	*/

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** Saving the response modes begins! ***\n");CHKERRQ(ierr);

	if (Weight->InvOutputWeightFlg) {
		ierr = MatMatMult(Weight->W_q_sqrt_inv,RSVD->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Res->U_hat);CHKERRQ(ierr);
	} else {
		ierr = MatDuplicate(RSVD->Y_hat,MAT_COPY_VALUES,&Res->U_hat);CHKERRQ(ierr);
	}

	if (RSVDt->SaveResultsOpt == 1) {

		ierr = PermuteMat(Res->U_hat, RSVDt);CHKERRQ(ierr);

		for (ik=0; ik<RSVDt->RSVD.k; ik++) {
			
			ierr = MatDenseGetSubMatrix(Res->U_hat,PETSC_DECIDE,PETSC_DECIDE,ik*RSVDt->RSVD.Nw_eff,(ik+1)*RSVDt->RSVD.Nw_eff,&Y1);CHKERRQ(ierr);
			ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d%s",dirs->FolderDir,"U_hat_k",(int) ik+1,"_allFreqs");CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
			ierr = MatView(Y1,fd);CHKERRQ(ierr);
			ierr = MatDenseRestoreSubMatrix(Res->U_hat,&Y1);CHKERRQ(ierr);

		}
	} else { // RSVDt->SaveResultsOpt == 2

		for (iw=0; iw<RSVDt->RSVD.Nw_eff; iw++) {
			
			ierr = MatDenseGetSubMatrix(Res->U_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&Y1);CHKERRQ(ierr);
			ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d%s",dirs->FolderDir,"U_hat_Freq",(int) iw,"_allK");CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
			ierr = MatView(Y1,fd);CHKERRQ(ierr);
			ierr = MatDenseRestoreSubMatrix(Res->U_hat,&Y1);CHKERRQ(ierr);

		}
	}

	ierr = MatDestroy(&Res->U_hat);CHKERRQ(ierr);

	/*
		Prints out the elapsed time and exits
	*/

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Saving the response modes elapsed time = %02d:%02d:%02d ***\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);	

	PetscFunctionReturn(0);

}



