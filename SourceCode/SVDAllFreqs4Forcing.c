
#include <slepcsvd.h>
#include <Variables.h>
#include <PermuteMat.h>

PetscErrorCode SVDAllFreqs4Forcing(RSVD_matrices *RSVD, RSVDt_vars *RSVDt, Weight_matrices *Weight, Resolvent_matrices *Res, Directories *dirs)
{
	/*
		Performs the economy SVD of matrices of size N \times k for Nw frequencies
	*/
	
	PetscErrorCode        ierr;
	PetscInt              iw, ik, hh, mm, ss;
	PetscReal             sigma;
	Mat                   Y, V;
	Vec                   v, Sig;
	SVD                   svd;
	PetscViewer           fd;
	PetscLogDouble        t1, t2;


	PetscFunctionBeginUser;

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** Reduced SVD begins! ***\n");CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&Res->V_hat);CHKERRQ(ierr);
	ierr = MatSetType(Res->V_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Res->V_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.Nb,RSVDt->RSVD.k*RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(Res->V_hat);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&Res->S_hat);CHKERRQ(ierr);
	ierr = MatSetType(Res->S_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Res->S_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.k,RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(Res->S_hat);CHKERRQ(ierr);

	for (iw=0; iw<RSVDt->RSVD.Nw_eff; iw++) {

		ierr = MatDenseGetSubMatrix(RSVD->Y_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&Y);CHKERRQ(ierr);
		ierr = SVDCreate(PETSC_COMM_WORLD,&svd);CHKERRQ(ierr);
		ierr = SVDSetOperators(svd,Y,NULL);CHKERRQ(ierr);
		ierr = SVDSetDimensions(svd,RSVDt->RSVD.k,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
		ierr = SVDSolve(svd);CHKERRQ(ierr);

		ierr = MatDenseGetSubMatrix(Res->V_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&V);CHKERRQ(ierr);
		ierr = MatDenseGetColumnVecWrite(Res->S_hat,iw,&Sig);CHKERRQ(ierr);
		for (ik=0; ik<RSVDt->RSVD.k; ik++) {
			ierr = MatDenseGetColumnVecWrite(V,ik,&v);CHKERRQ(ierr);
			ierr = SVDGetSingularTriplet(svd,ik,&sigma,v,NULL);CHKERRQ(ierr);
			ierr = VecSetValue(Sig,ik,sigma,INSERT_VALUES);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(V,ik,&v);CHKERRQ(ierr);
			ierr = VecAssemblyBegin(Sig);CHKERRQ(ierr);
			ierr = VecAssemblyEnd(Sig);CHKERRQ(ierr);
		}
		ierr = MatDenseRestoreColumnVecWrite(Res->S_hat,iw,&Sig);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Res->V_hat,&V);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(RSVD->Y_hat,&Y);CHKERRQ(ierr);
		ierr = SVDDestroy(&svd);CHKERRQ(ierr);
	}

	ierr = MatDestroy(&RSVD->Y_hat);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(Res->V_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Res->V_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	if (Weight->InvInputWeightFlg)  ierr = MatMatMult(Weight->W_f_sqrt_inv,Res->V_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&Res->V_hat);CHKERRQ(ierr);
	
	/*
		Prints out the elapsed time
	*/

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Reduced SVD elapsed time = %02d:%02d:%02d ***\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);

	/*
		Saves the forcing modes and gains
	*/

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** Saving the forcing modes and gains begins! ***\n");CHKERRQ(ierr);

	/*
		Saves resolvent modes for each mode separately (accross all frequencies) 
	*/

	ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->FolderDir,"S_hat");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = MatView(Res->S_hat,fd);CHKERRQ(ierr);

	/*
		Nw matrix of size N x k for response modes and the same for forcing modes (SaveResults == 1)
		OR 
		k matrix of size N x Nw for response modes and the same for forcing modes (SaveResults == 2)
	*/

	if (RSVDt->SaveResultsOpt == 1) {

		ierr = PermuteMat(Res->V_hat, RSVDt);CHKERRQ(ierr);

		for (ik=0; ik<RSVDt->RSVD.k; ik++) {

			ierr = MatDenseGetSubMatrix(Res->V_hat,PETSC_DECIDE,PETSC_DECIDE,ik*RSVDt->RSVD.Nw_eff,(ik+1)*RSVDt->RSVD.Nw_eff,&V);CHKERRQ(ierr);
			ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d%s",dirs->FolderDir,"V_hat_k",(int) ik+1,"_allFreqs");CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
			ierr = MatView(V,fd);CHKERRQ(ierr);
			ierr = MatDenseRestoreSubMatrix(Res->V_hat,&V);CHKERRQ(ierr);
		}
	} else { // RSVDt->SaveResultsOpt == 2

		for (iw=0; iw<RSVDt->RSVD.Nw_eff; iw++) {

			ierr = MatDenseGetSubMatrix(Res->V_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&V);CHKERRQ(ierr);
			ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d%s",dirs->FolderDir,"V_hat_Freq",(int) iw,"_allK");CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
			ierr = MatView(V,fd);CHKERRQ(ierr);
			ierr = MatDenseRestoreSubMatrix(Res->V_hat,&V);CHKERRQ(ierr);
		}
	}

	/*
		Prints out the elapsed time and exits
	*/

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Saving the forcing modes and gains elapsed time = %02d:%02d:%02d ***\n\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"One matrix of size %d (modes) x %d (frequencies) for gains\n", \
		(int) RSVDt->RSVD.k, (int) RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);
	if (RSVDt->SaveResultsOpt == 1) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"%d matrices of size %d x %d for response modes\n%d matrices of size %d x %d for forcing modes\n\n", \
				(int) RSVDt->RSVD.k, (int) RSVDt->RSVD.Nc, (int) RSVDt->RSVD.Nw_eff, \
				(int) RSVDt->RSVD.k, (int) RSVDt->RSVD.Nb, (int) RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);
	} else { // RSVDt->SaveResultsOpt == 2
		ierr = PetscPrintf(PETSC_COMM_WORLD,"%d matrices of size %d x %d for response modes\n%d matrices of size %d x %d for forcing modes\n\n", \
			(int) RSVDt->RSVD.Nw_eff, (int) RSVDt->RSVD.Nc, (int) RSVDt->RSVD.k, \
			(int) RSVDt->RSVD.Nw_eff, (int) RSVDt->RSVD.Nb, (int) RSVDt->RSVD.k);CHKERRQ(ierr);	
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"The results directory: %s\n\n",dirs->FolderDir);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}



