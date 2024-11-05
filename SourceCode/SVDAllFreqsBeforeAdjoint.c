
#include <slepcsvd.h>
#include <Variables.h>

PetscErrorCode SVDAllFreqsBeforeAdjoint(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt)
{
	/*
		Performs the economy SVD of matrices of size N \times k for Nw frequencies
		We perform SVD instead of QR to obtain U 
		This is generally more accurate than performing QR and recovering it later
	*/
	
	PetscErrorCode        ierr;
	PetscInt              iw, ik, hh, mm, ss;
	Mat                   Y, Y1, Y2;
	Vec                   U;
	SVD                   svd;
	PetscLogDouble        t1, t2;


	PetscFunctionBeginUser;

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** Reduced SVD begins! ***\n");CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&Y);CHKERRQ(ierr);
	ierr = MatSetType(Y,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Y,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.Nc,RSVDt->RSVD.k*RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(Y);CHKERRQ(ierr);

	for (iw=0; iw<RSVDt->RSVD.Nw_eff; iw++) {

		ierr = MatDenseGetSubMatrix(Y,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&Y1);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&Y2);CHKERRQ(ierr);
		ierr = SVDCreate(PETSC_COMM_WORLD,&svd);CHKERRQ(ierr);
		ierr = SVDSetOperators(svd,Y2,NULL);CHKERRQ(ierr);
		ierr = SVDSetDimensions(svd,RSVDt->RSVD.k,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
		ierr = SVDSolve(svd);CHKERRQ(ierr);

		for (ik=0; ik<RSVDt->RSVD.k; ik++) {
			ierr = MatDenseGetColumnVecWrite(Y1,ik,&U);CHKERRQ(ierr);
			ierr = SVDGetSingularTriplet(svd,ik,NULL,U,NULL);
			ierr = MatDenseRestoreColumnVecWrite(Y1,ik,&U);CHKERRQ(ierr);
		}
		ierr = MatDenseRestoreSubMatrix(RSVD_mat->Y_hat,&Y2);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Y,&Y1);CHKERRQ(ierr);
		
		ierr = SVDDestroy(&svd);CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(Y,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Y,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatCopy(Y,RSVD_mat->Y_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = MatDestroy(&Y);CHKERRQ(ierr);

	/*
		Printing out the elapsed time and exit
	*/

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Reduced SVD elapsed time = %02d:%02d:%02d ***\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}



