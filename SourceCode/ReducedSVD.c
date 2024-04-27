
#include <slepcsvd.h>
#include "Variables.h"
#include "RecoverU.h"

PetscErrorCode ReducedSVD(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, WS_matrices *WS_mat, Resolvent_matrices *Res_mat)
{
	/*
		Performs the economy SVD of matrices of size N \times k for Nw frequencies
	*/
	
	PetscErrorCode        ierr;
	PetscInt              iw, ik, hh, mm, ss;
	PetscReal             sigma;
	Mat                   Y_temp, U_temp, V_temp, U_til;
	Vec                   U, V, Sig;
	SVD                   svd;
	PetscLogDouble        t1, t2;


	PetscFunctionBeginUser;

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** Reduced SVD begins! ***\n");CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&U_til);CHKERRQ(ierr);
	ierr = MatSetType(U_til,MATMPIDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(U_til,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.k,RSVDt->RSVD.k*RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(U_til);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&Res_mat->V_hat);CHKERRQ(ierr);
	ierr = MatSetType(Res_mat->V_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Res_mat->V_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,RSVDt->RSVD.k*RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(Res_mat->V_hat);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&Res_mat->S_hat);CHKERRQ(ierr);
	ierr = MatSetType(Res_mat->S_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Res_mat->S_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.k,RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(Res_mat->S_hat);CHKERRQ(ierr);

	for (iw=0; iw<RSVDt->RSVD.Nw_eff; iw++) {

		ierr = MatDenseGetSubMatrix(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&Y_temp);CHKERRQ(ierr);
		ierr = SVDCreate(PETSC_COMM_WORLD,&svd);CHKERRQ(ierr);
		ierr = SVDSetOperators(svd,Y_temp,NULL);CHKERRQ(ierr);
		ierr = SVDSetDimensions(svd,RSVDt->RSVD.k,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
		ierr = SVDSolve(svd);CHKERRQ(ierr);

		ierr = MatDenseGetSubMatrix(U_til,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&U_temp);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(Res_mat->V_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&V_temp);CHKERRQ(ierr);
		ierr = MatDenseGetColumnVecWrite(Res_mat->S_hat,iw,&Sig);CHKERRQ(ierr);
		for (ik=0; ik<RSVDt->RSVD.k; ik++) {
			ierr = MatDenseGetColumnVecWrite(U_temp,ik,&U);CHKERRQ(ierr);
			ierr = MatDenseGetColumnVecWrite(V_temp,ik,&V);CHKERRQ(ierr);
			ierr = SVDGetSingularTriplet(svd,ik,&sigma,V,U);
			ierr = VecSetValue(Sig,ik,sigma,INSERT_VALUES);
			ierr = MatDenseRestoreColumnVecWrite(U_temp,ik,&U);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(V_temp,ik,&V);CHKERRQ(ierr);
			ierr = VecAssemblyBegin(Sig);CHKERRQ(ierr);
			ierr = VecAssemblyEnd(Sig);CHKERRQ(ierr);
		}
		ierr = MatDenseRestoreColumnVecWrite(Res_mat->S_hat,iw,&Sig);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(U_til,&U_temp);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Res_mat->V_hat,&V_temp);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(RSVD_mat->Y_hat,&Y_temp);CHKERRQ(ierr);
		ierr = SVDDestroy(&svd);CHKERRQ(ierr);
	}

	ierr = MatDestroy(&RSVD_mat->Y_hat);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(Res_mat->V_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Res_mat->V_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(U_til,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(U_til,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatConjugate(U_til);CHKERRQ(ierr);

	ierr = RecoverU(U_til,RSVDt,RSVD_mat,Res_mat);CHKERRQ(ierr);

	if (WS_mat->flg_Wout) {
		if (WS_mat->flg_diag) {
			ierr = MatDiagonalScale(Res_mat->U_hat,WS_mat->W_out_vec,NULL);CHKERRQ(ierr);
			ierr = MatDiagonalScale(Res_mat->V_hat,WS_mat->W_out_vec,NULL);CHKERRQ(ierr);
		}
		else {
			ierr = MatMatMult(WS_mat->W_out,Res_mat->U_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&Res_mat->U_hat);CHKERRQ(ierr);
			ierr = MatMatMult(WS_mat->W_out,Res_mat->V_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&Res_mat->V_hat);CHKERRQ(ierr);
		}
	}

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



