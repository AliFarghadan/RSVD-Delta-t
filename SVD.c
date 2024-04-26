
#include <SVD.h>
#include <MergeFreqs.h>
#include <SplitFreqs.h>
#include <PermuteMat.h>
#include <RecoverU.h>

PetscErrorCode SVD(RSVD_matrices *RSVD_mat, resolvent_matrices *res_mat, RSVDt_vars *RSVDt)
{
	
	PetscErrorCode        ierr;
	SVD                   svd;
	Mat                   Y_hat_re, Q_hat_re, Y_hat_tall, V_hat_tall, U_hat_tall, Q_hat_tall, U_til;
	Vec                   U, V;
	PetscReal             sigma;
	PetscInt              ik, N, k, Nw_eff, loc_rows;

	PetscFunctionBeginUser;

	Nw_eff  = (RSVDt->TS.flg_real_A) ? RSVDt->RSVD.Nw/2 : RSVDt->RSVD.Nw;
	ierr = MatGetSize(RSVD_mat->Y_hat,&N,&k);CHKERRQ(ierr);
	k   /= Nw_eff;

	ierr = MatCreate(PETSC_COMM_WORLD,&U_til);CHKERRQ(ierr);
	ierr = MatSetType(U_til,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(U_til,PETSC_DECIDE,PETSC_DECIDE,k,k);CHKERRQ(ierr);
	ierr = MatSetUp(U_til);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&Y_hat_re);CHKERRQ(ierr);
	ierr = MatSetType(Y_hat_re,MATDENSE);CHKERRQ(ierr);
	ierr = PermuteMat(RSVD_mat->Y_hat, Nw_eff, Y_hat_re);CHKERRQ(ierr);

	ierr = MergeFreqs(Y_hat_re,Nw_eff,&Y_hat_tall);CHKERRQ(ierr);
	ierr = MatGetLocalSize(Y_hat_tall,&loc_rows,NULL);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&V_hat_tall);CHKERRQ(ierr);
	ierr = MatSetType(V_hat_tall,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(V_hat_tall,loc_rows,PETSC_DECIDE,N*Nw_eff,k);CHKERRQ(ierr);
	ierr = MatSetUp(V_hat_tall);CHKERRQ(ierr);	

	ierr = SVDCreate(PETSC_COMM_WORLD,&svd);CHKERRQ(ierr);
	ierr = SVDSetOperators(svd,Y_hat_tall,NULL);CHKERRQ(ierr);
	ierr = SVDSetDimensions(svd,k,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = SVDSolve(svd);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &res_mat->S_hat);CHKERRQ(ierr);
	ierr = VecSetSizes(res_mat->S_hat, PETSC_DECIDE, RSVDt->RSVD.k);CHKERRQ(ierr);
	ierr = VecSetUp(res_mat->S_hat);CHKERRQ(ierr);

	for (ik=0; ik<k; ik++){
		ierr = MatDenseGetColumnVecWrite(U_til,ik,&U);CHKERRQ(ierr);
		ierr = MatDenseGetColumnVecWrite(V_hat_tall,ik,&V);CHKERRQ(ierr);
		ierr = SVDGetSingularTriplet(svd,ik,&sigma,V,U);
		ierr = VecSetValue(res_mat->S_hat,ik,sigma,INSERT_VALUES);
		ierr = MatDenseRestoreColumnVecWrite(U_til,ik,&U);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(V_hat_tall,ik,&V);CHKERRQ(ierr);		
	}

	ierr = VecAssemblyBegin(res_mat->S_hat);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(res_mat->S_hat);CHKERRQ(ierr);
	
	ierr = SplitFreqs(V_hat_tall,Nw_eff,&res_mat->V_hat);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&Q_hat_re);CHKERRQ(ierr);
	ierr = MatSetType(Q_hat_re,MATDENSE);CHKERRQ(ierr);
	ierr = PermuteMat(RSVD_mat->Q_hat,RSVDt->RSVD.Nw_out,Q_hat_re);CHKERRQ(ierr);
	ierr = MergeFreqs(Q_hat_re,RSVDt->RSVD.Nw_out,&Q_hat_tall);CHKERRQ(ierr);
	ierr = RecoverU(U_til,Q_hat_tall,&U_hat_tall);CHKERRQ(ierr);
	ierr = SplitFreqs(U_hat_tall,RSVDt->RSVD.Nw_out,&res_mat->U_hat);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
	
}

