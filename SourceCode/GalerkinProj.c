
#include <petscksp.h>
#include <Variables.h>

PetscErrorCode GalerkinProj(Mat Y_all, Mat M_tilde, Mat V, PetscInt ik, \
							RSVDt_vars *RSVDt, TS_removal_matrices *TSR, DFT_matrices *DFT_mat)
{
	/*
		Performs Galerkin projection to remove the transient response
	*/

	PetscErrorCode        ierr;
	KSP                   ksp;
	Mat                   Y_all_k,M_temp,Y_temp;
	Vec                   y1,y2,dft_iw,y_iw,b,yt;
	PetscInt              Nstore,iw;
	PetscReal             deltaT = RSVDt->TS.dt*RSVDt->TS.ResRatio;

	PetscFunctionBeginUser;

	Nstore = RSVDt->RSVD.Nw + 1;
	ierr = MatCreateVecs(V,&b,&yt);CHKERRQ(ierr);
	ierr = MatCreateVecs(V,NULL,&y1);CHKERRQ(ierr);
	ierr = MatCreateVecs(V,NULL,&y2);CHKERRQ(ierr);

	for (iw=0; iw<RSVDt->RSVD.Nw_eff; iw++) {

		ierr = MatDuplicate(M_tilde,MAT_COPY_VALUES,&M_temp);CHKERRQ(ierr);

		ierr = RSVDt->TS.DirAdj ? MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nstore,(ik+1)*Nstore-1,&Y_all_k) : \
				MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nstore+1,(ik+1)*Nstore,&Y_all_k);CHKERRQ(ierr);
		ierr = MatDenseGetColumnVecRead(DFT_mat->dft,iw,&dft_iw);CHKERRQ(ierr);
		ierr = MatDuplicate(Y_all_k,MAT_COPY_VALUES,&Y_temp);CHKERRQ(ierr);
		ierr = MatMult(Y_temp,dft_iw,y1);CHKERRQ(ierr);
		ierr = VecScale(y1, RSVDt->TS.ResRatio);CHKERRQ(ierr);
		ierr = MatDestroy(&Y_temp);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Y_all,&Y_all_k);CHKERRQ(ierr);

		ierr = RSVDt->TS.DirAdj ? MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nstore+1,(ik+1)*Nstore,&Y_all_k) : \
				MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nstore,(ik+1)*Nstore-1,&Y_all_k);CHKERRQ(ierr);
		ierr = MatDuplicate(Y_all_k,MAT_COPY_VALUES,&Y_temp);CHKERRQ(ierr);
		ierr = MatMult(Y_temp,dft_iw,y2);CHKERRQ(ierr);
		ierr = VecScale(y2, RSVDt->TS.ResRatio);CHKERRQ(ierr);
		ierr = MatDestroy(&Y_temp);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecRead(DFT_mat->dft,iw,&dft_iw);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Y_all,&Y_all_k);CHKERRQ(ierr);

		ierr = MatScale(M_temp, -1.);CHKERRQ(ierr);
		ierr = RSVDt->TS.DirAdj ? MatShift(M_temp,PetscExpComplex(PETSC_i * RSVDt->RSVD.w * iw * deltaT)) : \
			MatShift(M_temp,PetscExpComplex(-PETSC_i * RSVDt->RSVD.w * iw * deltaT));CHKERRQ(ierr);
		ierr = VecScale(y2, -1.);CHKERRQ(ierr);
		ierr = VecCopy(y1, yt);CHKERRQ(ierr);
		ierr = RSVDt->TS.DirAdj ? VecAYPX(yt,PetscExpComplex(PETSC_i * RSVDt->RSVD.w * iw * deltaT),y2) : \
				VecAYPX(yt,PetscExpComplex(-PETSC_i * RSVDt->RSVD.w * iw * deltaT),y2);CHKERRQ(ierr);
		ierr = MatMultHermitianTranspose(V,yt,b);CHKERRQ(ierr);

		ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
		ierr = KSPSetOperators(ksp,M_temp,M_temp);CHKERRQ(ierr);
		ierr = KSPSetType(ksp,KSPPGMRES);CHKERRQ(ierr);
		ierr = KSPSetTolerances(ksp,1.e-15,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
		ierr = KSPSolve(ksp,b,b);CHKERRQ(ierr);
		ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
		ierr = MatDestroy(&M_temp);CHKERRQ(ierr);
			
		ierr = MatMult(V,b,yt);CHKERRQ(ierr);
		ierr = VecAXPY(y1,-1.,yt);CHKERRQ(ierr);

		ierr = MatDenseGetSubMatrix(TSR->Y_hat_up,PETSC_DECIDE,PETSC_DECIDE,ik*RSVDt->RSVD.Nw_eff,(ik+1)*RSVDt->RSVD.Nw_eff,&Y_all_k);CHKERRQ(ierr);
		ierr = MatDenseGetColumnVecWrite(Y_all_k,iw,&y_iw);CHKERRQ(ierr);
		ierr = VecCopy(y1,y_iw);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(Y_all_k,iw,&y_iw);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(TSR->Y_hat_up,&Y_all_k);CHKERRQ(ierr);

	}

	ierr = VecDestroy(&y1);CHKERRQ(ierr);
	ierr = VecDestroy(&y2);CHKERRQ(ierr);
	ierr = VecDestroy(&b);CHKERRQ(ierr);
	ierr = VecDestroy(&yt);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
	
}

