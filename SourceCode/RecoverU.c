
#include <petscmat.h>
#include "Variables.h"

PetscErrorCode RecoverU(Mat U_til, RSVDt_vars *RSVDt, RSVD_matrices *RSVD_mat, Resolvent_matrices *Res_mat)
{
	/*
		Recovers \hat{U}  for all Nw matrices of size N \times k 
	*/

	PetscErrorCode        ierr;
	PetscInt              iw;
	Mat                   Q_temp, U_temp1, U_temp2;

	PetscFunctionBeginUser;

	ierr = MatCreate(PETSC_COMM_WORLD,&Res_mat->U_hat);CHKERRQ(ierr);
	ierr = MatSetType(Res_mat->U_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Res_mat->U_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,RSVDt->RSVD.k*RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(Res_mat->U_hat);CHKERRQ(ierr);

	for (iw=0; iw<RSVDt->RSVD.Nw_eff; iw++) {
		ierr = MatDenseGetSubMatrix(RSVD_mat->Q_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&Q_temp);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(U_til,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&U_temp1);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(Res_mat->U_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&U_temp2);CHKERRQ(ierr);
		ierr = MatMatMult(Q_temp,U_temp1,MAT_REUSE_MATRIX,PETSC_DEFAULT,&U_temp2);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Res_mat->U_hat,&U_temp2);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(U_til,&U_temp1);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(RSVD_mat->Q_hat,&Q_temp);CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(Res_mat->U_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Res_mat->U_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = MatDestroy(&RSVD_mat->Q_hat);CHKERRQ(ierr);
	ierr = MatDestroy(&U_til);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}



