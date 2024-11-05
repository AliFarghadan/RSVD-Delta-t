
#include <petscmat.h>
#include <Variables.h>
#include <QRDecomposition.h>
#include <GalerkinProj.h>

PetscErrorCode GalerkinProjAllModes(Mat Y_all, LNS_vars *LNS_mat, RSVDt_vars *RSVDt, \
								TS_removal_matrices *TSR, DFT_matrices *DFT_mat)
{
	/*
		Performs Galerkin projection to remove the transient response separately for each test vector
	*/

	PetscErrorCode        ierr;
	Mat                   V,Y_all_k,M_tilde;
	PetscInt              Nstore,k,ik;

	PetscFunctionBeginUser;

	ierr = MatGetSize(TSR->M_tilde_all,&Nstore,&k);CHKERRQ(ierr);
	k   /= Nstore;

	ierr = MatCreate(PETSC_COMM_WORLD,&TSR->Y_hat_up);CHKERRQ(ierr);
	ierr = MatSetType(TSR->Y_hat_up,MATMPIDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(TSR->Y_hat_up,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,RSVDt->RSVD.Nw_eff*k);CHKERRQ(ierr);
	ierr = MatSetUp(TSR->Y_hat_up);CHKERRQ(ierr);

	for (ik=0; ik<k; ik++) {

		ierr = MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nstore,(ik+1)*Nstore,&Y_all_k);CHKERRQ(ierr);
		ierr = MatDuplicate(Y_all_k,MAT_COPY_VALUES,&V);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Y_all,&Y_all_k);CHKERRQ(ierr);
		ierr = QRDecomposition(V);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(TSR->M_tilde_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nstore,(ik+1)*Nstore,&M_tilde);CHKERRQ(ierr);
		ierr = GalerkinProj(Y_all,M_tilde,V,ik,RSVDt,TSR,DFT_mat);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(TSR->M_tilde_all,&M_tilde);CHKERRQ(ierr);

		ierr = MatDestroy(&V);CHKERRQ(ierr);
	}

	ierr = MatDestroy(&TSR->M_tilde_all);CHKERRQ(ierr);
	ierr = MatDestroy(&Y_all);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
	
}

