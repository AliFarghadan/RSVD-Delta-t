
#include <petscmat.h>
#include "Variables.h"
#include "QRDecomposition.h"
#include "TSExpADeltaT.h"

PetscErrorCode CreateMtilde(Mat Y_all, LNS_vars *LNS_mat, RSVDt_vars *RSVDt, TS_removal_matrices *TSR)
{
	/*
		Creates M_tilde_all = U_i'exp(A\Delta t)U_i for i = 1,2,...,k test vectors in a loop
		U_i is the trial basis (matrix) for the ith test vector built from Nw+1 snapshots in time: U_i = QR(Y_all_i)
	*/  

	PetscErrorCode        ierr;
	Mat                   Y_all_k,V,MV;
	PetscInt              Nstore,ik;

	PetscFunctionBeginUser;

	Nstore = RSVDt->RSVD.Nw + 1;

	ierr = MatCreate(PETSC_COMM_WORLD,&TSR->M_tilde_all);CHKERRQ(ierr);
	ierr = MatSetType(TSR->M_tilde_all,MATMPIDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(TSR->M_tilde_all,PETSC_DECIDE,PETSC_DECIDE,Nstore,Nstore*RSVDt->RSVD.k);CHKERRQ(ierr);
	ierr = MatSetUp(TSR->M_tilde_all);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&V);CHKERRQ(ierr);
	ierr = MatSetType(V,MATMPIDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(V,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,Nstore);CHKERRQ(ierr);
	ierr = MatSetUp(V);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(V,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(V,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&MV);CHKERRQ(ierr);
	ierr = MatSetType(MV,MATMPIDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(MV,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,Nstore);CHKERRQ(ierr);
	ierr = MatSetUp(MV);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(MV,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(MV,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	for (ik=0; ik<RSVDt->RSVD.k; ik++) {

		ierr = MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nstore,(ik+1)*Nstore,&Y_all_k);CHKERRQ(ierr);
		ierr = MatCopy(Y_all_k,V,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Y_all,&Y_all_k);CHKERRQ(ierr);
		ierr = QRDecomposition(V);CHKERRQ(ierr);
		ierr = MatCopy(V,MV,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		ierr = TSExpADeltaT(LNS_mat, RSVDt, MV);CHKERRQ(ierr);
		ierr = MatConjugate(V);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(TSR->M_tilde_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nstore,(ik+1)*Nstore,&Y_all_k);CHKERRQ(ierr);
		ierr = MatTransposeMatMult(V,MV,MAT_REUSE_MATRIX,PETSC_DEFAULT,&Y_all_k);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(TSR->M_tilde_all,&Y_all_k);CHKERRQ(ierr);

	}

	ierr = MatDestroy(&MV);CHKERRQ(ierr);
	ierr = MatDestroy(&V);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
	
}



