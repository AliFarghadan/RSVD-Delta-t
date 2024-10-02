
#include <petscmat.h>
#include "Variables.h"

PetscErrorCode CreateDFTiDFTMats(RSVDt_vars *RSVDt, DFT_matrices *DFT_mat, LNS_vars *LNS_mat)
{
	/*
		Creates the DFT and inverse DFT matrices
	*/  

	PetscErrorCode        ierr;
	PetscInt              iw, it, Nw, Nt, start, end;
	PetscScalar           alpha, v;
	Vec                   V, J;

	PetscFunctionBeginUser;

	ierr = MatCreate(PETSC_COMM_WORLD,&DFT_mat->dft);CHKERRQ(ierr);
	ierr = MatSetType(DFT_mat->dft,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(DFT_mat->dft,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.Nw,RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(DFT_mat->dft);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&DFT_mat->idft);CHKERRQ(ierr);
	ierr = MatSetType(DFT_mat->idft,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(DFT_mat->idft,PETSC_DECIDE,PETSC_DECIDE,2*RSVDt->TS.Ns,RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(DFT_mat->idft);CHKERRQ(ierr);
	ierr = MatGetSize(DFT_mat->idft,&Nt,&Nw);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &J);CHKERRQ(ierr);
	ierr = VecSetSizes(J, PETSC_DECIDE, Nt);CHKERRQ(ierr);
	ierr = VecSetUp(J);CHKERRQ(ierr);

	VecGetOwnershipRange(J,&start,&end);
	for (it=start; it<end; it++) {
		v = it;
		VecSetValues(J,1,&it,&v,INSERT_VALUES);
	}

	ierr = VecAssemblyBegin(J);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(J);CHKERRQ(ierr);

	alpha = -2*PETSC_PI*PETSC_i/Nt;

	for (iw=0; iw<Nw; iw++) {
		ierr = MatDenseGetColumnVecWrite(DFT_mat->idft,iw,&V);CHKERRQ(ierr);
		if (RSVDt->TS.RealOperator) {
			v    = alpha*iw;
		} else if (PetscFmodReal(Nw,2) == 0)  {
			v    = (iw<=Nw/2-1) ? alpha*iw : alpha*((iw-Nw)+Nt);
		} else {
			v    = (iw<=(Nw-1)/2) ? alpha*iw : alpha*((iw-Nw)+Nt);
		}
		ierr = VecSet(V,v);CHKERRQ(ierr);
		ierr = VecPointwiseMult(V,V,J);CHKERRQ(ierr);
		ierr = VecExp(V);CHKERRQ(ierr);
		ierr = VecConjugate(V);CHKERRQ(ierr);
		ierr = VecScale(V, 2./Nt);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(DFT_mat->idft,iw,&V);CHKERRQ(ierr);
	}

	ierr = VecDestroy(&J);CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &J);CHKERRQ(ierr);
	ierr = RSVDt->TS.RealOperator ? VecSetSizes(J,PETSC_DECIDE,Nw*2) : \
						VecSetSizes(J,PETSC_DECIDE,Nw);CHKERRQ(ierr);
	ierr = VecSetUp(J);CHKERRQ(ierr);

	VecGetOwnershipRange(J,&start,&end);
	for (it=start; it<end; it++) {
		v = it;
		VecSetValues(J,1,&it,&v,INSERT_VALUES);
	}

	ierr = VecAssemblyBegin(J);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(J);CHKERRQ(ierr);

	alpha = RSVDt->TS.RealOperator ? -2*PETSC_PI*PETSC_i/(2*Nw) : -2*PETSC_PI*PETSC_i/Nw;
	
	for (iw=0; iw<Nw; iw++) {
		ierr = MatDenseGetColumnVecWrite(DFT_mat->dft,iw,&V);CHKERRQ(ierr);
		v    = alpha*iw;
		ierr = VecSet(V,v);CHKERRQ(ierr);
		ierr = VecPointwiseMult(V,V,J);CHKERRQ(ierr);
		ierr = VecExp(V);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(DFT_mat->dft,iw,&V);CHKERRQ(ierr);
	}

	ierr = VecDestroy(&J);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(DFT_mat->dft,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(DFT_mat->dft,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(DFT_mat->idft,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(DFT_mat->idft,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatTranspose(DFT_mat->idft,MAT_INPLACE_MATRIX,&DFT_mat->idft);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

