
#include <CreateDFTiDFTMats.h>

PetscErrorCode CreateDFTiDFTMats(RSVDt_vars *RSVDt, DFT_matrices *DFT_mat, LNS_params *LNS_mat, PetscBool flg_real_A)
{

	PetscErrorCode        ierr;
	PetscInt              iw, it, Nw, Nw_out, Nt, ib, Np, Nb, start, end, st_ind;
	PetscScalar           alpha, v;
	Vec                   V_temp, J;

	PetscFunctionBeginUser;

	/*
		Create the dft, inv_dft, and inv_dft_LNS matrices
	*/	

	ierr = MatCreate(PETSC_COMM_WORLD,&DFT_mat->dft_dir);CHKERRQ(ierr);
	ierr = MatSetType(DFT_mat->dft_dir,MATDENSE);CHKERRQ(ierr);
	if (RSVDt->TS.flg_real_A) {
		ierr = MatSetSizes(DFT_mat->dft_dir,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.Nw_out,RSVDt->RSVD.Nw_out/2);CHKERRQ(ierr);
	} else {
		ierr = MatSetSizes(DFT_mat->dft_dir,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.Nw_out,RSVDt->RSVD.Nw_out);CHKERRQ(ierr);
	}
	ierr = MatSetUp(DFT_mat->dft_dir);CHKERRQ(ierr); 

	ierr = MatCreate(PETSC_COMM_WORLD,&DFT_mat.inv_dft_dir);CHKERRQ(ierr);
	ierr = MatSetType(DFT_mat->inv_dft_dir,MATDENSE);CHKERRQ(ierr);
	if (RSVDt->TS.flg_real_A) {
		ierr = MatSetSizes(DFT_mat->inv_dft_dir,PETSC_DECIDE,PETSC_DECIDE,2*RSVDt->TS.Nt,RSVDt->RSVD.Nw/2);CHKERRQ(ierr);
	} else {
		ierr = MatSetSizes(DFT_mat->inv_dft_dir,PETSC_DECIDE,PETSC_DECIDE,2*RSVDt->TS.Nt,RSVDt->RSVD.Nw);CHKERRQ(ierr);
	}
	ierr = MatSetUp(DFT_mat->inv_dft_dir);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&DFT_mat.dft_adj);CHKERRQ(ierr);
	ierr = MatSetType(DFT_mat->dft_adj,MATDENSE);CHKERRQ(ierr);
	if (RSVDt->TS.flg_real_A) {
		ierr = MatSetSizes(DFT_mat->dft_adj,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.Nw,RSVDt->RSVD.Nw/2);CHKERRQ(ierr);
	} else {
		ierr = MatSetSizes(DFT_mat->dft_adj,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.Nw,RSVDt->RSVD.Nw);CHKERRQ(ierr);
	}
	ierr = MatSetUp(DFT_mat->dft_adj);CHKERRQ(ierr); 

	ierr = MatCreate(PETSC_COMM_WORLD,&DFT_mat->inv_dft_adj);CHKERRQ(ierr);
	ierr = MatSetType(DFT_mat->inv_dft_adj,MATDENSE);CHKERRQ(ierr);
	if (RSVDt->TS.flg_real_A) {
		ierr = MatSetSizes(DFT_mat->inv_dft_adj,PETSC_DECIDE,PETSC_DECIDE,2*RSVDt->TS.Nt,RSVDt->RSVD.Nw_out/2);CHKERRQ(ierr);
	} else {
		ierr = MatSetSizes(DFT_mat->inv_dft_adj,PETSC_DECIDE,PETSC_DECIDE,2*RSVDt->TS.Nt,RSVDt->RSVD.Nw_out);CHKERRQ(ierr);
	}
	ierr = MatSetUp(DFT_mat->inv_dft_adj);CHKERRQ(ierr);	

	ierr = MatCreate(PETSC_COMM_WORLD,&LNS_mat->LNS_DFT.inv_dft_LNS);CHKERRQ(ierr);
	ierr = MatSetType(LNS_mat->LNS_DFT.inv_dft_LNS,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(LNS_mat->LNS_DFT.inv_dft_LNS,PETSC_DECIDE,PETSC_DECIDE,2*RSVDt->TS.Nt,LNS_mat->RSVDt.LNS.Nb);CHKERRQ(ierr);	
	ierr = MatSetUp(LNS_mat->LNS_DFT.inv_dft_LNS);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&LNS_mat->LNS_DFT.dft_LNS);CHKERRQ(ierr);
	ierr = MatSetType(LNS_mat->LNS_DFT.dft_LNS,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(LNS_mat->LNS_DFT.dft_LNS,PETSC_DECIDE,PETSC_DECIDE,LNS_mat->RSVDt.LNS.Np,LNS_mat->RSVDt.LNS.Nb);CHKERRQ(ierr);
	ierr = MatSetUp(LNS_mat->LNS_DFT.dft_LNS);CHKERRQ(ierr);	

	ierr = MatGetSize(DFT_mat->inv_dft_dir,&Nt,&Nw);CHKERRQ(ierr);
	ierr = MatGetSize(DFT_mat->inv_dft_adj,NULL,&Nw_out);CHKERRQ(ierr);
	ierr = MatGetSize(LNS_mat->LNS_DFT.dft_LNS,&Np,&Nb);CHKERRQ(ierr);

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
		ierr = MatDenseGetColumnVecWrite(DFT_mat->inv_dft_dir,iw,&V_temp);CHKERRQ(ierr);
		if (flg_real_A) {
			v    = alpha*iw;
		} else if (PetscFmodReal(Nw,2) == 0)  {
			v    = (iw<=Nw/2-1) ? alpha*iw : alpha*((iw-Nw)+Nt);
		} else {
			v    = (iw<=(Nw-1)/2) ? alpha*iw : alpha*((iw-Nw)+Nt);
		}
		ierr = VecSet(V_temp,v);CHKERRQ(ierr);
		ierr = VecPointwiseMult(V_temp,V_temp,J);CHKERRQ(ierr);
		ierr = VecExp(V_temp);CHKERRQ(ierr);
		ierr = VecConjugate(V_temp);CHKERRQ(ierr);
		ierr = VecScale(V_temp, 2./Nt);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(DFT_mat->inv_dft_dir,iw,&V_temp);CHKERRQ(ierr);
	}

	for (iw=0; iw<Nw_out; iw++) {
		ierr = MatDenseGetColumnVecWrite(DFT_mat->inv_dft_adj,iw,&V_temp);CHKERRQ(ierr);
		if (flg_real_A) {
			v    = alpha*iw;
		} else if (PetscFmodReal(Nw_out,2) == 0)  {
			v    = (iw<=Nw_out/2-1) ? alpha*iw : alpha*((iw-Nw_out)+Nt);
		} else {
			v    = (iw<=(Nw_out-1)/2) ? alpha*iw : alpha*((iw-Nw_out)+Nt);
		}
		ierr = VecSet(V_temp,v);CHKERRQ(ierr);
		ierr = VecPointwiseMult(V_temp,V_temp,J);CHKERRQ(ierr);
		ierr = VecExp(V_temp);CHKERRQ(ierr);
		ierr = VecConjugate(V_temp);CHKERRQ(ierr);
		ierr = VecScale(V_temp, 2./Nt);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(DFT_mat->inv_dft_adj,iw,&V_temp);CHKERRQ(ierr);
	}	

	for (ib=0; ib<Nb; ib++) {
	ierr = MatDenseGetColumnVecWrite(LNS_mat->LNS_DFT.inv_dft_LNS,ib,&V_temp);CHKERRQ(ierr);
		if (LNS_mat->RSVDt.TS.flg_real_A_temp) {
			v    = alpha*ib;
		} else if (PetscFmodReal(Nb,2) == 0)  {
			v    = (ib<=Nb/2-1) ? alpha*ib : alpha*((ib-Nb)+Nt);
		} else {
			v    = (ib<=(Nb-1)/2) ? alpha*ib : alpha*((ib-Nb)+Nt);
		}	
		ierr = VecSet(V_temp,v);CHKERRQ(ierr);
		ierr = VecPointwiseMult(V_temp,V_temp,J);CHKERRQ(ierr);
		ierr = VecExp(V_temp);CHKERRQ(ierr);	
		ierr = LNS_mat->RSVDt.TS.flg_real_A_temp && ib>0 ? VecScale(V_temp, 2.) : \
							VecScale(V_temp, 1.);CHKERRQ(ierr);
		ierr = VecConjugate(V_temp);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(LNS_mat->LNS_DFT.inv_dft_LNS,ib,&V_temp);CHKERRQ(ierr);
	}

	ierr = VecDestroy(&J);CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &J);CHKERRQ(ierr);
	ierr = flg_real_A ? VecSetSizes(J,PETSC_DECIDE,Nw_out*2) : \
						VecSetSizes(J,PETSC_DECIDE,Nw_out);CHKERRQ(ierr);
	ierr = VecSetUp(J);CHKERRQ(ierr);

	VecGetOwnershipRange(J,&start,&end);
	for (it=start; it<end; it++) {
		v = it;
		VecSetValues(J,1,&it,&v,INSERT_VALUES);
	}

	ierr = VecAssemblyBegin(J);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(J);CHKERRQ(ierr);

	alpha = flg_real_A ? -2*PETSC_PI*PETSC_i/(2*Nw_out) : -2*PETSC_PI*PETSC_i/Nw_out;
	
	for (iw=0; iw<Nw_out; iw++) {
		ierr = MatDenseGetColumnVecWrite(DFT_mat->dft_dir,iw,&V_temp);CHKERRQ(ierr);
		v    = alpha*iw;
		ierr = VecSet(V_temp,v);CHKERRQ(ierr);
		ierr = VecPointwiseMult(V_temp,V_temp,J);CHKERRQ(ierr);
		ierr = VecExp(V_temp);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(DFT_mat->dft_dir,iw,&V_temp);CHKERRQ(ierr);
	}

	ierr = VecDestroy(&J);CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &J);CHKERRQ(ierr);
	ierr = flg_real_A ? VecSetSizes(J,PETSC_DECIDE,Nw*2) : \
						VecSetSizes(J,PETSC_DECIDE,Nw);CHKERRQ(ierr);
	ierr = VecSetUp(J);CHKERRQ(ierr);

	VecGetOwnershipRange(J,&start,&end);
	for (it=start; it<end; it++) {
		v = it;
		VecSetValues(J,1,&it,&v,INSERT_VALUES);
	}

	ierr = VecAssemblyBegin(J);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(J);CHKERRQ(ierr);

	alpha = flg_real_A ? -2*PETSC_PI*PETSC_i/(2*Nw) : -2*PETSC_PI*PETSC_i/Nw;
	
	for (iw=0; iw<Nw; iw++) {
		ierr = MatDenseGetColumnVecWrite(DFT_mat->dft_adj,iw,&V_temp);CHKERRQ(ierr);
		v    = alpha*iw;
		ierr = VecSet(V_temp,v);CHKERRQ(ierr);
		ierr = VecPointwiseMult(V_temp,V_temp,J);CHKERRQ(ierr);
		ierr = VecExp(V_temp);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(DFT_mat->dft_adj,iw,&V_temp);CHKERRQ(ierr);
	}

	ierr = VecDestroy(&J);CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &J);CHKERRQ(ierr);
	ierr =  flg_real_A ? VecSetSizes(J,PETSC_DECIDE,Np*2) : \
						VecSetSizes(J,PETSC_DECIDE,Np);CHKERRQ(ierr);
	ierr = VecSetUp(J);CHKERRQ(ierr);

	VecGetOwnershipRange(J,&start,&end);
	for (it=start; it<end; it++) {
		v = it;
		VecSetValues(J,1,&it,&v,INSERT_VALUES);
	}

	ierr = VecAssemblyBegin(J);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(J);CHKERRQ(ierr);

	alpha = -2*PETSC_PI*PETSC_i/Np;
	
	for (ib=0; ib<Nb; ib++) {
		ierr = MatDenseGetColumnVecWrite(LNS_mat->LNS_DFT.dft_LNS,ib,&V_temp);CHKERRQ(ierr);
		if (LNS_mat->RSVDt.TS.flg_real_A_temp) {
			v    = alpha*ib;
		} else if (PetscFmodReal(Nb,2) == 0)  {
			v    = (ib<=Nb/2-1) ? alpha*ib : alpha*((ib-Nb)+Np);
		} else {
			v    = (ib<=(Nb-1)/2) ? alpha*ib : alpha*((ib-Nb)+Np);
		}
		ierr = VecSet(V_temp,v);CHKERRQ(ierr);
		ierr = VecPointwiseMult(V_temp,V_temp,J);CHKERRQ(ierr);
		ierr = VecExp(V_temp);CHKERRQ(ierr);
		ierr = VecScale(V_temp, 1./Np);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(LNS_mat->LNS_DFT.dft_LNS,ib,&V_temp);CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(DFT_mat->dft_dir,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(DFT_mat->dft_dir,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(DFT_mat->dft_adj,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(DFT_mat->dft_adj,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);	
	ierr = MatAssemblyBegin(DFT_mat->inv_dft_dir,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(DFT_mat->inv_dft_dir,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(DFT_mat->inv_dft_adj,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(DFT_mat->inv_dft_adj,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);	
	ierr = MatAssemblyBegin(LNS_mat->LNS_DFT.inv_dft_LNS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(LNS_mat->LNS_DFT.inv_dft_LNS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(LNS_mat->LNS_DFT.dft_LNS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(LNS_mat->LNS_DFT.dft_LNS,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);	
	ierr = MatTranspose(DFT_mat->inv_dft_dir,MAT_INPLACE_MATRIX,&DFT_mat->inv_dft_dir);CHKERRQ(ierr);
	ierr = MatTranspose(DFT_mat->inv_dft_adj,MAT_INPLACE_MATRIX,&DFT_mat->inv_dft_adj);CHKERRQ(ierr);
	ierr = MatTranspose(LNS_mat->LNS_DFT.inv_dft_LNS,MAT_INPLACE_MATRIX,&LNS_mat->LNS_DFT.inv_dft_LNS);CHKERRQ(ierr);

	ierr = VecDestroy(&J);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

