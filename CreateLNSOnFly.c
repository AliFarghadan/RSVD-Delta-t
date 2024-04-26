
#include <CreateLNSOnFly.h>

PetscErrorCode CreateLNSOnFly(LNS_params *LNS_mat, PetscInt jt_cyc, Mat A)
{

	PetscErrorCode        ierr;
	Vec                   V_read;
	const PetscScalar    *A_array;
	PetscLogDouble        t1, t2;

	PetscFunctionBeginUser;

	ierr = MatDenseGetColumnVecRead(LNS_mat->LNS_DFT.inv_dft_LNS,jt_cyc,&V_read);CHKERRQ(ierr);
	if (LNS_mat->RSVDt.TS.flg_real_A_temp) {
		ierr = MatMult(LNS_mat->A_hat,V_read,LNS_mat->A_nnz);CHKERRQ(ierr);	
		ierr = VecRealPart(LNS_mat->A_nnz);CHKERRQ(ierr);
	} else {
		ierr = MatMult(LNS_mat->A_hat,V_read,LNS_mat->A_nnz);CHKERRQ(ierr);
	}
	ierr = MatDenseRestoreColumnVecRead(LNS_mat->LNS_DFT.inv_dft_LNS,jt_cyc,&V_read);CHKERRQ(ierr);

	ierr = VecGetArrayRead(LNS_mat->A_nnz, &A_array);CHKERRQ(ierr);
	ierr = MatUpdateMPIAIJWithArray(A, A_array);CHKERRQ(ierr);
	ierr = VecRestoreArrayRead(LNS_mat->A_nnz, &A_array);CHKERRQ(ierr);	

	PetscFunctionReturn(0);
	
}

