
#include <StoreQ.h>

PetscErrorCode StoreQ(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt)
{
	
	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->Q_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVD_mat->Q_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(RSVD_mat->Q_hat,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,RSVDt->RSVD.Nw_out*RSVDt->RSVD.k);CHKERRQ(ierr); // Nw_eff to be added
	ierr = MatSetUp(RSVD_mat->Q_hat);CHKERRQ(ierr);	
	ierr = MatCopy(RSVD_mat->Y_hat,RSVD_mat->Q_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}
