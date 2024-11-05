
#include <petscmat.h>
#include <Variables.h>
#include <CreateTSMats.h>
#include <DestroyTSMats.h>
#include <TSTransRK4.h>

PetscErrorCode TSExpADeltaT(LNS_vars *LNS_mat, RSVDt_vars *RSVDt, Mat U)
{
	/*
		Performs RK4 time stepping of e^(A \Delta t)U_i, for i = 1, 2, ..., k, 
		where U_i is the ith column of U, the trial basis
	*/  

	PetscErrorCode       ierr;
	PetscInt             it, is, Nstore;
	Vec                  y1;
	TS_matrices          TS_mat;

	PetscFunctionBeginUser;

	ierr = MatGetSize(U,NULL,&Nstore);CHKERRQ(ierr);
	ierr = CreateTSMats(&TS_mat,RSVDt->RSVD.N);CHKERRQ(ierr);

	for (is=0; is<Nstore; is++) {

		ierr = MatDenseGetColumnVecWrite(U,is,&y1);CHKERRQ(ierr);

		for (it=0; it<RSVDt->TS.ResRatio; it++) {
			ierr = TSTransRK4(LNS_mat, &TS_mat, RSVDt, it, y1);CHKERRQ(ierr);
		}

		ierr = MatDenseRestoreColumnVecWrite(U,is,&y1);CHKERRQ(ierr);
	}

	ierr = DestroyTSMats(&TS_mat);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}



