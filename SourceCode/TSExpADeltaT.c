
#include <petscmat.h>
#include <Variables.h>
#include <CreateTSMats.h>
#include <DestroyTSMats.h>
#include <TSTransRK4.h>

PetscErrorCode TSExpADeltaT(LNS_vars *LNS, RSVDt_vars *RSVDt, Mat U)
{
	/*
		Performs RK4 time stepping of e^(A \Delta t)U_i, for i = 1, 2, ..., k, 
		where U_i is the ith column of U, the trial basis
	*/  

	PetscErrorCode       ierr;
	PetscInt             it, i, Nstore;
	Vec                  y1;
	TS_matrices          TS;

	PetscFunctionBeginUser;

	ierr = MatGetSize(U,NULL,&Nstore);CHKERRQ(ierr);
	ierr = CreateTSMats(&TS,RSVDt->RSVD.N);CHKERRQ(ierr);

	for (i=0; i<Nstore; i++) {

		ierr = MatDenseGetColumnVecWrite(U,i,&y1);CHKERRQ(ierr);

		for (it=0; it<RSVDt->TS.ResRatio; it++) {
			ierr = TSTransRK4(LNS, &TS, RSVDt, y1);CHKERRQ(ierr);
		}

		ierr = MatDenseRestoreColumnVecWrite(U,i,&y1);CHKERRQ(ierr);
	}

	ierr = DestroyTSMats(&TS);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}



