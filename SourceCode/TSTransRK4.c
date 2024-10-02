
#include <petscmat.h>
#include "Variables.h"

PetscErrorCode TSTransRK4(LNS_vars *LNS_mat, TS_matrices *TS_mat, RSVDt_vars *RSVDt, PetscInt i, Vec y)
{
	/*
		Performs one RK4 iteration assuming zero forcing 
	*/
	
	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	if (RSVDt->TS.DirAdj) {
		ierr = MatMult(LNS_mat->A,y,TS_mat->k1);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k1,y);CHKERRQ(ierr);
		ierr = MatMult(LNS_mat->A,TS_mat->y_temp,TS_mat->k2);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k2,y);CHKERRQ(ierr);
		ierr = MatMult(LNS_mat->A,TS_mat->y_temp,TS_mat->k3);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt,TS_mat->k3,y);CHKERRQ(ierr);
		ierr = MatMult(LNS_mat->A,TS_mat->y_temp,TS_mat->k4);CHKERRQ(ierr);
	} else{
		ierr = MatMultHermitianTranspose(LNS_mat->A,y,TS_mat->k1);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k1,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTranspose(LNS_mat->A,TS_mat->y_temp,TS_mat->k2);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k2,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTranspose(LNS_mat->A,TS_mat->y_temp,TS_mat->k3);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt,TS_mat->k3,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTranspose(LNS_mat->A,TS_mat->y_temp,TS_mat->k4);CHKERRQ(ierr);
	}

	ierr = VecAXPY(y,RSVDt->TS.dt/6,TS_mat->k1);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6*2,TS_mat->k2);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6*2,TS_mat->k3);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6,TS_mat->k4);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
	
}

