
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode TSTransRK4(LNS_vars *LNS, TS_matrices *TS, RSVDt_vars *RSVDt, Vec y)
{
	/*
		Performs one RK4 iteration assuming zero forcing 
	*/
	
	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	if (RSVDt->TS.DirAdj) {
		ierr = MatMult(LNS->A,y,TS->k1);CHKERRQ(ierr);
		ierr = VecWAXPY(TS->y_temp,RSVDt->TS.dt/2,TS->k1,y);CHKERRQ(ierr);
		ierr = MatMult(LNS->A,TS->y_temp,TS->k2);CHKERRQ(ierr);
		ierr = VecWAXPY(TS->y_temp,RSVDt->TS.dt/2,TS->k2,y);CHKERRQ(ierr);
		ierr = MatMult(LNS->A,TS->y_temp,TS->k3);CHKERRQ(ierr);
		ierr = VecWAXPY(TS->y_temp,RSVDt->TS.dt,TS->k3,y);CHKERRQ(ierr);
		ierr = MatMult(LNS->A,TS->y_temp,TS->k4);CHKERRQ(ierr);
	} else{
		ierr = MatMultHermitianTranspose(LNS->A,y,TS->k1);CHKERRQ(ierr);
		ierr = VecWAXPY(TS->y_temp,RSVDt->TS.dt/2,TS->k1,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTranspose(LNS->A,TS->y_temp,TS->k2);CHKERRQ(ierr);
		ierr = VecWAXPY(TS->y_temp,RSVDt->TS.dt/2,TS->k2,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTranspose(LNS->A,TS->y_temp,TS->k3);CHKERRQ(ierr);
		ierr = VecWAXPY(TS->y_temp,RSVDt->TS.dt,TS->k3,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTranspose(LNS->A,TS->y_temp,TS->k4);CHKERRQ(ierr);
	}

	ierr = VecAXPY(y,RSVDt->TS.dt/6,TS->k1);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6*2,TS->k2);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6*2,TS->k3);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6,TS->k4);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
	
}

