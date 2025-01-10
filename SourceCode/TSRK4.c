
#include <petscmat.h>
#include <Variables.h>
#include <CreateForcingOnFly.h>

PetscErrorCode TSRK4(Mat F_hat, DFT_matrices *DFT, LNS_vars *LNS, \
					TS_matrices *TS, RSVDt_vars *RSVDt, PetscInt i, Vec y)
{
	/*
		Performs one RK4 iteration after creating forcing terms on fly
	*/
	
	PetscErrorCode        ierr;
	PetscInt              jt_cyc;

	PetscFunctionBeginUser;

	ierr = RSVDt->TS.DirAdj ? MatMultAdd(LNS->A,y,TS->F3,TS->k1) : \
				MatMultHermitianTransposeAdd(LNS->A,y,TS->F3,TS->k1);CHKERRQ(ierr);

	jt_cyc = RSVDt->TS.DirAdj ? PetscFmodReal(2*(i-1)+1,2*RSVDt->TS.Ns) : \
					PetscFmodReal(1e6*RSVDt->TS.Ns-1-(2*(i-1)), 2*RSVDt->TS.Ns);
	ierr = CreateForcingOnFly(F_hat,DFT,jt_cyc,TS->F2);CHKERRQ(ierr);

	jt_cyc = RSVDt->TS.DirAdj ? PetscFmodReal(2*(i-1)+2,2*RSVDt->TS.Ns) : \
					PetscFmodReal(1e6*RSVDt->TS.Ns-1-(2*(i-1)+1), 2*RSVDt->TS.Ns);
	ierr = CreateForcingOnFly(F_hat,DFT,jt_cyc,TS->F3);CHKERRQ(ierr);

	if (RSVDt->TS.DirAdj) {
		ierr = VecWAXPY(TS->y_temp,RSVDt->TS.dt/2,TS->k1,y);CHKERRQ(ierr);
		ierr = MatMultAdd(LNS->A,TS->y_temp,TS->F2,TS->k2);CHKERRQ(ierr);
		ierr = VecWAXPY(TS->y_temp,RSVDt->TS.dt/2,TS->k2,y);CHKERRQ(ierr);
		ierr = MatMultAdd(LNS->A,TS->y_temp,TS->F2,TS->k3);CHKERRQ(ierr);
		ierr = VecWAXPY(TS->y_temp,RSVDt->TS.dt,TS->k3,y);CHKERRQ(ierr);
		ierr = MatMultAdd(LNS->A,TS->y_temp,TS->F3,TS->k4);CHKERRQ(ierr);
	} else {
		ierr = VecWAXPY(TS->y_temp,RSVDt->TS.dt/2,TS->k1,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTransposeAdd(LNS->A,TS->y_temp,TS->F2,TS->k2);CHKERRQ(ierr);
		ierr = VecWAXPY(TS->y_temp,RSVDt->TS.dt/2,TS->k2,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTransposeAdd(LNS->A,TS->y_temp,TS->F2,TS->k3);CHKERRQ(ierr);
		ierr = VecWAXPY(TS->y_temp,RSVDt->TS.dt,TS->k3,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTransposeAdd(LNS->A,TS->y_temp,TS->F3,TS->k4);CHKERRQ(ierr);
	}
	ierr = VecAXPY(y,RSVDt->TS.dt/6,TS->k1);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6*2,TS->k2);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6*2,TS->k3);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6,TS->k4);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

