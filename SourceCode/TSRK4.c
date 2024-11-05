
#include <petscmat.h>
#include <Variables.h>
#include <CreateForcingOnFly.h>

PetscErrorCode TSRK4(Mat F_hat, DFT_matrices *DFT_mat, LNS_vars *LNS_mat, \
					TS_matrices *TS_mat, RSVDt_vars *RSVDt, PetscInt i, Vec y)
{
	/*
		Performs one RK4 iteration after creating forcing terms on fly
	*/
	
	PetscErrorCode        ierr;
	PetscInt              jt_cyc;

	PetscFunctionBeginUser;

	ierr = RSVDt->TS.DirAdj ? MatMultAdd(LNS_mat->A,y,TS_mat->F3,TS_mat->k1) : \
				MatMultHermitianTransposeAdd(LNS_mat->A,y,TS_mat->F3,TS_mat->k1);CHKERRQ(ierr);

	jt_cyc = RSVDt->TS.DirAdj ? PetscFmodReal(2*(i-1)+1,2*RSVDt->TS.Ns) : \
					PetscFmodReal(1e6*RSVDt->TS.Ns-1-(2*(i-1)), 2*RSVDt->TS.Ns);
	ierr = CreateForcingOnFly(F_hat,DFT_mat,jt_cyc,TS_mat->F2);CHKERRQ(ierr);

	jt_cyc = RSVDt->TS.DirAdj ? PetscFmodReal(2*(i-1)+2,2*RSVDt->TS.Ns) : \
					PetscFmodReal(1e6*RSVDt->TS.Ns-1-(2*(i-1)+1), 2*RSVDt->TS.Ns);
	ierr = CreateForcingOnFly(F_hat,DFT_mat,jt_cyc,TS_mat->F3);CHKERRQ(ierr);

	if (RSVDt->TS.DirAdj) {
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k1,y);CHKERRQ(ierr);
		ierr = MatMultAdd(LNS_mat->A,TS_mat->y_temp,TS_mat->F2,TS_mat->k2);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k2,y);CHKERRQ(ierr);
		ierr = MatMultAdd(LNS_mat->A,TS_mat->y_temp,TS_mat->F2,TS_mat->k3);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt,TS_mat->k3,y);CHKERRQ(ierr);
		ierr = MatMultAdd(LNS_mat->A,TS_mat->y_temp,TS_mat->F3,TS_mat->k4);CHKERRQ(ierr);
	} else {
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k1,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTransposeAdd(LNS_mat->A,TS_mat->y_temp,TS_mat->F2,TS_mat->k2);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k2,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTransposeAdd(LNS_mat->A,TS_mat->y_temp,TS_mat->F2,TS_mat->k3);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt,TS_mat->k3,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTransposeAdd(LNS_mat->A,TS_mat->y_temp,TS_mat->F3,TS_mat->k4);CHKERRQ(ierr);
	}
	ierr = VecAXPY(y,RSVDt->TS.dt/6,TS_mat->k1);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6*2,TS_mat->k2);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6*2,TS_mat->k3);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6,TS_mat->k4);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

