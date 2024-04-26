
#include <TSTransRK4.h>

PetscErrorCode TSTransRK4(LNS_params *LNS_mat, \
					TS_matrices *TS_mat, RSVDt_vars *RSVDt, PetscInt i, Vec y)
{
	
	PetscErrorCode       ierr;
	PetscInt             jt_cyc;
	PetscLogDouble       t1, t2;

	PetscFunctionBeginUser;
	
	if (RSVDt->TS.flg_dir_adj) {
		ierr = MatMult(TS_mat->A3,y,TS_mat->k1);CHKERRQ(ierr);
	} else{
		ierr = MatMultHermitianTranspose(TS_mat->A3,y,TS_mat->k1);CHKERRQ(ierr);
	}

	jt_cyc = RSVDt->TS.flg_dir_adj ? PetscFmodReal(2*(i-1)+1,2*RSVDt->TS.Nt) : \
					PetscFmodReal(1e6*RSVDt->TS.Nt-1-(2*(i-1)+1), 2*RSVDt->TS.Nt);
	ierr = CreateStreamingLNS(LNS_mat,jt_cyc,TS_mat->A2);CHKERRQ(ierr);

	jt_cyc = RSVDt->TS.flg_dir_adj ? PetscFmodReal(2*(i-1)+2,2*RSVDt->TS.Nt) : \
					PetscFmodReal(1e6*RSVDt->TS.Nt-1-(2*(i-1)+2), 2*RSVDt->TS.Nt);
	ierr = CreateStreamingLNS(LNS_mat,jt_cyc,TS_mat->A3);CHKERRQ(ierr);

	if (RSVDt->TS.flg_dir_adj) {
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k1,y);CHKERRQ(ierr);
		ierr = MatMult(TS_mat->A2,TS_mat->y_temp,TS_mat->k2);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k2,y);CHKERRQ(ierr);
		ierr = MatMult(TS_mat->A2,TS_mat->y_temp,TS_mat->k3);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt,TS_mat->k3,y);CHKERRQ(ierr);
		ierr = MatMult(TS_mat->A3,TS_mat->y_temp,TS_mat->k4);CHKERRQ(ierr);
	} else{
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k1,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTranspose(TS_mat->A2,TS_mat->y_temp,TS_mat->k2);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt/2,TS_mat->k2,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTranspose(TS_mat->A2,TS_mat->y_temp,TS_mat->k3);CHKERRQ(ierr);
		ierr = VecWAXPY(TS_mat->y_temp,RSVDt->TS.dt,TS_mat->k3,y);CHKERRQ(ierr);
		ierr = MatMultHermitianTranspose(TS_mat->A3,TS_mat->y_temp,TS_mat->k4);CHKERRQ(ierr);	
	}

	ierr = VecAXPY(y,RSVDt->TS.dt/6,TS_mat->k1);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6*2,TS_mat->k2);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6*2,TS_mat->k3);CHKERRQ(ierr);
	ierr = VecAXPY(y,RSVDt->TS.dt/6,TS_mat->k4);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
	
}

