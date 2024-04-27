
#include <petscmat.h>
#include "Variables.h"

PetscErrorCode SetupTimeFreqGrid(RSVDt_vars *RSVDt)
{
	/*
		Initializes the time integration variables including time step and length of time integration
	*/

	PetscErrorCode        ierr;
	PetscReal             T_ss,dt_w;

	PetscFunctionBeginUser;

	RSVDt->RSVD.w       *= (RSVDt->RSVD.twopi_flg) ? 2*PETSC_PI : 1;
	T_ss                     = 2*PETSC_PI/RSVDt->RSVD.w;
	dt_w                     = T_ss/RSVDt->RSVD.Nw;
	RSVDt->TS.dt             = dt_w/PetscCeilReal(dt_w/RSVDt->TS.dt);
	RSVDt->TS.Ns             = round(T_ss/RSVDt->TS.dt);
	RSVDt->TS.Nt             = round(RSVDt->TS.TransLen/RSVDt->TS.dt);
	RSVDt->TS.ResRatio       = RSVDt->TS.Ns/RSVDt->RSVD.Nw;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"@ w = %g\n@ TransLen = %g\n@ dt = %g\n@ Ns = %d \
							\n@ Nt = %d\n@ ResRatio = %d\n",\
							RSVDt->RSVD.w,RSVDt->TS.TransLen, \
							RSVDt->TS.dt,(int)RSVDt->TS.Ns,(int)RSVDt->TS.Nt,\
							(int)RSVDt->TS.ResRatio);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}




