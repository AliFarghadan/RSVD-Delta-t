
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode SetupTimeFreqGrid(RSVDt_vars *RSVDt, TransRun_vars *TR_vars)
{
	/*
		Initializes the time integration variables including time step and length of time integration
	*/

	PetscErrorCode        ierr=0;
	PetscReal             T_ss,dt_w;

	PetscFunctionBeginUser;

	T_ss                = 2*PETSC_PI/RSVDt->RSVD.w;
	dt_w                = T_ss/RSVDt->RSVD.Nw;
	RSVDt->TS.dt        = dt_w/PetscCeilReal(dt_w/RSVDt->TS.dt);
	RSVDt->TS.Ns        = round(T_ss/RSVDt->TS.dt);
	RSVDt->TS.Nt        = round(RSVDt->TS.TransientLength/RSVDt->TS.dt);
	RSVDt->TS.ResRatio  = RSVDt->TS.Ns/RSVDt->RSVD.Nw;

	if (!TR_vars->TransRun) {

		ierr = PetscPrintf(PETSC_COMM_WORLD,\
			"RSVD parameters:\n@ k        = %d\n@ q        = %d\n@ w        = %g\n\n",\
			(int)RSVDt->RSVD.k,(int)RSVDt->RSVD.q,RSVDt->RSVD.w);CHKERRQ(ierr);

		ierr = PetscPrintf(PETSC_COMM_WORLD,"Time-stepping parameters:\n"
			"@ TransientLength = %g\n@ dt              = %g\n@ Ns              = %d\n@ Nt              = %d\n",\
			RSVDt->TS.TransientLength, RSVDt->TS.dt,(int)RSVDt->TS.Ns,(int)RSVDt->TS.Nt);CHKERRQ(ierr);
		
		if (RSVDt->TS.TransientRemoval) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Transient removal strategy is requested\n\n");CHKERRQ(ierr);
		} else {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Transient removal strategy is NOT requested\n\n");CHKERRQ(ierr);
		}
	}

	PetscFunctionReturn(0);
	
}




