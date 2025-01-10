
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode SetupTimeFreqGrid(RSVDt_vars *RSVDt, TransRun_vars *TR)
{
	/*
		Initializes the time integration variables including time step, time integration length, etc.
	*/

	PetscErrorCode        ierr=0;
	PetscReal             T_ss,dt_w,Nw_tmp,w_min,w_max;

	PetscFunctionBeginUser;

	Nw_tmp              = RSVDt->RSVD.Nw;
	RSVDt->RSVD.Nw     += RSVDt->TS.RealOperator && PetscFmodReal(RSVDt->RSVD.Nw, 2) == 1 ? 1 : 0;
	RSVDt->RSVD.Nw_eff  = RSVDt->TS.RealOperator ? RSVDt->RSVD.Nw/2 : RSVDt->RSVD.Nw;
	T_ss                = 2*PETSC_PI/RSVDt->RSVD.w;
	dt_w                = T_ss/RSVDt->RSVD.Nw;
	RSVDt->TS.dt        = dt_w/PetscCeilReal(dt_w/RSVDt->TS.dt);
	RSVDt->TS.Ns        = round(T_ss/RSVDt->TS.dt);
	RSVDt->TS.Nt        = round(RSVDt->TS.TransientLength/RSVDt->TS.dt);
	RSVDt->TS.ResRatio  = RSVDt->TS.Ns/RSVDt->RSVD.Nw;
	w_min               = RSVDt->TS.RealOperator ? 0 : -RSVDt->RSVD.Nw/2 * RSVDt->RSVD.w;
	w_max               = PetscFmodReal(RSVDt->RSVD.Nw, 2) == 1 ? RSVDt->RSVD.Nw/2 * RSVDt->RSVD.w : \
								(RSVDt->RSVD.Nw/2 - 1) * RSVDt->RSVD.w;

	if (!TR->TransRun && RSVDt->Display) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,\
			"RSVD parameters:\n@ k     = %d\n@ q     = %d\n@ w     = %g\n@ Nw    = %d",\
			(int)RSVDt->RSVD.k,(int)RSVDt->RSVD.q,RSVDt->RSVD.w,(int)Nw_tmp);CHKERRQ(ierr);

		ierr = RSVDt->TS.RealOperator ? PetscPrintf(PETSC_COMM_WORLD," --> The LNS operator is real-valued --> ") : \
		PetscPrintf(PETSC_COMM_WORLD," --> The LNS operator is complex-valued --> ");CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"The effective number of frequencies, Nw = %d\n",(int)RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);

		ierr = PetscPrintf(PETSC_COMM_WORLD,"@ w_min = %g\n@ w_max = %g\n\n",w_min,w_max);CHKERRQ(ierr);

		ierr = PetscPrintf(PETSC_COMM_WORLD,"Time-stepping parameters:\n"
			"@ TransientLength   = %g\n@ SteadyStatePeriod = %g\n@ dt                = %g\n",\
			RSVDt->TS.TransientLength, T_ss, RSVDt->TS.dt);CHKERRQ(ierr);
		
		if (RSVDt->TS.TransientRemoval) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Transient removal strategy is requested\n\n");CHKERRQ(ierr);
		} else {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Transient removal strategy is NOT requested\n\n");CHKERRQ(ierr);
		}
	}

	PetscFunctionReturn(0);
	
}




