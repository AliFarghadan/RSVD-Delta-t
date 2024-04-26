
#include <SetupTimeFreqGrid.h>

PetscErrorCode SetupTimeFreqGrid(RSVDt_vars *RSVDt)
{

	PetscErrorCode        ierr;
	PetscReal             T_ss,dt_loc,dt_w_comm,Nw_Nw_out;
	PetscInt              Nt_loc,Nt_transient_loc;

	PetscFunctionBeginUser;

	// Nw_Nw_out                = RSVDt->RSVD.Nw * RSVDt->RSVD.Nw_out;
	Nw_Nw_out                = RSVDt->RSVD.Nw;
	T_ss                     = 2*PETSC_PI/RSVDt->RSVD.w_min;
	dt_w_comm                = T_ss/Nw_Nw_out;
	// dt_w                     = T_ss/RSVDt->RSVD.Nw;
	// dt_w_out                 = T_ss/RSVDt->RSVD.Nw_out;
	dt_loc                   = dt_w_comm/PetscCeilReal(dt_w_comm/RSVDt->TS.dt_max);
	RSVDt->TS.dt             = dt_loc;
	Nt_loc                   = round(T_ss/dt_loc);
	RSVDt->TS.Nt             = Nt_loc;
	Nt_transient_loc         = round(RSVDt->TS.t_transient/dt_loc);
	RSVDt->TS.Nt_transient   = Nt_transient_loc;
	RSVDt->TS.time_ratio     = Nt_loc/RSVDt->RSVD.Nw;
	RSVDt->TS.time_ratio_out = Nt_loc/RSVDt->RSVD.Nw_out;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"@ w_min = %g\n@ t_transient = %g\n@ dt = %g\n@ Nt = %d \
							\n@ Nt_transient = %d\n@ time_ratio = %d\n@ time_ratio_out = %d\n",\
							RSVDt->RSVD.w_min,RSVDt->TS.t_transient, \
							RSVDt->TS.dt,(int)RSVDt->TS.Nt,(int)RSVDt->TS.Nt_transient,\
							(int)RSVDt->TS.time_ratio,(int)RSVDt->TS.time_ratio_out);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

