
#include <SaveSnapshots.h>

PetscErrorCode SaveSnapshots(Vec y, PetscInt i, RSVDt_vars *RSVDt, Mat Y_all)
{

	PetscErrorCode        ierr;
	PetscInt              pos,gap,Nw_eff;
	Vec                   Y_temp;

	PetscFunctionBeginUser;

	gap    = (RSVDt->TS.flg_dir_adj) ? RSVDt->TS.time_ratio_out : RSVDt->TS.time_ratio;
	Nw_eff = (RSVDt->TS.flg_dir_adj) ? RSVDt->RSVD.Nw_out : RSVDt->RSVD.Nw;

	if (i > RSVDt->TS.Nt_transient) {
		pos = i;
		if (PetscFmodReal(pos,gap) == 0) {
			pos  = pos/gap;
			pos  = PetscFmodReal(pos, Nw_eff);
			pos  = (RSVDt->TS.flg_dir_adj) ? pos : Nw_eff-pos;
			pos  = PetscFmodReal(pos, Nw_eff);
			pos += (RSVDt->TS.flg_eff_trans) ? 1 : 0;
			ierr = MatDenseGetColumnVecWrite(Y_all,pos,&Y_temp);CHKERRQ(ierr);
			ierr = VecCopy(y,Y_temp);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(Y_all,pos,&Y_temp);CHKERRQ(ierr);
		}
	}

	PetscFunctionReturn(0);
	
}
