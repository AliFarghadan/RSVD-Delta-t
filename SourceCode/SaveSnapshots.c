
#include <petscmat.h>
#include "Variables.h"

PetscErrorCode SaveSnapshots(Vec y, PetscInt i, RSVDt_vars *RSVDt, Mat Y_all)
{
	/*
		Stores the response snapshots at equal time intervals
	*/

	PetscErrorCode        ierr;
	PetscInt              pos, Nstore;
	Vec                   Y_temp;

	PetscFunctionBeginUser;

	Nstore = RSVDt->TS.TransientRemoval ? RSVDt->RSVD.Nw+1 : RSVDt->RSVD.Nw;

	if (i > RSVDt->TS.Nt) {
		pos = i - RSVDt->TS.Nt;
		if (PetscFmodReal(pos,RSVDt->TS.ResRatio) == 0) {
			pos  = pos/RSVDt->TS.ResRatio - 1;
			pos  = PetscFmodReal(pos, Nstore);
			pos  = RSVDt->TS.DirAdj ? pos : Nstore-1-pos;
			pos  = PetscFmodReal(pos, Nstore);
			ierr = MatDenseGetColumnVecWrite(Y_all,pos,&Y_temp);CHKERRQ(ierr);
			ierr = VecCopy(y,Y_temp);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(Y_all,pos,&Y_temp);CHKERRQ(ierr);
		}
	}

	PetscFunctionReturn(0);
	
}



