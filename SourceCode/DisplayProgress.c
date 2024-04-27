
#include <petscmat.h>
#include "Variables.h"

PetscErrorCode DisplayProgress(PetscInt i, PetscInt rend, PetscInt ik, PetscInt *prg_cnt, PetscLogDouble t1, RSVDt_vars *RSVDt, PetscBool inLoop)
{
	/*
		Displays the progress percantage, elapsed time, and estimated remaining time
	*/

	PetscErrorCode        ierr;
	PetscReal             T_rem;
	PetscLogDouble        t2;
	PetscInt              hh, mm, ss, hh_rem, mm_rem, ss_rem;

	PetscFunctionBeginUser;

	if (ik == 0 && inLoop && (double)i/rend >= 0.01*(*prg_cnt)) {
		if (RSVDt->display == 2) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Progress percentage: %d%%\n",(int)*prg_cnt);CHKERRQ(ierr);
			ierr = PetscTime(&t2);CHKERRQ(ierr);
			hh   = (t2-t1)/3600;
			mm   = (t2-t1-3600*hh)/60;
			ss   = t2-t1-3600*hh-mm*60;
			if (*prg_cnt < 100) ierr = PetscPrintf(PETSC_COMM_WORLD,"Elapsed time = %02d:%02d:%02d\n", \
							(int)hh, (int)mm, (int)ss);CHKERRQ(ierr);
		}
		*prg_cnt += 10;
	}

	if (!inLoop) {
		if (RSVDt->display >= 1) {
			ierr   = PetscTime(&t2);CHKERRQ(ierr);
			T_rem  = (RSVDt->RSVD.k-1-ik)*(t2-t1);
			hh     = (t2-t1)/3600;
			mm     = (t2-t1-3600*hh)/60;
			ss     = t2-t1-3600*hh-mm*60;
			hh_rem = T_rem/3600;
			mm_rem = (T_rem-3600*hh_rem)/60;
			ss_rem = T_rem-3600*hh_rem-mm_rem*60;

			if (ik<RSVDt->RSVD.k-1) {
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Elapsed time = %02d:%02d:%02d (k = %d out of %d), Estimated time remaining = %02d:%02d:%02d\n", \
									(int)hh, (int)mm, (int)ss, (int)ik+1, (int)RSVDt->RSVD.k, (int)hh_rem, (int)mm_rem, (int)ss_rem);CHKERRQ(ierr);
			} else {
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Elapsed time = %02d:%02d:%02d (k = %d out of %d)\n", \
									(int)hh, (int)mm, (int)ss, (int)ik+1, (int)RSVDt->RSVD.k);CHKERRQ(ierr);    
			}
		}
	}

	PetscFunctionReturn(0);

}

