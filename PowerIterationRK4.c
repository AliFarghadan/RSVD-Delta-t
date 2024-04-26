
#include <PowerIterationRK4.h>
#include <ForwardActionRK4.h>
#include <AdjointActionRK4.h>
#include <QR.h>

PetscErrorCode PowerIterationRK4(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, LNS_params *LNS_mat, \
			TS_matrices *TS_mat, DFT_matrices *DFT_mat, Sing_matrices *Sing_mat, Directories *dirs)
{

	PetscErrorCode        ierr;
	PetscInt              iq;

	PetscFunctionBeginUser;	

	for (iq=0; iq<RSVDt->RSVD.q; iq++) {

		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n ******* Inside power iteration, %d/%d *******",(int)iq+1,(int)RSVDt.RSVD.q);CHKERRQ(ierr);

		if (RSVDt->TS.flg_real_A) {
			ierr = QR(&RSVD_mat, RSVDt->RSVD.Nw/2);CHKERRQ(ierr);
		} else {
			ierr = QR(&RSVD_mat, RSVDt->RSVD.Nw);CHKERRQ(ierr);
		}

		ierr = ForwardActionRK4(RSVD_mat, RSVDt, LNS_mat, TS_mat, DFT_mat, Sing_mat, dirs);CHKERRQ(ierr);

		if (RSVDt->TS.flg_real_A) {
			ierr = QR(&RSVD_mat, RSVDt->RSVD.Nw/2);CHKERRQ(ierr);
		} else {
			ierr = QR(&RSVD_mat, RSVDt->RSVD.Nw);CHKERRQ(ierr);
		}
		
		ierr = AdjointActionRK4(RSVD_mat, RSVDt, LNS_mat, TS_mat, DFT_mat, Sing_mat, dirs);CHKERRQ(ierr);
		
	}

	PetscFunctionReturn(0);
	
}

