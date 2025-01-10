
#include <petscmat.h>
#include <Variables.h>
#include <QRAllFreqs.h>
#include <AdjointActionRK4.h>
#include <DirectActionRK4.h>

PetscErrorCode PowerIterationRK4(RSVD_matrices *RSVD, RSVDt_vars *RSVDt, LNS_vars *LNS, \
			DFT_matrices *DFT, Weight_matrices *Weight, Directories *dirs, TS_removal_matrices *TSR)
{
	/*
		Performs power iteration for q times 
	*/

	PetscErrorCode        ierr;
	PetscInt              iq;

	PetscFunctionBeginUser;

	for (iq=0; iq<RSVDt->RSVD.q; iq++) {

		if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n******** Inside power iteration, %d/%d *******\n",(int)iq+1,(int)RSVDt->RSVD.q);CHKERRQ(ierr);

		ierr = QRAllFreqs(RSVD, RSVDt);CHKERRQ(ierr);
		ierr = AdjointActionRK4(RSVD, RSVDt, LNS, DFT, Weight, dirs, TSR);CHKERRQ(ierr);
		ierr = QRAllFreqs(RSVD, RSVDt);CHKERRQ(ierr);
		ierr = DirectActionRK4(RSVD, RSVDt, LNS, DFT, Weight, dirs, TSR);CHKERRQ(ierr);
		
	}

	if (RSVDt->RSVD.q > 0 && RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n******** Power iteration DONE! *************\n");CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}



