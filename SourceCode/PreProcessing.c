
#include <petscmat.h>
#include "Variables.h"
#include "ReadUserInput.h"
#include "LNSStructure.h"
#include "SetupTimeFreqGrid.h"
#include "CreateDFTiDFTMats.h"
#include "CreateRandomMat.h"
#include "ReadWSMats.h"

PetscErrorCode PreProcessing(RSVDt_vars *RSVDt, WS_matrices *WS_mat, LNS_vars *LNS_mat, RSVD_matrices *RSVD_mat, \
					Resolvent_matrices *Res_mat, TransientRun_vars *TR_vars, DFT_matrices *DFT_mat, Directories *dirs)
{
	/*
		Read user inputs and create required matrices before running the algorithm
	*/

	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	/*
		Reads in user input parameters
	*/

	ierr = ReadUserInput(RSVDt, WS_mat, LNS_mat, Res_mat, TR_vars, dirs);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n********************************************\n"
			"*************** Problem info ***************\n********************************************\n");CHKERRQ(ierr);

	ierr = RSVDt->TS.flg_real_A ? PetscPrintf(PETSC_COMM_WORLD,"\n*** The LNS operator is real-valued! ***\n\n") : \
					PetscPrintf(PETSC_COMM_WORLD,"\n*** The LNS operator is complex-valued! ***\n\n");CHKERRQ(ierr);

	/*
		Reads in the LNS operator (+ discounting if desired)
	*/

	ierr = LNSStructure(LNS_mat, RSVDt, dirs);CHKERRQ(ierr);

	/*
		Initializes the time-stepping variables
	*/

	ierr = SetupTimeFreqGrid(RSVDt);CHKERRQ(ierr);

	/*
		Creates the discrete Fourier transform (DFT) and inverse DFT matrices
	*/

	ierr = CreateDFTiDFTMats(RSVDt, DFT_mat, LNS_mat);CHKERRQ(ierr);

	/*
		Generates a random input matrix (or read in if desired)
	*/

	ierr = CreateRandomMat(RSVD_mat, RSVDt, dirs);CHKERRQ(ierr);

	/*
		Read weight and spatial matrices (if applicable)
	*/

	ierr = ReadWSMats(WS_mat, dirs);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}


