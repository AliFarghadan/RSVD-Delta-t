
#include <petscmat.h>
#include "Variables.h"
#include "ReadUserInput.h"
#include "LNSStructure.h"
#include "SetupTimeFreqGrid.h"
#include "CreateDFTiDFTMats.h"
#include "CreateRandomMat.h"
#include "CreateResultsDir.h"
#include "ReadWSMats.h"

PetscErrorCode PreProcessing(RSVDt_vars *RSVDt, WS_matrices *WS_mat, LNS_vars *LNS_mat, RSVD_matrices *RSVD_mat, \
					Resolvent_matrices *Res_mat, TransRun_vars *TR_vars, DFT_matrices *DFT_mat, Directories *dirs)
{
	/*
		Reads user inputs and creates required matrices before running the algorithm
	*/

	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	/*
		Reads in user input parameters
	*/

	ierr = ReadUserInput(RSVDt, WS_mat, LNS_mat, Res_mat, TR_vars, dirs);

	if (!TR_vars->TransRun) ierr = PetscPrintf(PETSC_COMM_WORLD,"********************************************\n"
			"*************** Problem info ***************\n********************************************\n\n");CHKERRQ(ierr);

	/*
		Initializes the time-stepping variables
	*/

	ierr = SetupTimeFreqGrid(RSVDt, TR_vars);CHKERRQ(ierr);

	/*
		Loads the LNS operator (+ discounting if desired)
	*/

	ierr = LNSStructure(LNS_mat, RSVDt, TR_vars, dirs);CHKERRQ(ierr);

	/*
		Creates the discrete Fourier transform (DFT) and inverse DFT (iDFT) matrices
	*/

	ierr = CreateDFTiDFTMats(RSVDt, DFT_mat, LNS_mat);CHKERRQ(ierr);
	
	/*
		Generates a random input matrix (or read in if desired)
	*/

	if (!TR_vars->TransRun) ierr = CreateRandomMat(RSVD_mat, RSVDt, dirs);CHKERRQ(ierr);

	/*
		Loads the weight and spatial matrices (if applicable)
	*/

	if (!TR_vars->TransRun) ierr = ReadWSMats(WS_mat, dirs);CHKERRQ(ierr);

	/*
		Creates a folder for the results
	*/

	ierr = !TR_vars->TransRun ? CreateResultsDir(dirs, "ResolventModes_") : \
				CreateResultsDir(dirs, "TransientSnapshots_");CHKERRQ(ierr); 	

	PetscFunctionReturn(0);

}


