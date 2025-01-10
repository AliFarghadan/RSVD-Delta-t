
#include <petscmat.h>
#include <Variables.h>
#include <ReadUserInput.h>
#include <LNSStructure.h>
#include <SetupTimeFreqGrid.h>
#include <CreateDFTiDFTMats.h>
#include <CreateRandomMat.h>
#include <CreateResultsDir.h>
#include <ReadWeightMats.h>
#include <SaveInputVarsCopy.h>

PetscErrorCode PreProcessing(RSVDt_vars *RSVDt, Weight_matrices *Weight_mat, LNS_vars *LNS, RSVD_matrices *RSVD, \
					TransRun_vars *TR, DFT_matrices *DFT, Directories *dirs)
{
	/*
		Reads user inputs and creates required matrices before running the algorithm
	*/

	PetscErrorCode        ierr;

	PetscFunctionBeginUser;
				
	/*
		Reads in user input parameters
	*/

	ierr = ReadUserInput(RSVDt, Weight_mat, LNS, TR, dirs);CHKERRQ(ierr);

	/*
		Creates a folder for the results
	*/

	ierr = !TR->TransRun ? CreateResultsDir(dirs, "ResolventModes_") : \
				CreateResultsDir(dirs, "TransientSnapshots_");CHKERRQ(ierr); 	

	/*
		Saves a copy of the input variables in the results folder
	*/

	ierr = SaveInputVarsCopy(dirs);CHKERRQ(ierr); 

	/*
		Loads the LNS operator (+ discounting if desired)
	*/

	ierr = LNSStructure(LNS, RSVDt, TR, dirs);CHKERRQ(ierr);

	/*
		Initializes the time-stepping variables
	*/

	ierr = SetupTimeFreqGrid(RSVDt, TR);CHKERRQ(ierr);

	/*
		Creates the discrete Fourier transform (DFT) and inverse DFT (iDFT) matrices
	*/

	ierr = CreateDFTiDFTMats(RSVDt, DFT);CHKERRQ(ierr);

	/*
		Loads the weight and spatial matrices (if applicable)
	*/

	if (!TR->TransRun) ierr = ReadWeightMats(RSVDt, Weight_mat, dirs);CHKERRQ(ierr);

	/*
		Generates a random input matrix (or read in if desired)
	*/

	if (!TR->TransRun) ierr = CreateRandomMat(RSVD, RSVDt, dirs);CHKERRQ(ierr);	

	PetscFunctionReturn(0);

}


