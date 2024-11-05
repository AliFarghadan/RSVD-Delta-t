
#include <petscmat.h>
#include <Variables.h>
#include <CreateMtilde.h>
#include <GalerkinProjAllModes.h>
#include <ReversePermuteMat.h>

PetscErrorCode TransientRemovalStrategy(Mat Y_all, LNS_vars *LNS_mat, RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, TS_removal_matrices *TSR, DFT_matrices *DFT_mat)
{

	/*
		This code follows the efficient transient removal strategy, explained in our paper
	*/

	PetscErrorCode        ierr;
	PetscInt              hh,mm,ss;
	PetscLogDouble        t1, t2;

	PetscFunctionBeginUser;

	ierr = PetscTime(&t1);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Trasnient removal strategy begins! ***\n"
					"The process takes about as long as one period integration\n");CHKERRQ(ierr);

	/*
		Creates M_tilde for all test vectors
	*/

	ierr = CreateMtilde(Y_all,LNS_mat,RSVDt,TSR);CHKERRQ(ierr);

	/*
		Updates the steady-state solutions by removing the transient part using Galerkin projection
	*/

	ierr = GalerkinProjAllModes(Y_all,LNS_mat,RSVDt,TSR,DFT_mat);CHKERRQ(ierr);

	ierr = ReversePermuteMat(TSR->Y_hat_up, RSVDt);CHKERRQ(ierr);
	ierr = MatDuplicate(TSR->Y_hat_up,MAT_COPY_VALUES,&RSVD_mat->Y_hat);CHKERRQ(ierr);
	ierr = MatDestroy(&TSR->Y_hat_up);CHKERRQ(ierr);

	/*
		Prints out the elapsed time and exits
	*/

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Trasnient removal strategy elapsed time = %02d:%02d:%02d\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}



