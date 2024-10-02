
#include <petscmat.h>
#include "Variables.h"
#include "CreateTSMats.h"
#include "PermuteMat.h"
#include "CreateForcingOnFly.h"
#include "TSRK4.h"
#include "SaveSnapshots.h"
#include "DisplayProgress.h"
#include "DestroyTSMats.h"
#include "TransientRemovalStrategy.h"
#include "DFT.h"


PetscErrorCode TSActionRK4(RSVD_matrices *RSVD_mat, DFT_matrices *DFT_mat, \
		LNS_vars *LNS_mat, RSVDt_vars *RSVDt, Directories *dirs, TS_removal_matrices *TSR)
{
	/*
		Time stepping of the LNS equations to obtain either action of direct or adjoint resolvent operator
	*/

	PetscErrorCode       ierr;
	TS_matrices          TS_mat;
	Mat                  Y_all,F_temp,Y_all_k,F_hat_k;
	Vec                  y,y_temp;
	PetscInt             N,k,Nw,Nstore,i,ik,rend,prg_cnt=0,hh,mm,ss;
	PetscLogDouble       t0,t1=0,t2;
	
	PetscFunctionBeginUser;

	ierr = RSVDt->TS.DirAdj ? PetscPrintf(PETSC_COMM_WORLD,"\n*** Direct action begins! ***\n") : \
				PetscPrintf(PETSC_COMM_WORLD,"\n*** Adjoint action begins! ***\n");CHKERRQ(ierr);

	/*
		Defines the time-stepping parameters
	*/

	N      = RSVDt->RSVD.N;
	k      = RSVDt->RSVD.k;
	Nw     = RSVDt->RSVD.Nw;
	rend   = RSVDt->TS.Nt + RSVDt->TS.Ns;
	rend  += RSVDt->TS.TransientRemoval ? RSVDt->TS.ResRatio : 0;
	Nstore = RSVDt->TS.TransientRemoval ? Nw+1 : Nw;
	
	/*
		Uses the solution from the previous action as forcing for the current action
	*/

	ierr = MatDuplicate(RSVD_mat->Y_hat,MAT_COPY_VALUES,&RSVD_mat->F_hat);CHKERRQ(ierr);
	ierr = MatDestroy(&RSVD_mat->Y_hat);CHKERRQ(ierr);

	// // zero the first few columns --test-- !!!!!!!!
	for (ik=0; ik<0; ik++) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Skipping the zeroth forcing!!! ***\n");CHKERRQ(ierr);
		ierr = MatDenseGetColumnVecWrite(RSVD_mat->F_hat,ik,&y_temp);CHKERRQ(ierr);
		ierr = VecScale(y_temp, 0.0);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(RSVD_mat->F_hat,ik,&y_temp);CHKERRQ(ierr);
	}

	/*
		Creates all required matrices
	*/

	ierr = VecCreate(PETSC_COMM_WORLD,&y);CHKERRQ(ierr);
	ierr = VecSetType(y,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(y,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = CreateTSMats(&TS_mat,N);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&Y_all_k);CHKERRQ(ierr);
	ierr = MatSetType(Y_all_k,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Y_all_k,PETSC_DECIDE,PETSC_DECIDE,N,Nstore);CHKERRQ(ierr);
	ierr = MatSetUp(Y_all_k);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&Y_all);CHKERRQ(ierr);
	ierr = MatSetType(Y_all,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Y_all,PETSC_DECIDE,PETSC_DECIDE,N,Nstore*k);CHKERRQ(ierr);
	ierr = MatSetUp(Y_all);CHKERRQ(ierr);

	/*
		Permutes the forcing matrix
	*/

	ierr = PermuteMat(RSVD_mat->F_hat, RSVDt);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&F_hat_k);CHKERRQ(ierr);
	ierr = MatSetType(F_hat_k,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(F_hat_k,PETSC_DECIDE,PETSC_DECIDE,N,RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);
	ierr = MatSetUp(F_hat_k);CHKERRQ(ierr);

	/*
		Extracts the forcing submatrix for each mode and loops over all modes
	*/  

	ierr = PetscTime(&t0);CHKERRQ(ierr);

	for (ik=0; ik<k; ik++) {

		ierr = VecZeroEntries(y);CHKERRQ(ierr);

		ierr = MatDenseGetSubMatrix(RSVD_mat->F_hat,PETSC_DECIDE,PETSC_DECIDE,ik*RSVDt->RSVD.Nw_eff,(ik+1)*RSVDt->RSVD.Nw_eff,&F_temp);CHKERRQ(ierr);
		ierr = MatCopy(F_temp,F_hat_k,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		ierr = MatAssemblyBegin(F_hat_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(F_hat_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(RSVD_mat->F_hat,&F_temp);CHKERRQ(ierr);

		/*
			Creates the initial forcing
		*/

		ierr = CreateForcingOnFly(F_hat_k,DFT_mat,0,TS_mat.F3);CHKERRQ(ierr);

		/*
			Measures the time-stepping wall-time
		*/

		ierr = PetscTime(&t1);CHKERRQ(ierr);

		/*
			Performs time stepping
			Monitoring the progress by -display feature
		*/

		for (i=1; i<=rend; i++) {

			ierr = TSRK4(F_hat_k,DFT_mat,LNS_mat,&TS_mat,RSVDt,i,y);CHKERRQ(ierr);
			ierr = SaveSnapshots(y,i,RSVDt,Y_all_k);CHKERRQ(ierr);
			ierr = DisplayProgress(i,rend,ik,&prg_cnt,t1,RSVDt,1);CHKERRQ(ierr);

		}

		ierr = DisplayProgress(i,rend,ik,&prg_cnt,t1,RSVDt,0);CHKERRQ(ierr);

		ierr = MatAssemblyBegin(Y_all_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(Y_all_k,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatDenseGetSubMatrix(Y_all,PETSC_DECIDE,PETSC_DECIDE,ik*Nstore,(ik+1)*Nstore,&F_temp);CHKERRQ(ierr);
		ierr = MatCopy(Y_all_k,F_temp,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(Y_all,&F_temp);CHKERRQ(ierr);

	}

	/*
		Removes the forcing input as it is no longer needed after the repsonse is obtained
	*/

	ierr = VecDestroy(&y);CHKERRQ(ierr);
	ierr = MatDestroy(&RSVD_mat->F_hat);CHKERRQ(ierr);
	ierr = DestroyTSMats(&TS_mat);CHKERRQ(ierr);

	/*
		Takes the response to the frequency domain
		Efficient transient removal is performed before discrete Fourier transform (DFT) if desired
	*/

	ierr = RSVDt->TS.TransientRemoval ? TransientRemovalStrategy(Y_all,LNS_mat,RSVD_mat,RSVDt,TSR,DFT_mat) : \
											DFT(Y_all,DFT_mat,RSVD_mat,RSVDt);CHKERRQ(ierr);
	
	ierr = RSVDt->TS.DirAdj ? PetscPrintf(PETSC_COMM_WORLD,"*** Direct action ") : PetscPrintf(PETSC_COMM_WORLD,"*** Adjoint action ");CHKERRQ(ierr);

	/*
		Printing out the elapsed time and exit
	*/

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t0)/3600;
	mm   = (t2-t0-3600*hh)/60;
	ss   = t2-t0-3600*hh-mm*60;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"elapsed time = %02d:%02d:%02d ***\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}


