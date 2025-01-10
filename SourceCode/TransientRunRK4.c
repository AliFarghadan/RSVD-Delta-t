
#include <petscmat.h>
#include <Variables.h>
#include <CreateTSMats.h>
#include <TSTransRK4.h>
#include <TransientRemovalEst.h>
#include <DestroyTSMats.h>
#include <SaveInputVarsCopy.h>

PetscErrorCode TransientRunRK4(TransRun_vars *TR, RSVDt_vars *RSVDt, LNS_vars *LNS, DFT_matrices *DFT, Directories *dirs)
{
	/*
		Transient simulation using RK4
		This function runs only if 'TransRun' is True
		Solves the initial value problem for the differential equation \dot{q} - Aq = 0 in the time domain
		This solver starts from a specified or random initial condition and computes the solution trajectory of q(t) over time
	*/

	PetscErrorCode        ierr;
	TS_matrices           TS_mat;
	Mat                   Q_all;
	Vec                   q0,trans_norm,q_temp,qss;
	PetscInt              Nt_saved,i,rend,pos=0,Nt_period,Ns,Nt_saved_delta,hh,mm,ss;
	PetscReal             norm,TSS,deltaT;
	PetscRandom           r;
	PetscViewer           fd;
	PetscLogDouble        t1, t2;

	PetscFunctionBeginUser;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n********************************************\n"
			"************** Transient run ***************\n********************************************\n\n");CHKERRQ(ierr);

	/*
		Define time-stepping (TS) parameters
	*/

	Ns               = RSVDt->RSVD.Nw;
	TSS              = 2*PETSC_PI/RSVDt->RSVD.w;
	deltaT           = TSS/Ns;
	RSVDt->TS.dt     = deltaT/PetscCeilReal(deltaT/RSVDt->TS.dt);
	Nt_period        = round(TSS/RSVDt->TS.dt);
	Nt_saved_delta   = Nt_period/Ns;
	rend             = Nt_period*TR->TransPeriods;
	Nt_saved         = PetscFloorReal(rend/TR->TransSaveMod)+1;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Total integration time = %g = %d x %g with dt = %g\n", \
					RSVDt->TS.dt*rend, (int) TR->TransPeriods, TSS, RSVDt->TS.dt);CHKERRQ(ierr);
	if (TR->TransRemovalEst) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Transient removal strategy estimation is requested\n\n");CHKERRQ(ierr);
	} else {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Transient removal strategy estimation is NOT requested\n\n");CHKERRQ(ierr);
	}

	/*
		Creates the required matrices/vecs
	*/

	ierr = TR->TransSave ? PetscPrintf(PETSC_COMM_WORLD,"Total snapshots to be saved = %d\n", (int) Nt_saved) : \
										PetscPrintf(PETSC_COMM_WORLD,"Saving snapshots is not requested\n");CHKERRQ(ierr);  

	if (TR->TransRemovalEst) {
		ierr = MatCreate(PETSC_COMM_WORLD,&Q_all);CHKERRQ(ierr);
		ierr = MatSetType(Q_all,MATDENSE);CHKERRQ(ierr);
		ierr = MatSetSizes(Q_all,PETSC_DECIDE,PETSC_DECIDE,RSVDt->RSVD.N,Ns);CHKERRQ(ierr);
		ierr = MatSetUp(Q_all);CHKERRQ(ierr);

		ierr = VecCreate(PETSC_COMM_WORLD, &qss);CHKERRQ(ierr);
		ierr = VecSetSizes(qss, PETSC_DECIDE, RSVDt->RSVD.N);CHKERRQ(ierr);
		ierr = VecSetUp(qss);CHKERRQ(ierr);
	}

	ierr = VecCreate(PETSC_COMM_WORLD,&trans_norm);CHKERRQ(ierr);
	ierr = VecSetType(trans_norm,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(trans_norm,PETSC_DECIDE,Nt_saved);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&q0);CHKERRQ(ierr);
	ierr = VecSetType(q0,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(q0,PETSC_DECIDE,RSVDt->RSVD.N);CHKERRQ(ierr);

	ierr = CreateTSMats(&TS_mat,RSVDt->RSVD.N);CHKERRQ(ierr);

	/*
		Reads in the initial condition or generates it randomly
		The initial condition is normalized unless the initial vector is specified
	*/

	if (TR->TransICFlg) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading the initial vector\n");CHKERRQ(ierr);
		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%s",dirs->RootDir,dirs->ResultsDir,dirs->TransICDir);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = VecLoad(q0,fd);CHKERRQ(ierr);
	} else {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Random vector is generated\n\n");CHKERRQ(ierr);
		ierr = PetscRandomCreate(PETSC_COMM_WORLD,&r);CHKERRQ(ierr);
		ierr = PetscRandomSetSeed(r, RSVDt->RSVD.RandSeed);CHKERRQ(ierr);
		ierr = PetscRandomSeed(r);CHKERRQ(ierr);
		ierr = VecSetRandom(q0,r);CHKERRQ(ierr);
		ierr = VecNormalize(q0,NULL);CHKERRQ(ierr);
	}

	if (TR->TransRemovalEst) {
		ierr = VecCopy(q0, qss);CHKERRQ(ierr);
		ierr = VecScale(qss, -1.);CHKERRQ(ierr);
	}

	/*
		Initializes the norm of transient snapshots
		The first saved snapshot is the random or the IC snapshot (if applicable)
	*/ 

	ierr = VecNorm(q0,NORM_2,&norm);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of the initial snapshot: %f @ t = %d\n", (double)norm, 0);CHKERRQ(ierr);
	ierr = VecSetValues(trans_norm,1,&pos,(PetscScalar*)&norm,INSERT_VALUES);CHKERRQ(ierr);

	if (TR->TransRemovalEst) {
		ierr = MatDenseGetColumnVecWrite(Q_all,0,&q_temp);CHKERRQ(ierr);
		ierr = VecCopy(q0,q_temp);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(Q_all,0,&q_temp);CHKERRQ(ierr);
	}

	if (TR->TransSave) {
		ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d",dirs->FolderDir,"q_transient_",(int) pos);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
		ierr = VecView(q0,fd);CHKERRQ(ierr);
	}

	/*
		Performs time stepping
	*/

	ierr = PetscTime(&t1);CHKERRQ(ierr);

	for (i=1; i<=rend; i++) {

		ierr = TSTransRK4(LNS,&TS_mat,RSVDt,q0);CHKERRQ(ierr);

		/*
			Computes the norm of snapshot and check the convergence/divergence status
			Saves the snapshot if desired
		*/

		if (PetscFmodReal(i,TR->TransSaveMod) == 0) { 

			pos  = i/TR->TransSaveMod;

			ierr = VecNorm(q0,NORM_2,&norm);CHKERRQ(ierr);

			if (PetscIsNanReal(norm)) {
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Stopped! The norm of the last snapshot is NAN!\n");CHKERRQ(ierr);
				break;
			}

			if (norm > TR->TransDivVal) {
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Diverged! The norm of the last snapshot = %g exceeds the divergence threshold = %g\n",\
																	norm, TR->TransDivVal);CHKERRQ(ierr);
				break;
			}

			if (norm < TR->TransConVal) {
				ierr = PetscPrintf(PETSC_COMM_WORLD,"Converged! The norm of the last snapshot = %g is smaller than the convergence threshold = %g\n",\
																	norm, TR->TransConVal);CHKERRQ(ierr);
				break;
			}

			if (TR->TransSave) {
				ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d",dirs->FolderDir,"q_transient_",(int) pos);CHKERRQ(ierr);
				ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
				ierr = VecView(q0,fd);CHKERRQ(ierr);
			}

			ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of the snapshot: %f @ t = %f\n",norm,i*RSVDt->TS.dt);CHKERRQ(ierr);
			ierr = VecSetValues(trans_norm,1,&pos,(PetscScalar*)&norm,INSERT_VALUES);CHKERRQ(ierr);

		}

		/*
			Estimates the transient removal performance (if applicable)
			The estimate occurs at the end of each period
		*/

		if (TR->TransRemovalEst) {
			if (PetscFmodReal(i,Nt_saved_delta) == 0) {
				pos  = i/Nt_saved_delta;
				pos  = PetscFmodReal(pos,Ns);
				if (pos == 0) ierr = TransientRemovalEst(Q_all, qss, i/Nt_period, RSVDt, DFT, dirs);CHKERRQ(ierr);
				ierr = MatDenseGetColumnVecWrite(Q_all,pos,&q_temp);CHKERRQ(ierr);
				ierr = VecCopy(q0,q_temp);CHKERRQ(ierr);
				ierr = MatDenseRestoreColumnVecWrite(Q_all,pos,&q_temp);CHKERRQ(ierr);
			}
		}
	}

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Transient simulation elapsed time = %02d:%02d:%02d ***\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);

	/*
		Saves the norm of snapshots throughout the intergation (every "TransSaveMod" number)
		Saves the last snapshot
	*/ 

	ierr = PetscPrintf(PETSC_COMM_WORLD,"\nThe results directory is: %s\n\n",dirs->FolderDir);CHKERRQ(ierr);	

	ierr = VecAssemblyBegin(trans_norm);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(trans_norm);CHKERRQ(ierr);
	ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->FolderDir,"q_transient_norms");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = VecView(trans_norm,fd);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->FolderDir,"q_transient_last_snapshot");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = VecView(q0,fd);CHKERRQ(ierr);

	/*
		Removes the variables from memory
		Saves a copy of the input variables in the results folder and exit
	*/

	ierr = VecDestroy(&q0);CHKERRQ(ierr);
	if (TR->TransRemovalEst) ierr = MatDestroy(&Q_all);CHKERRQ(ierr);
	if (TR->TransRemovalEst) ierr = VecDestroy(&qss);CHKERRQ(ierr);
	ierr = VecDestroy(&trans_norm);CHKERRQ(ierr);
	ierr = DestroyTSMats(&TS_mat);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"DONE :))\n\n");CHKERRQ(ierr);

	PetscFunctionReturn(0);

}

