
#include <petscksp.h>
#include "Variables.h"
#include "QR_simple.h"

PetscErrorCode EfficientTransientRemovalEstimation(Mat Q_transient, Vec qss, PetscInt period_index, RSVDt_vars *RSVDt, DFT_matrices *DFT_mat, Directories *dirs)
{
	/*
		Once the matrix of transient responses is created, we can test the transinet removal performance.
		This function assumes Q_all is periodic and saved in the time domain.
	*/  
	
	PetscErrorCode        ierr;
	PetscReal             w, norm, norm_diff, deltaT;
	PetscInt              N, Ns, is, is1, iw;
	Vec                   Initial_transient_norm, Updated_transient_norm, q_temp, q1, q2, b, qdiff;
	Mat                   Q_all, Q_transient_hat;
	PetscViewer           fd;

	PetscFunctionBeginUser;

	/*
		Display the information of the transient removal strategy
	*/

	ierr = MatGetSize(Q_transient,&N,&Ns);CHKERRQ(ierr);
	deltaT  = 2*PETSC_PI/RSVDt->RSVD.w/Ns;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"**** deltaT = %g ****\n", (double) deltaT);
	Ns  -= period_index == 1 ? 1 : 0; // skip the initial snapshot as it is a zero vector
	ierr = MatCreate(PETSC_COMM_WORLD,&Q_all);CHKERRQ(ierr);
	ierr = MatSetType(Q_all,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Q_all,PETSC_DECIDE,PETSC_DECIDE,N,Ns);CHKERRQ(ierr);
	ierr = MatSetUp(Q_all);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"**** Subspace dim = %d ****\n", (int) Ns);

	ierr = VecCreate(PETSC_COMM_WORLD, &Initial_transient_norm);CHKERRQ(ierr);
	ierr = VecSetSizes(Initial_transient_norm, PETSC_DECIDE, RSVDt->RSVD.Nw/2);CHKERRQ(ierr);
	ierr = VecSetUp(Initial_transient_norm);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &Updated_transient_norm);CHKERRQ(ierr);
	ierr = VecSetSizes(Updated_transient_norm, PETSC_DECIDE, RSVDt->RSVD.Nw/2);CHKERRQ(ierr);
	ierr = VecSetUp(Updated_transient_norm);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &q_temp);CHKERRQ(ierr);
	ierr = VecSetSizes(q_temp, PETSC_DECIDE, N);CHKERRQ(ierr);
	ierr = VecSetUp(q_temp);CHKERRQ(ierr);

	ierr = MatCreateVecs(Q_all,&b,&qdiff);CHKERRQ(ierr);

	/*
		Obtain Q_transient_hat from Q_transient
		Check how the transient norm decays from a period to the next for all frequencies
	*/

	ierr = PetscPrintf(PETSC_COMM_WORLD,"**** period index = %d ****\n", (int) period_index);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"**** Displaying transient norms for (positive) frequencies ****\n");CHKERRQ(ierr);

	ierr = MatMatMult(Q_transient,DFT_mat->dft,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Q_transient_hat);CHKERRQ(ierr);
	ierr = MatScale(Q_transient_hat, 1./RSVDt->RSVD.Nw);CHKERRQ(ierr);

	for (iw=0; iw<RSVDt->RSVD.Nw/2; iw++) {
		ierr = MatDenseGetColumnVecWrite(Q_transient_hat,iw,&q1);CHKERRQ(ierr);
		ierr = VecNorm(q1,NORM_2,&norm);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"**** Norm of the ratio norm(qt)/norm(qs): %g @ w = %g ****\n", \
									(double) norm, (double) RSVDt->RSVD.w * iw);CHKERRQ(ierr);
		ierr = VecSetValues(Initial_transient_norm,1,&iw,(PetscScalar*)&norm,INSERT_VALUES);CHKERRQ(ierr);
		ierr = VecAssemblyBegin(Initial_transient_norm);CHKERRQ(ierr);
		ierr = VecAssemblyEnd(Initial_transient_norm);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(Q_transient_hat,iw,&q1);CHKERRQ(ierr);
	}

	ierr = MatDestroy(&Q_transient_hat);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"**** Displaying updated transient norms for (positive) frequencies ****\n");CHKERRQ(ierr);

	for (iw=0; iw<RSVDt->RSVD.Nw/2; iw++) {

		/*
			Build the matrix of synthetic snapshots 
			y = y_steadystate e^(i w \Delta t) + Q_transient
		*/

		w = RSVDt->RSVD.w * iw;

		for (is=0; is<Ns; is++) {
			is1  = period_index == 1 ? is+1 : is; // skip the initial snapshot as it is a zero vector
			ierr = MatDenseGetColumnVecWrite(Q_all,is,&q1);CHKERRQ(ierr);
			ierr = MatDenseGetColumnVecRead(Q_transient,is1,&q2);CHKERRQ(ierr);

			ierr = VecCopy(q2, q1);CHKERRQ(ierr);
			ierr = VecCopy(qss, q_temp);CHKERRQ(ierr);
			ierr = VecScale(q_temp, PetscExpComplex(PETSC_i * w * is * deltaT));CHKERRQ(ierr);
			ierr = VecAXPY(q1,1.,q_temp);CHKERRQ(ierr);
			
			ierr = MatDenseRestoreColumnVecRead(Q_transient,is1,&q2);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(Q_all,is,&q1);CHKERRQ(ierr);
		}

		/*
			Build the orthogonal subspace from matrix of snapshots (U = QR(Q_all))
			Compute qdiff = y_steadystate - U(U'*y_steadystate)
			Save the norm(qdiff)/norm(qss) --> norm closer to 0 means more accuracy (less transient)
		*/

		ierr = VecNorm(qss,NORM_2,&norm);CHKERRQ(ierr);
		ierr = QR_simple(Q_all);CHKERRQ(ierr);
		ierr = MatMultHermitianTranspose(Q_all,qss,b);CHKERRQ(ierr);
		ierr = MatMult(Q_all,b,qdiff);CHKERRQ(ierr);
		ierr = VecAYPX(qdiff,-1,qss);CHKERRQ(ierr);
		ierr = VecNorm(qdiff,NORM_2,&norm_diff);CHKERRQ(ierr);
		norm = norm_diff/norm;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"**** Norm of the ratio norm(qt_up)/norm(qs): %g @ w = %g ****\n",norm, w);CHKERRQ(ierr);
		ierr = VecSetValues(Updated_transient_norm,1,&iw,(PetscScalar*)&norm,INSERT_VALUES);CHKERRQ(ierr);
		ierr = VecAssemblyBegin(Updated_transient_norm);CHKERRQ(ierr);
		ierr = VecAssemblyEnd(Updated_transient_norm);CHKERRQ(ierr);
	
	}

	ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d",dirs->FolderDir,"Initial_transient_norm_period_",(int) period_index);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = VecView(Initial_transient_norm,fd);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d",dirs->FolderDir,"Updated_transient_norm_period_",(int) period_index);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = VecView(Updated_transient_norm,fd);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

