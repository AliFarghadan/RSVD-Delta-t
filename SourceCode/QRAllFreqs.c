
#include <slepcbv.h>
#include "Variables.h"

PetscErrorCode QRAllFreqs(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt)
{
	/*
		Performs QR decomposition for all Nw matrices of size N \times k
	*/  

	PetscErrorCode        ierr;
	PetscInt              iw,hh,mm,ss;
	Mat                   Y_temp, Q_temp;
	BV                    Q;
	PetscLogDouble        t1, t2;

	PetscFunctionBeginUser;

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** QR decomposition begins! ***\n");CHKERRQ(ierr);

	for (iw=0; iw<RSVDt->RSVD.Nw_eff; iw++) {
		ierr = MatDenseGetSubMatrix(RSVD_mat->Y_hat,PETSC_DECIDE,PETSC_DECIDE,iw*RSVDt->RSVD.k,(iw+1)*RSVDt->RSVD.k,&Y_temp);CHKERRQ(ierr);
		ierr = BVCreateFromMat(Y_temp,&Q);CHKERRQ(ierr);
		ierr = BVSetType(Q,BVVECS);CHKERRQ(ierr);
		ierr = BVOrthogonalize(Q,NULL);CHKERRQ(ierr);
		ierr = BVCreateMat(Q,&Q_temp);CHKERRQ(ierr);
		ierr = MatCopy(Q_temp,Y_temp,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
		ierr = MatDenseRestoreSubMatrix(RSVD_mat->Y_hat,&Y_temp);CHKERRQ(ierr);
		ierr = BVDestroy(&Q);CHKERRQ(ierr);
		ierr = MatDestroy(&Q_temp);CHKERRQ(ierr);
	}

	/*
		Printing out the elapsed time and exit
	*/

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"*** QR decomposition elapsed time = %02d:%02d:%02d ***\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}



