
#include <petscmat.h>
#include "Variables.h"

PetscErrorCode LNSStructure(LNS_vars *LNS_mat, RSVDt_vars *RSVDt, TransRun_vars *TR_vars, Directories *dirs)
{
	/*
		Reads the LNS operator and applies discounting if desired
	*/

	PetscErrorCode        ierr;
	PetscLogDouble        t1, t2;
	PetscViewer           fd;
	Mat                   A;
	PetscReal             norm;

	PetscFunctionBeginUser;

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Loading the operator: %s ***\n",dirs->OperatorDir);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->RootDir,dirs->OperatorDir);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&LNS_mat->A);CHKERRQ(ierr);
	ierr = MatSetType(LNS_mat->A,MATMPIAIJ);CHKERRQ(ierr);
	ierr = MatLoad(LNS_mat->A,fd);CHKERRQ(ierr);
	ierr = MatGetSize(LNS_mat->A,&RSVDt->RSVD.N,NULL);CHKERRQ(ierr);

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"*** The operator is loaded in %f seconds (N = %d) ***\n", t2-t1, (int)RSVDt->RSVD.N);CHKERRQ(ierr);

	/*
		Check whether the operator is real-valued or complex-valued. This determines the number of frequencies.
		If real-valued: only positive frequencies are considered, otherwise both negative and positive 
	*/

	ierr = MatDuplicate(LNS_mat->A,MAT_COPY_VALUES,&A);CHKERRQ(ierr);
	ierr = MatImaginaryPart(A);CHKERRQ(ierr);
	ierr = MatNorm(A,NORM_FROBENIUS,&norm);CHKERRQ(ierr);
	RSVDt->TS.RealOperator = norm < 1e-12 ? PETSC_TRUE : PETSC_FALSE;
	if (!TR_vars->TransRun) ierr = RSVDt->TS.RealOperator ? PetscPrintf(PETSC_COMM_WORLD,"\nThe LNS operator is real-valued!\n") : \
		PetscPrintf(PETSC_COMM_WORLD,"\nThe LNS operator is complex-valued!\n");CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	RSVDt->RSVD.Nw_eff = RSVDt->TS.RealOperator ? RSVDt->RSVD.Nw/2 : RSVDt->RSVD.Nw;
	if (!TR_vars->TransRun) ierr = RSVDt->TS.RealOperator ? PetscPrintf(PETSC_COMM_WORLD,\
					"The effective number of frequencies, Nw = %d\n",(int)RSVDt->RSVD.Nw_eff) : \
		PetscPrintf(PETSC_COMM_WORLD,"The effective number of frequencies, Nw = %d\n",(int)RSVDt->RSVD.Nw_eff);CHKERRQ(ierr);

	/*
		Discounting
	*/

	if (LNS_mat->RSVDt.Disc.DiscFlg) {
		ierr = MatShift(LNS_mat->A,-LNS_mat->RSVDt.Disc.beta);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"\n---- Discounting with beta = %g ----\n\n", LNS_mat->RSVDt.Disc.beta);CHKERRQ(ierr);
	}

	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}

