
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode LNSStructure(LNS_vars *LNS, RSVDt_vars *RSVDt, TransRun_vars *TR, Directories *dirs)
{
	/*
		Reads the LNS operator and applies discounting if desired
	*/

	PetscErrorCode        ierr;
	PetscLogDouble        t1, t2;
	PetscViewer           fd;
	Mat                   A;
	PetscReal             norm;
	PetscInt              hh,mm,ss;

	PetscFunctionBeginUser;

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"Loading the operator: %s%s\n",dirs->RootDir, dirs->OperatorDir);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->RootDir, dirs->OperatorDir);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&LNS->A);CHKERRQ(ierr);
	ierr = MatSetType(LNS->A,MATMPIAIJ);CHKERRQ(ierr);
	ierr = MatLoad(LNS->A,fd);CHKERRQ(ierr);
	ierr = MatGetSize(LNS->A,&RSVDt->RSVD.N,NULL);CHKERRQ(ierr);

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"Loaded operator elapsed time = %02d:%02d:%02d (N = %d)\n\n", (int)hh, (int)mm, (int)ss, (int)RSVDt->RSVD.N);CHKERRQ(ierr);

	/*
		Checks whether the operator is real-valued or complex-valued. This determines the number of frequencies.
		If real-valued: only positive (and zero) frequencies are considered, otherwise both negative and positive 
	*/

	ierr = MatDuplicate(LNS->A,MAT_COPY_VALUES,&A);CHKERRQ(ierr);
	ierr = MatImaginaryPart(A);CHKERRQ(ierr);
	ierr = MatNorm(A,NORM_FROBENIUS,&norm);CHKERRQ(ierr);
	RSVDt->TS.RealOperator = norm < 1e-12 ? PETSC_TRUE : PETSC_FALSE;
	ierr = MatDestroy(&A);CHKERRQ(ierr);

	/*
		Applies discounting and exits
	*/

	if (LNS->RSVDt.Disc.DiscFlg) {
		ierr = MatShift(LNS->A,-LNS->RSVDt.Disc.beta);CHKERRQ(ierr);
		if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n---- Discounting with beta = %g ----\n\n", LNS->RSVDt.Disc.beta);CHKERRQ(ierr);
	}

	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}

