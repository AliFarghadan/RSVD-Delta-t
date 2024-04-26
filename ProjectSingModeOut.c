
#include <ProjectSingModeOut.h>
#include <PermuteMat.h>
#include <MergeFreqs.h>
#include <SplitFreqs.h>
#include <ReversePermuteMat.h>

PetscErrorCode ProjectSingModeOut(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, Mat SingVec)
{

	PetscErrorCode       ierr;
	PetscInt             ik;
	PetscScalar          coeff;
	Vec                  V_vec, Y_temp;
	Mat                  Y_hat_tall, Y_hat_re;

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetBool(NULL,NULL,"-Singular_flg",&Singular_flg,&Singular_flg);CHKERRQ(ierr);

	/* Compute y = y - v(v^Ty) */

	if (Singular_flg) {

		ierr = MatCreate(PETSC_COMM_WORLD,&Y_hat_re);CHKERRQ(ierr);
		ierr = MatSetType(Y_hat_re,MATDENSE);CHKERRQ(ierr);
		ierr = PermuteMat(RSVD_mat->Y_hat, RSVDt->RSVD.Nw, Y_hat_re);CHKERRQ(ierr);
		ierr = MergeFreqs(Y_hat_re,RSVDt->RSVD.Nw,&Y_hat_tall);CHKERRQ(ierr);
		ierr = MatDenseGetColumnVecRead(SingVec,0,&V_vec);CHKERRQ(ierr);
		for (ik=0; ik<RSVDt->RSVD.k; ik++) {
			ierr = MatDenseGetColumnVecWrite(Y_hat_tall,ik,&Y_temp);CHKERRQ(ierr);
			ierr = VecDot(Y_temp,V_vec,&coeff);CHKERRQ(ierr);
			ierr = VecAXPY(Y_temp,-coeff,V_vec);CHKERRQ(ierr);
			ierr = MatDenseRestoreColumnVecWrite(Y_hat_tall,ik,&Y_temp);CHKERRQ(ierr);
		}
		ierr = MatDenseRestoreColumnVecRead(SingVec,0,&V_vec);CHKERRQ(ierr);

		ierr = SplitFreqs(Y_hat_tall,RSVDt->RSVD.Nw,&Y_hat_re);CHKERRQ(ierr);
		ierr = MatCreate(PETSC_COMM_WORLD,&RSVD_mat->Y_hat);CHKERRQ(ierr);
		ierr = MatSetType(RSVD_mat->Y_hat,MATDENSE);CHKERRQ(ierr);
		ierr = ReversePermuteMat(Y_hat_re,RSVDt->RSVD.Nw,RSVD_mat->Y_hat);CHKERRQ(ierr);

	}

	PetscFunctionReturn(0);
	
}
