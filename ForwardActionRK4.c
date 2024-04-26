
#include <ForwardActionRK4.h>
#include <TSActionRK4.h>
#include <ProjectSingModeOut.h>
#include <ApplyWSMats.h>

PetscErrorCode ForwardActionRK4(RSVD_matrices *RSVD_mat, RSVDt_vars *RSVDt, LNS_params *LNS_mat, \
			TS_matrices *TS_mat, DFT_matrices *DFT_mat, Sing_matrices *Sing_mat, Directories *dirs)
{

	PetscErrorCode        ierr;

	PetscFunctionBeginUser;
		
	RSVDt.TS.flg_dir_adj = 1;

	ierr = ProjectSingModeOut(RSVD_mat, RSVDt, Sing_mat->SingU);CHKERRQ(ierr);
	ierr = ApplyWSMats(RSVD_mat, RSVDt, WS_mat, 1)
	ierr = TSActionRK4(RSVD_mat, DFT_mat, LNS_mat, TS_mat, RSVDt, dirs);CHKERRQ(ierr);
	ierr = ProjectSingModeOut(RSVD_mat, RSVDt, Sing_mat->SingV);CHKERRQ(ierr);	
	ierr = ApplyWSMats(RSVD_mat, RSVDt, WS_mat, 0)

	PetscFunctionReturn(0);
	
}

