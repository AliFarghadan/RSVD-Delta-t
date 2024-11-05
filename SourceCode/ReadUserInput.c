
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode ReadUserInput(RSVDt_vars *RSVDt, Weight_matrices *Weight_mat, LNS_vars *LNS_mat, \
							Resolvent_matrices *Res_mat, TransRun_vars *TR_vars, Directories *dirs)
{
	/*
		Reads in user input parameters
	*/

	PetscErrorCode        ierr;
	PetscBool             flg_set;
	char                  filename[PETSC_MAX_PATH_LEN];

	PetscFunctionBeginUser;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** Reading user inputs ***\n");CHKERRQ(ierr);

	ierr = PetscOptionsGetString(NULL, NULL,"-inputs",(char*)&filename,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify inputs variables via -inputs");CHKERRQ(ierr);
	ierr = PetscOptionsInsertFileYAML(PETSC_COMM_WORLD, NULL, filename, PETSC_FALSE);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,NULL,"-TransientLength",&RSVDt->TS.TransientLength,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'TransientLength'");CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&RSVDt->TS.dt,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'dt'");CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,NULL,"-w",&RSVDt->RSVD.w,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'dt'");CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-InputForcingFlg",&RSVDt->RSVD.InputForcingFlg,&flg_set);CHKERRQ(ierr);	
	if (!flg_set) {
		RSVDt->RSVD.InputForcingFlg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'InputForcingFlg' variable not found. Setting 'InputForcingFlg' to default value: %d\n", (int) RSVDt->RSVD.InputForcingFlg);
	} else {
		if (RSVDt->RSVD.InputForcingFlg) {
			ierr = PetscOptionsGetString(NULL,NULL,"-InputForcingDir",(char*)&dirs->InputForcingDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
			if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'InputForcingDir'");CHKERRQ(ierr);
		}
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-InvInputWeightFlg",&Weight_mat->InvInputWeightFlg,&flg_set);CHKERRQ(ierr);	
	if (!flg_set) {
		Weight_mat->InvInputWeightFlg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'InvInputWeightFlg' variable not found. Setting 'InputWeightFlg' to default value: %d\n", (int) Weight_mat->InvInputWeightFlg);
	} else {
		if (Weight_mat->InvInputWeightFlg) {
			ierr = PetscOptionsGetString(NULL,NULL,"-InvInputWeightDir",(char*)&dirs->InvInputWeightDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
			if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'InvInputWeightDir'");CHKERRQ(ierr);
		}
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-OutputWeightFlg",&Weight_mat->OutputWeightFlg,&flg_set);CHKERRQ(ierr);	
	if (!flg_set) {
		Weight_mat->OutputWeightFlg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'OutputWeightFlg' variable not found. Setting 'OutputWeightFlg' to default value: %d\n", (int) Weight_mat->OutputWeightFlg);
	} else {
		if (Weight_mat->OutputWeightFlg) {
			ierr = PetscOptionsGetString(NULL,NULL,"-OutputWeightDir",(char*)&dirs->OutputWeightDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
			if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'OutputWeightDir'");CHKERRQ(ierr);
		}
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-InvOutputWeightFlg",&Weight_mat->InvOutputWeightFlg,&flg_set);CHKERRQ(ierr);	
	if (!flg_set) {
		Weight_mat->InvOutputWeightFlg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'InvOutputWeightFlg' variable not found. Setting 'InvOutputWeightFlg' to default value: %d\n", (int) Weight_mat->InvOutputWeightFlg);
	} else {
		if (Weight_mat->InvOutputWeightFlg) {
			ierr = PetscOptionsGetString(NULL,NULL,"-InvOutputWeightDir",(char*)&dirs->InvOutputWeightDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
			if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'InvOutputWeightDir'");CHKERRQ(ierr);
		}
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-InputMatrixFlg",&Weight_mat->InputMatrixFlg,&flg_set);CHKERRQ(ierr);	
	if (!flg_set) {
		Weight_mat->InputMatrixFlg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'InputMatrixFlg' variable not found. Setting 'InputMatrixFlg' to default value: %d\n", (int) Weight_mat->InputMatrixFlg);
	} else {
		if (Weight_mat->InputMatrixFlg) {
			ierr = PetscOptionsGetString(NULL,NULL,"-InputMatrixDir",(char*)&dirs->InputMatrixDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
			if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'InputMatrixDir'");CHKERRQ(ierr);
		}
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-OutputMatrixFlg",&Weight_mat->OutputMatrixFlg,&flg_set);CHKERRQ(ierr);	
	if (!flg_set) {
		Weight_mat->OutputMatrixFlg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'OutputMatrixFlg' variable not found. Setting 'OutputMatrixFlg' to default value: %d\n", (int) Weight_mat->OutputMatrixFlg);
	} else {
		if (Weight_mat->OutputMatrixFlg) {
			ierr = PetscOptionsGetString(NULL,NULL,"-OutputMatrixDir",(char*)&dirs->OutputMatrixDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
			if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'OutputMatrixDir'");CHKERRQ(ierr);
		}
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-TwoPI",&RSVDt->RSVD.TwoPI,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVDt->RSVD.TwoPI = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TwoPI' variable not found. Setting 'TwoPI' to default value: %d\n", (int) RSVDt->RSVD.TwoPI);
	}
	RSVDt->RSVD.w *= (RSVDt->RSVD.TwoPI) ? 2*PETSC_PI : 1;
	ierr = PetscOptionsGetBool(NULL,NULL,"-TransRun",&TR_vars->TransRun,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		TR_vars->TransRun = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransRun' variable not found. Setting 'TransRun' to default value: %d\n", (int) TR_vars->TransRun);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-TransientRemoval",&RSVDt->TS.TransientRemoval,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVDt->TS.TransientRemoval = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransientRemoval' variable not found. Setting 'TransientRemoval' to default value: %d\n", (int) RSVDt->TS.TransientRemoval);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-DiscFlg",&LNS_mat->RSVDt.Disc.DiscFlg,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		LNS_mat->RSVDt.Disc.DiscFlg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'DiscFlg' variable not found. Setting 'DiscFlg' to default value: %d\n", (int) LNS_mat->RSVDt.Disc.DiscFlg);
	}
	ierr = PetscOptionsGetInt(NULL,NULL,"-SaveResultsOpt",&RSVDt->SaveResultsOpt,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVDt->SaveResultsOpt = 1;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'SaveResultsOpt' variable not found. Setting 'SaveResultsOpt' to default value: %d\n", (int) RSVDt->SaveResultsOpt);
	}	
	ierr = PetscOptionsGetString(NULL, NULL,"-RootDir",(char*)&dirs->RootDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'RootDir'");CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL,"-ResultsDir",(char*)&dirs->ResultsDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'ResultsDir'");CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL,"-OperatorDir",(char*)&dirs->OperatorDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'OperatorDir'");CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-k",&RSVDt->RSVD.k,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVDt->RSVD.k = 3;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'k' variable not found. Setting 'k' to default value: %d\n", (int) RSVDt->RSVD.k);
	}
	ierr = PetscOptionsGetInt(NULL,NULL,"-q",&RSVDt->RSVD.q,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVDt->RSVD.q = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'q' variable not found. Setting 'q' to default value: %d\n", (int) RSVDt->RSVD.q);
	}
	ierr = PetscOptionsGetInt(NULL,NULL,"-Nw",&RSVDt->RSVD.Nw,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'Nw'");CHKERRQ(ierr);	
	ierr = PetscOptionsGetInt(NULL,NULL,"-RandSeed",&RSVDt->RSVD.RandSeed,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVDt->RSVD.RandSeed = 1373;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'RandSeed' variable not found. Setting 'RandSeed' to default value: %d\n", (int) RSVDt->RSVD.RandSeed);
	}	
	ierr = PetscOptionsGetInt(NULL,NULL,"-Display",&RSVDt->Display,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVDt->Display = 2;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'Display' variable not found. Setting 'Display' to default value: %d\n", (int) RSVDt->Display);
	}
	ierr = PetscOptionsGetReal(NULL,NULL,"-beta",&LNS_mat->RSVDt.Disc.beta,&flg_set);CHKERRQ(ierr);
	if (!flg_set && LNS_mat->RSVDt.Disc.DiscFlg) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Discounting flag is on! Either set Discounting to zero or specify 'beta'");
	if (LNS_mat->RSVDt.Disc.beta < 0 && LNS_mat->RSVDt.Disc.DiscFlg) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"'beta' must be positive, current value: %g", LNS_mat->RSVDt.Disc.beta);
	ierr = PetscOptionsGetBool(NULL,NULL,"-TransICFlg",&TR_vars->TransICFlg,&flg_set);CHKERRQ(ierr);
	if (!flg_set && TR_vars->TransRun) {
		TR_vars->TransICFlg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransICFlg' variable not found. Setting 'TransICFlg' to default value: %d\n", (int) TR_vars->TransICFlg);
	} else {
		ierr = PetscOptionsGetString(NULL,NULL,"-TransICDir",(char*)&dirs->TransICDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
		if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'TransICDir'");CHKERRQ(ierr);
	}
	ierr = PetscOptionsGetReal(NULL,NULL,"-TransDivVal",&TR_vars->TransDivVal,&flg_set);CHKERRQ(ierr);
	if (!flg_set && TR_vars->TransRun) {
		TR_vars->TransDivVal = 1e3;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransDivVal' variable not found. Setting 'TransDivVal' to default value: %g\n", TR_vars->TransDivVal);
	}
	ierr = PetscOptionsGetReal(NULL,NULL,"-TransConVal",&TR_vars->TransConVal,&flg_set);CHKERRQ(ierr);
	if (!flg_set && TR_vars->TransRun) {
		TR_vars->TransConVal = 1e-8;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransConVal' variable not found. Setting 'TransConVal' to default value: %g\n", TR_vars->TransConVal);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-TransSave",&TR_vars->TransSave,&flg_set);CHKERRQ(ierr);
	if (!flg_set && TR_vars->TransRun) {
		TR_vars->TransSave = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransSave' variable not found. Setting 'TransSave' to default value: %d\n", (int) TR_vars->TransSave);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-TransRemovalEst",&TR_vars->TransRemovalEst,&flg_set);CHKERRQ(ierr);
	if (!flg_set && TR_vars->TransRun) {
		TR_vars->TransRemovalEst = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransRemovalEst' variable not found. Setting 'TransRemovalEst' to default value: %d\n", (int) TR_vars->TransRemovalEst);
	}
	ierr = PetscOptionsGetInt(NULL,NULL,"-TransSaveMod",&TR_vars->TransSaveMod,&flg_set);CHKERRQ(ierr);
	if (!flg_set && TR_vars->TransRun) {
		TR_vars->TransSaveMod = 100;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransSaveMod' variable not found. Setting 'TransSaveMod' to default value: %d\n", (int) TR_vars->TransSaveMod);
	}
	ierr = PetscOptionsGetInt(NULL,NULL,"-TransPeriods",&TR_vars->TransPeriods,&flg_set);CHKERRQ(ierr);
	if (!flg_set && TR_vars->TransRun) {
		TR_vars->TransPeriods = 1;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransPeriods' variable not found. Setting 'TransPeriods' to default value: %d\n", (int) TR_vars->TransPeriods);
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);

	PetscFunctionReturn(0);

}

