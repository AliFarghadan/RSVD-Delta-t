
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode ReadUserInput(RSVDt_vars *RSVDt, Weight_matrices *Weight, LNS_vars *LNS, \
							TransRun_vars *TR, Directories *dirs)
{
	/*
		Reads in user input parameters
	*/

	PetscErrorCode        ierr;
	PetscBool             flg_set;
	char                  filename[PETSC_MAX_PATH_LEN];

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetString(NULL, NULL,"-inputs",(char*)&filename,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify inputs variables via -inputs");CHKERRQ(ierr);
	ierr = PetscOptionsInsertFileYAML(PETSC_COMM_WORLD, NULL, filename, PETSC_FALSE);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-TransRun",&TR->TransRun,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		TR->TransRun = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransRun' variable not found. Setting 'TransRun' to default value: %d\n", (int) TR->TransRun);
	}	
	ierr = PetscOptionsGetString(NULL, NULL,"-RootDir",(char*)&dirs->RootDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'RootDir'");CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL,"-ResultsDir",(char*)&dirs->ResultsDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'ResultsDir'");CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL,"-OperatorDir",(char*)&dirs->OperatorDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'OperatorDir'");CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL,NULL,"-RandSeed",&RSVDt->RSVD.RandSeed,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVDt->RSVD.RandSeed = 1373;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'RandSeed' variable not found. Setting 'RandSeed' to default value: %d\n", (int) RSVDt->RSVD.RandSeed);
	}	
	ierr = PetscOptionsGetInt(NULL,NULL,"-Display",&RSVDt->Display,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVDt->Display = 2;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'Display' variable not found. Setting 'Display' to default value: %d\n", (int) RSVDt->Display);
	} else if (RSVDt->Display > 2 || RSVDt->Display < 0) {
		RSVDt->Display = 2;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'Display' must be 0, 1 or 2. Setting 'Display' to default value: %d\n", (int) RSVDt->Display);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-DiscFlg",&LNS->RSVDt.Disc.DiscFlg,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		LNS->RSVDt.Disc.DiscFlg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'DiscFlg' variable not found. Setting 'DiscFlg' to default value: %d\n", (int) LNS->RSVDt.Disc.DiscFlg);
	}	
	ierr = PetscOptionsGetReal(NULL,NULL,"-beta",&LNS->RSVDt.Disc.beta,&flg_set);CHKERRQ(ierr);
	if (!flg_set && LNS->RSVDt.Disc.DiscFlg) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Discounting flag is on! Either set 'DiscFlg' to zero or specify 'beta'");
	if (LNS->RSVDt.Disc.beta < 0 && LNS->RSVDt.Disc.DiscFlg) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"'beta' must be positive, current value: %g", LNS->RSVDt.Disc.beta);
	ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&RSVDt->TS.dt,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'dt'");CHKERRQ(ierr);
	} else if (RSVDt->TS.dt < 0) {
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"'dt' must be non-negative, current value: %g\n", RSVDt->TS.dt);CHKERRQ(ierr);
	}
	ierr = PetscOptionsGetReal(NULL,NULL,"-w",&RSVDt->RSVD.w,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'w'");CHKERRQ(ierr);
	} else if (RSVDt->RSVD.w <= 0) {
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"'w' must be positive, current value: %g\n", RSVDt->RSVD.w);CHKERRQ(ierr);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-TwoPI",&RSVDt->RSVD.TwoPI,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVDt->RSVD.TwoPI = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TwoPI' variable not found. Setting 'TwoPI' to default value: %d\n", (int) RSVDt->RSVD.TwoPI);
	}
	RSVDt->RSVD.w *= (RSVDt->RSVD.TwoPI) ? 2*PETSC_PI : 1;	
	ierr = PetscOptionsGetInt(NULL,NULL,"-Nw",&RSVDt->RSVD.Nw,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'Nw'");CHKERRQ(ierr);	
	} else if (RSVDt->RSVD.Nw < 1) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"'Nw' must be a positive integer, current value: %d\n", (int) RSVDt->RSVD.Nw);CHKERRQ(ierr);
	}
	if (!TR->TransRun) {	
		ierr = PetscOptionsGetReal(NULL,NULL,"-TransientLength",&RSVDt->TS.TransientLength,&flg_set);CHKERRQ(ierr);
		if (!flg_set) {
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'TransientLength'");CHKERRQ(ierr);
		} else if (RSVDt->TS.TransientLength < 0) {
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"'TransientLength' must be non-negative, current value: %g\n", RSVDt->TS.TransientLength);CHKERRQ(ierr);
		}
		ierr = PetscOptionsGetBool(NULL,NULL,"-InputForcingFlg",&RSVDt->RSVD.InputForcingFlg,&flg_set);CHKERRQ(ierr);	
		if (!flg_set) {
			RSVDt->RSVD.InputForcingFlg = 0;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'InputForcingFlg' variable not found. Setting 'InputForcingFlg' to default value: %d\n", (int) RSVDt->RSVD.InputForcingFlg);
		} else if (RSVDt->RSVD.InputForcingFlg) {
			ierr = PetscOptionsGetString(NULL,NULL,"-InputForcingDir",(char*)&dirs->InputForcingDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
			if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'InputForcingDir'");CHKERRQ(ierr);
		}
		ierr = PetscOptionsGetBool(NULL,NULL,"-InvInputWeightFlg",&Weight->InvInputWeightFlg,&flg_set);CHKERRQ(ierr);	
		if (!flg_set) {
			Weight->InvInputWeightFlg = 0;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'InvInputWeightFlg' variable not found. Setting 'InputWeightFlg' to default value: %d\n", (int) Weight->InvInputWeightFlg);
		} else if (Weight->InvInputWeightFlg) {
			ierr = PetscOptionsGetString(NULL,NULL,"-InvInputWeightDir",(char*)&dirs->InvInputWeightDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
			if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'InvInputWeightDir'");CHKERRQ(ierr);
		}
		ierr = PetscOptionsGetBool(NULL,NULL,"-OutputWeightFlg",&Weight->OutputWeightFlg,&flg_set);CHKERRQ(ierr);	
		if (!flg_set) {
			Weight->OutputWeightFlg = 0;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'OutputWeightFlg' variable not found. Setting 'OutputWeightFlg' to default value: %d\n", (int) Weight->OutputWeightFlg);
		} else if (Weight->OutputWeightFlg) {
			ierr = PetscOptionsGetString(NULL,NULL,"-OutputWeightDir",(char*)&dirs->OutputWeightDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
			if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'OutputWeightDir'");CHKERRQ(ierr);
		}
		ierr = PetscOptionsGetBool(NULL,NULL,"-InvOutputWeightFlg",&Weight->InvOutputWeightFlg,&flg_set);CHKERRQ(ierr);	
		if (!flg_set) {
			Weight->InvOutputWeightFlg = 0;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'InvOutputWeightFlg' variable not found. Setting 'InvOutputWeightFlg' to default value: %d\n", (int) Weight->InvOutputWeightFlg);
		} else if (Weight->InvOutputWeightFlg) {
			ierr = PetscOptionsGetString(NULL,NULL,"-InvOutputWeightDir",(char*)&dirs->InvOutputWeightDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
			if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'InvOutputWeightDir'");CHKERRQ(ierr);
		}
		ierr = PetscOptionsGetBool(NULL,NULL,"-InputMatrixFlg",&Weight->InputMatrixFlg,&flg_set);CHKERRQ(ierr);	
		if (!flg_set) {
			Weight->InputMatrixFlg = 0;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'InputMatrixFlg' variable not found. Setting 'InputMatrixFlg' to default value: %d\n", (int) Weight->InputMatrixFlg);
		} else if (Weight->InputMatrixFlg) {
			ierr = PetscOptionsGetString(NULL,NULL,"-InputMatrixDir",(char*)&dirs->InputMatrixDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
			if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'InputMatrixDir'");CHKERRQ(ierr);
		}
		ierr = PetscOptionsGetBool(NULL,NULL,"-OutputMatrixFlg",&Weight->OutputMatrixFlg,&flg_set);CHKERRQ(ierr);	
		if (!flg_set) {
			Weight->OutputMatrixFlg = 0;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'OutputMatrixFlg' variable not found. Setting 'OutputMatrixFlg' to default value: %d\n", (int) Weight->OutputMatrixFlg);
		} else if (Weight->OutputMatrixFlg) {
			ierr = PetscOptionsGetString(NULL,NULL,"-OutputMatrixDir",(char*)&dirs->OutputMatrixDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
			if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'OutputMatrixDir'");CHKERRQ(ierr);
		}
		ierr = PetscOptionsGetBool(NULL,NULL,"-TransientRemoval",&RSVDt->TS.TransientRemoval,&flg_set);CHKERRQ(ierr);
		if (!flg_set) {
			RSVDt->TS.TransientRemoval = 0;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransientRemoval' variable not found. Setting 'TransientRemoval' to default value: %d\n", (int) RSVDt->TS.TransientRemoval);
		}
		ierr = PetscOptionsGetInt(NULL,NULL,"-SaveResultsOpt",&RSVDt->SaveResultsOpt,&flg_set);CHKERRQ(ierr);
		if (!flg_set) {
			RSVDt->SaveResultsOpt = 1;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'SaveResultsOpt' variable not found. Setting 'SaveResultsOpt' to default value: %d\n", (int) RSVDt->SaveResultsOpt);
		} else if (RSVDt->SaveResultsOpt > 2 || RSVDt->SaveResultsOpt < 1) {
			RSVDt->SaveResultsOpt = 1;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'SaveResultsOpt' must be 1 or 2. Setting 'SaveResultsOpt' to default value: %d\n", (int) RSVDt->SaveResultsOpt);
		}
		ierr = PetscOptionsGetInt(NULL,NULL,"-k",&RSVDt->RSVD.k,&flg_set);CHKERRQ(ierr);
		if (!flg_set) {
			RSVDt->RSVD.k = 3;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'k' variable not found. Setting 'k' to default value: %d\n", (int) RSVDt->RSVD.k);
		} else if (RSVDt->RSVD.k < 1) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'k' must be a positive integer, current value: %d\n", (int) RSVDt->RSVD.k);CHKERRQ(ierr);
			RSVDt->RSVD.k = 3;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting 'k' to default value: %d\n", (int) RSVDt->RSVD.k);CHKERRQ(ierr);
		}
		ierr = PetscOptionsGetInt(NULL,NULL,"-q",&RSVDt->RSVD.q,&flg_set);CHKERRQ(ierr);
		if (!flg_set) {
			RSVDt->RSVD.q = 0;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'q' variable not found. Setting 'q' to default value: %d\n", (int) RSVDt->RSVD.q);
		} else if (RSVDt->RSVD.q < 0) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'q' must be a non-negative integer, current value: %d\n", (int) RSVDt->RSVD.q);CHKERRQ(ierr);
			RSVDt->RSVD.q = 0;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting 'q' to default value: %d\n", (int) RSVDt->RSVD.q);CHKERRQ(ierr);
		}
		if (RSVDt->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n********************************************\n"
				"*************** Problem info ***************\n********************************************\n\n");CHKERRQ(ierr);		
	} else { // Transient run variables
		ierr = PetscOptionsGetReal(NULL,NULL,"-TransDivVal",&TR->TransDivVal,&flg_set);CHKERRQ(ierr);
		if (!flg_set) {
			TR->TransDivVal = 1e3;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransDivVal' variable not found. Setting 'TransDivVal' to default value: %g\n", TR->TransDivVal);
		} else if (TR->TransDivVal <= 0) {
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Warning: 'TransDivVal' must be positive, current value: %g\n", TR->TransDivVal);CHKERRQ(ierr);
			TR->TransDivVal = 1e3;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting 'TransDivVal' to default value: %g\n", TR->TransDivVal);CHKERRQ(ierr);
		}	
		ierr = PetscOptionsGetReal(NULL,NULL,"-TransConVal",&TR->TransConVal,&flg_set);CHKERRQ(ierr);
		if (!flg_set) {
			TR->TransConVal = 1e-8;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransConVal' variable not found. Setting 'TransConVal' to default value: %g\n", TR->TransConVal);
		} else if (TR->TransConVal <= 0) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransConVal' must be positive, current value: %g\n", TR->TransConVal);CHKERRQ(ierr);
			TR->TransConVal = 1e-8;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting 'TransConVal' to default value: %g\n", TR->TransConVal);CHKERRQ(ierr);
		}
		ierr = PetscOptionsGetBool(NULL,NULL,"-TransSave",&TR->TransSave,&flg_set);CHKERRQ(ierr);
		if (!flg_set) {
			TR->TransSave = 0;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransSave' variable not found. Setting 'TransSave' to default value: %d\n", (int) TR->TransSave);
		}
		ierr = PetscOptionsGetBool(NULL,NULL,"-TransRemovalEst",&TR->TransRemovalEst,&flg_set);CHKERRQ(ierr);
		if (!flg_set) {
			TR->TransRemovalEst = 0;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransRemovalEst' variable not found. Setting 'TransRemovalEst' to default value: %d\n", (int) TR->TransRemovalEst);
		}
		ierr = PetscOptionsGetInt(NULL,NULL,"-TransSaveMod",&TR->TransSaveMod,&flg_set);CHKERRQ(ierr);
		if (!flg_set) {
			TR->TransSaveMod = 100;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransSaveMod' variable not found. Setting 'TransSaveMod' to default value: %d\n", (int) TR->TransSaveMod);
		} else if (TR->TransSaveMod < 1) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransSaveMod' must be positive, current value: %d\n", (int)TR->TransSaveMod);CHKERRQ(ierr);
			TR->TransSaveMod = 100;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting 'TransSaveMod' to default value: %d\n", (int)TR->TransSaveMod);CHKERRQ(ierr);
		}
		ierr = PetscOptionsGetInt(NULL,NULL,"-TransPeriods",&TR->TransPeriods,&flg_set);CHKERRQ(ierr);
		if (!flg_set) {
			TR->TransPeriods = 1;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransPeriods' variable not found. Setting 'TransPeriods' to default value: %d\n", (int) TR->TransPeriods);
		}  else if (TR->TransPeriods < 1) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransPeriods' must be positive, current value: %d\n", (int)TR->TransPeriods);CHKERRQ(ierr);
			TR->TransPeriods = 1;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting 'TransPeriods' to default value: %d\n", (int)TR->TransPeriods);CHKERRQ(ierr);
		}
		ierr = PetscOptionsGetBool(NULL,NULL,"-TransICFlg",&TR->TransICFlg,&flg_set);CHKERRQ(ierr);
		if (!flg_set) {
			TR->TransICFlg = 0;
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransICFlg' variable not found. Setting 'TransICFlg' to default value: %d\n", (int) TR->TransICFlg);
		} else if (TR->TransICFlg) {
			ierr = PetscOptionsGetString(NULL,NULL,"-TransICDir",(char*)&dirs->TransICDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
			if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Transient IC flag is on! Either set 'TransICFlg' to zero or specify 'TransICDir'");CHKERRQ(ierr);
		}
	}

	PetscFunctionReturn(0);

}

