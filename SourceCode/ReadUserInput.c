
#include <petscmat.h>
#include "Variables.h"

PetscErrorCode ReadUserInput(RSVDt_vars *RSVDt, WS_matrices *WS_mat, LNS_vars *LNS_mat, \
							Resolvent_matrices *Res_mat, TransientRun_vars *TR_vars, Directories *dirs)
{
	/*
		Reads in user input parameters
	*/

	PetscErrorCode        ierr;
	PetscBool             flg_set;

	PetscFunctionBeginUser;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** Read in user inputs ***\n");CHKERRQ(ierr);

	ierr = PetscOptionsGetReal(NULL,NULL,"-TransLen",&RSVDt->TS.TransLen,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVDt->TS.TransLen = 100;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransLen' variable not found. Setting 'TransLen' to default value: %g\n", RSVDt->TS.TransLen);
	}
	ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&RSVDt->TS.dt,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'dt'");CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,NULL,"-w",&RSVDt->RSVD.w,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'w'");CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-Twopi_flg",&RSVDt->RSVD.twopi_flg,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVDt->RSVD.twopi_flg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'Twopi_flg' variable not found. Setting 'Twopi_flg' to default value: %d\n", (int) RSVDt->RSVD.twopi_flg);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-TransRun",&TR_vars->TransientRun,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		TR_vars->TransientRun = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransRun' variable not found. Setting 'TransRun' to default value: %d\n", (int) TR_vars->TransientRun);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-real_A",&RSVDt->TS.flg_real_A,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVDt->TS.flg_real_A = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'real_A' variable not found. Setting 'real_A' to default value: %d\n", (int) RSVDt->TS.flg_real_A);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-TransRemoval",&RSVDt->TS.TransRemoval,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVDt->TS.TransRemoval = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransRemoval' variable not found. Setting 'TransRemoval' to default value: %d\n", (int) RSVDt->TS.TransRemoval);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-Discounting",&LNS_mat->RSVDt.disc.flg_disc,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		LNS_mat->RSVDt.disc.flg_disc = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'Discounting' variable not found. Setting 'Discounting' to default value: %d\n", (int) LNS_mat->RSVDt.disc.flg_disc);
	}
	ierr = PetscOptionsGetString(NULL, NULL,"-root_dir",(char*)&dirs->root_dir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL,"-input_folder",(char*)&dirs->input,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL,"-results_dir",(char*)&dirs->output,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL,"-Mat_filename",(char*)&dirs->Mat_filename,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
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
	ierr = PetscOptionsGetInt(NULL,NULL,"-display",&RSVDt->display,&flg_set);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,NULL,"-Beta",&LNS_mat->RSVDt.disc.Beta,&flg_set);CHKERRQ(ierr);
	if (!flg_set && LNS_mat->RSVDt.disc.flg_disc) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Discounting flag is on! Either set Discounting to zero or specify 'Beta'");CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,NULL,"-DivVal",&TR_vars->DivVal,&flg_set);CHKERRQ(ierr);
	if (!flg_set && TR_vars->TransientRun) {
		TR_vars->DivVal = 1e3;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'DivVal' variable not found. Setting 'DivVal' to default value: %g\n", TR_vars->DivVal);
	}
	ierr = PetscOptionsGetReal(NULL,NULL,"-ConVal",&TR_vars->ConVal,&flg_set);CHKERRQ(ierr);
	if (!flg_set && TR_vars->TransientRun) {
		TR_vars->ConVal = 1e-8;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'ConVal' variable not found. Setting 'ConVal' to default value: %g\n", TR_vars->ConVal);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-TransientSave",&TR_vars->TransientSave,&flg_set);CHKERRQ(ierr);
	if (!flg_set && TR_vars->TransientRun) {
		TR_vars->TransientSave = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TransientSave' variable not found. Setting 'TransientSave' to default value: %d\n", (int) TR_vars->TransientSave);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-EfficientTransientTest",&TR_vars->EffTransTest,&flg_set);CHKERRQ(ierr);
	if (!flg_set && TR_vars->TransientRun) {
		TR_vars->EffTransTest = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'EfficientTransientTest' variable not found. Setting 'EfficientTransientTest' to default value: %d\n", (int) TR_vars->EffTransTest);
	}
	ierr = PetscOptionsGetInt(NULL,NULL,"-Rseed",&TR_vars->Rseed,&flg_set);CHKERRQ(ierr);
	if (!flg_set && TR_vars->TransientRun) {
		TR_vars->Rseed = 1373;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'Rseed' variable not found. Setting 'Rseed' to default value: %d\n", (int) TR_vars->Rseed);
	}
	ierr = PetscOptionsGetInt(NULL,NULL,"-mod",&TR_vars->mod,&flg_set);CHKERRQ(ierr);
	if (!flg_set && TR_vars->TransientRun) {
		TR_vars->mod = 10;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'mod' variable not found. Setting 'mod' to default value: %d\n", (int) TR_vars->mod);
	}
	ierr = PetscOptionsGetInt(NULL,NULL,"-Periods",&TR_vars->Periods,&flg_set);CHKERRQ(ierr);
	if (!flg_set && TR_vars->TransientRun) {
		TR_vars->Periods = 1;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'Periods' variable not found. Setting 'Periods' to default value: %d\n", (int) TR_vars->Periods);
	}   
	RSVDt->RSVD.Nw_eff = RSVDt->TS.flg_real_A ? RSVDt->RSVD.Nw/2 : RSVDt->RSVD.Nw;

	PetscFunctionReturn(0);

}

// PetscErrorCode ReadUserInput(RSVDt_vars *RSVDt, WS_matrices *WS_mat, LNS_vars *LNS_mat, \
//                              Resolvent_matrices *Res_mat, TransientRun_vars *TR_vars, Directories *dirs)
// {

//  /*
//      Reads in user input parameters
//  */

//  PetscErrorCode        ierr;
//  char                  filename[PETSC_MAX_PATH_LEN];
//  PetscBool             input_flg;
//  FILE                 *file;
//  char                  line[PETSC_MAX_PATH_LEN];
//  char                  option[PETSC_MAX_PATH_LEN];
//  char                  value_str[PETSC_MAX_PATH_LEN];

//  PetscFunctionBeginUser;

//  ierr = PetscOptionsGetString(NULL, NULL,"-inputs",(char*)&filename,PETSC_MAX_PATH_LEN,&input_flg);CHKERRQ(ierr);
//  if (!input_flg) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify inputs variables via -inputs");CHKERRQ(ierr);
	
//  file = fopen(filename, "r");
//  if (file == NULL) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Wrong input filename");CHKERRQ(ierr);

//  while (fgets(line, PETSC_MAX_PATH_LEN, file) != NULL) {
//         if (line[0] == '#' || line[0] == '/') {
//             continue;
//         } else if (sscanf(line, "%s = %s", option, value_str) == 2) {
//          if (strcmp(option, "k") == 0) RSVDt->RSVD.k = atoi(value_str);
//          else if (strcmp(option, "q") == 0) RSVDt->RSVD.q = atoi(value_str);
//          else if (strcmp(option, "Nw") == 0) RSVDt->RSVD.Nw = atoi(value_str);
//          else if (strcmp(option, "display") == 0) RSVDt->display = atoi(value_str);
//          else if (strcmp(option, "diagW") == 0) WS_mat->flg_diag = atoi(value_str);
//          else if (strcmp(option, "Discounting") == 0) LNS_mat->RSVDt.disc.flg_disc = atoi(value_str);
//          else if (strcmp(option, "Twopi_flg") == 0) RSVDt->RSVD.twopi_flg = atoi(value_str);
//          else if (strcmp(option, "TransRun") == 0) TR_vars->TransientRun = atoi(value_str);
//          else if (strcmp(option, "real_A") == 0) RSVDt->TS.flg_real_A = atoi(value_str);
//          else if (strcmp(option, "TransRemoval") == 0) RSVDt->TS.flg_TransRemovals = atoi(value_str);
//          else if (strcmp(option, "all_in_one") == 0) Res_mat->all_in_one_flg = atoi(value_str);
//          else if (strcmp(option, "dt_max") == 0) RSVDt->TS.dt_max = atof(value_str);
//          else if (strcmp(option, "w") == 0) RSVDt->RSVD.w = atof(value_str);
//          else if (strcmp(option, "Beta") == 0) LNS_mat->RSVDt.disc.Beta = atof(value_str);
//          else if (strcmp(option, "root_dir") == 0) strcpy(dirs->root_dir, value_str);
//          else if (strcmp(option, "input_folder") == 0) strcpy(dirs->input, value_str);
//          else if (strcmp(option, "output_folder") == 0) strcpy(dirs->output, value_str);
//          else if (strcmp(option, "Mat_filename") == 0) strcpy(dirs->Mat_filename, value_str);
//          else if (strcmp(option, "TransientSave") == 0) TR_vars->TransientSave = atoi(value_str);
//          else if (strcmp(option, "TransientRun_dir") == 0) TR_vars->TransientRun_dir = atoi(value_str);
//          else if (strcmp(option, "EfficientTransientTest") == 0) TR_vars->EffTransTest_flg = atoi(value_str);
//          else if (strcmp(option, "k_trans") == 0) TR_vars->k_trans = atoi(value_str);
//          else if (strcmp(option, "mod") == 0) TR_vars->mod = atoi(value_str);
//          else if (strcmp(option, "Rseed") == 0) TR_vars->Rseed = atoi(value_str);
//          else if (strcmp(option, "Periods") == 0) TR_vars->Periods = atoi(value_str);
//      }
//  }

//  fclose(file);

//  RSVDt->RSVD.Nw_eff = RSVDt->TS.flg_real_A ? RSVDt->RSVD.Nw/2 : RSVDt->RSVD.Nw;

//  PetscFunctionReturn(0);

// }

