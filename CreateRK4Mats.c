
#include <CreateRK4Mats.h>

PetscErrorCode CreateRK4Mats(TS_matrices *TS_mat, PetscInt N)
{

	PetscErrorCode      ierr;

	PetscFunctionBeginUser;

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->F1);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->F1,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->F1,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->F2);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->F2,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->F2,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->F3);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->F3,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->F3,PETSC_DECIDE,N);CHKERRQ(ierr); 

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->k1);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->k1,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->k1,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->k2);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->k2,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->k2,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->k3);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->k3,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->k3,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->k4);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->k4,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->k4,PETSC_DECIDE,N);CHKERRQ(ierr);	

	ierr = VecCreate(PETSC_COMM_WORLD,&TS_mat->y_temp);CHKERRQ(ierr);
	ierr = VecSetType(TS_mat->y_temp,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS_mat->y_temp,PETSC_DECIDE,N);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

