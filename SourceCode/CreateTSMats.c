
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode CreateTSMats(TS_matrices *TS, PetscInt N)
{
	/*
		Creates the required RK4 time-stepping matrices 
	*/  

	PetscErrorCode        ierr;

	PetscFunctionBeginUser;

	ierr = VecCreate(PETSC_COMM_WORLD,&TS->F1);CHKERRQ(ierr);
	ierr = VecSetType(TS->F1,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS->F1,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS->F2);CHKERRQ(ierr);
	ierr = VecSetType(TS->F2,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS->F2,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS->F3);CHKERRQ(ierr);
	ierr = VecSetType(TS->F3,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS->F3,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS->k1);CHKERRQ(ierr);
	ierr = VecSetType(TS->k1,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS->k1,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS->k2);CHKERRQ(ierr);
	ierr = VecSetType(TS->k2,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS->k2,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS->k3);CHKERRQ(ierr);
	ierr = VecSetType(TS->k3,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS->k3,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS->k4);CHKERRQ(ierr);
	ierr = VecSetType(TS->k4,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS->k4,PETSC_DECIDE,N);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&TS->y_temp);CHKERRQ(ierr);
	ierr = VecSetType(TS->y_temp,VECMPI);CHKERRQ(ierr);
	ierr = VecSetSizes(TS->y_temp,PETSC_DECIDE,N);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}



