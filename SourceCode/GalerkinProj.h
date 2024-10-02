
#ifndef GALERKINPROJ_H
#define GALERKINPROJ_H

PetscErrorCode GalerkinProj(Mat, Mat, Mat, PetscInt, RSVDt_vars*, TS_removal_matrices*, DFT_matrices*);
#endif