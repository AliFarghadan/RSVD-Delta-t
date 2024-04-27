
#ifndef REMOVETRANSIENTPROJECTION_H
#define REMOVETRANSIENTPROJECTION_H

PetscErrorCode RemoveTransientProjection(Mat, Mat, Mat, PetscInt, RSVDt_vars*, TS_removal_matrices*, DFT_matrices*);
#endif