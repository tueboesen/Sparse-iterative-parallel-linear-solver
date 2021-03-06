# Sparse Iterative parallel linear solver

## Introduction

This is a sparse iterative linear solver I developed while doing my Phd at the HydroGeophysics Group, Aarhus University, Denmark.
The code is developed in Fortran/OpenMP, and is designed with focus on speed especially in parallel.
The code uses a block-splitting of the linear system without any overlap. 
It is intended to be used on matrices which have already been reordered using a reverse Cuthill-Mckee reordering.


## What is in this pack

| Filename    | Description |
| ----------- | ----------- |
| Example.f90 | Contains a small example, that initialize the parallelization, loads a matrix and solves the linear system.       |
| Solver.f90  | Contains the iterative sparse solver. |
| SparseMatrix.f90  | Contains various sparse matrix functions needed by the iterative solver. |
| Paralllelization.f90  | Contains the parallelization framework. |
| Misc.f90  | Contains a few helper routines needed by the other modules. |
| Load.f90  | Contains a minimalistic structure for loading and setting the various settings for the example. |
| Matrix folder  | Contains the matrix used in the example saved in compressed sparse row format.|
| Technical_Report.pdf  | A technical report which introduces both the parallelization framework and the iterative solver.|
| Settings.ini  | Settings for the parallelization framework and the iterative solver.|

## How to set it up

Open either the parallelization solution file or project file, or build your own project by using Example.f90 as the main program and the rest as dependencies.

Once the project has been built, place all files in the same folder and run it like this:

-> solver.exe settings.ini 001 outputfilename        (where 001 is the name of the sparseMatrix in the matrix folder)
