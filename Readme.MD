# Sparse Iterative parallel linear solver

## Introduction

This is a sparse iterative linear solver developed at the HydroGeophysics Group, Aarhus University, Denmark.
The code is developed in Fortran/OpenMP, and is designed with focus on speed especially in parallel.
The code uses a block-splitting of the linear system without any overlap. 
It is intended to be used on matrices with have already been reordered using a reverse Cuthill-Mckee reordering.


## What is in this pack

Example.f90              !Contains a small example, that initialize the parallelization, loads a matrix and solves the linear system

Solver.f90               !Contains the iterative sparse solver 

SparseMatrix.f90         !Contains various sparse matrix functions needed by the iterative solver 

Paralllelization.f90     !Contains the parallelization framework

Misc.f90                 !Contains various small helping routines needed by the other modules   

Load.f90                 !Contains a minimalistic structure for loading and setting the various settings for the example

Matrix folder            !Contains a matrix saved in compressed sparse row format

Technical_Report.pdf     !A technical report which introduces parallelization and explains why and how the parallelization framework was created. Furthermore it contains information about the sparse iterative solver also created by the Hydrogeophysics Group at Aarhus University.

Settings.ini             !Settings for the parallelization framework and the iterative solver

## How to set it up

Open either the parallelization solution file or project file, or built your own project by using Example.f90 as the main program and the rest as dependencies.

once the project have been built, place all files in the same folder and run it like this:

-> solver.exe settings.ini 001 outputfilename        (Where 001 is the name of the sparseMatrix in the matrix folder)
