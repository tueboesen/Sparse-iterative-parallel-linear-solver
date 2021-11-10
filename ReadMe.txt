------------------------------------
------------Introduction------------
------------------------------------
This is a sparse iterative linear solver developed at the HydroGeophysics Group, Aarhus University, Denmark.
The code is developed in Fortran/OpenMP, and is designed with focus on speed especially in parallel.
The code uses a block-splitting of the linear system without any overlap. 
It is intended to be used on matrices with have already been reordered using a reverse Cuthill-Mckee reordering.


------------------------------------
-------What is in this pack---------
------------------------------------
Example.f90              !Contains a small example, that initialize the parallelization, loads a matrix and solves the linear system
Solver.f90               !Contains the iterative sparse solver 
SparseMatrix.f90         !Contains various sparse matrix functions needed by the iterative solver 
Paralllelization.f90     !Contains the parallelization framework
Misc.f90                 !Contains various small helping routines needed by the other modules   
Load.f90                 !Contains a minimalistic structure for loading and setting the various settings for the example
Matrix.zip               !Contains a matrix saved in compressed sparse row format
Readme.txt               !Contains information needed for running the example
Technical_Report.pdf     !A technical report which introduces parallelization and explains why and how the parallelization framework was created. Furthermore it contains information about the sparse iterative solver also created by the Hydrogeophysics Group at Aarhus University.
Settings.ini             !Settings for the parallelization framework and the iterative solver
------------------------------------
---------How to set it up-----------
------------------------------------
Open either the parallelization solution file or project file, or built your own project by using Example.f90 as the main program and the rest as dependencies.

once the project have been built, place all files in the same folder and run it like this:

-> solver.exe settings.ini 001 outputfilename        (Where 001 is the name of the sparseMatrix in the matrix folder)

------------------------------------
-------------License----------------
------------------------------------
This code is provided freely under the MIT License.

Copyright (c) [2017] [HydroGeophysics Group, Aarhus University, Denmark]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.