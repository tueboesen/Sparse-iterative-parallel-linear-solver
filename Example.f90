    !MIT License
    !
    !Copyright (c) [2017] [HydroGeophysics Group, Aarhus University, Denmark]
    !
    !Permission is hereby granted, free of charge, to any person obtaining a copy
    !of this software and associated documentation files (the "Software"), to deal
    !in the Software without restriction, including without limitation the rights
    !to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    !copies of the Software, and to permit persons to whom the Software is
    !furnished to do so, subject to the following conditions:
    !
    !The above copyright notice and this permission notice shall be included in all
    !copies or substantial portions of the Software.
    !
    !THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    !IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    !FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    !AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    !LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    !OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    !SOFTWARE.

    !This is a Simple example that shows the sparse Iterative Solver from Aarhus Hydrogeogphysical Group in action.
    !This example loads a linear system from txt files in CSR format, it loads an ini-file with various settings, and then it solves the linear system
    program Solver_example
    use mParallel
    use mSolver
    use mSMat
    use mInput
    implicit none
    integer :: AffinityMode,NCPUs,UseNested,StartThread,NCPUsLow,NCPUInner,NCPUOuter,NCPUInnerLow,AffinityPattern
    integer :: n,i,nCols,nRows,UseRCM
    real*8, pointer :: rhs(:),x(:),rhs_backup(:),Residual(:)
    character*256 :: SettingFileName,SparseMatrixFiles,OutFile
    CHARACTER(len=256) :: arg
    type(TSparseMat) :: A
    type(Tsettings)  :: Set
    !Solver vars
    integer :: MaxIter, NmodPar, NoRhs, Solverinformation, solvertype,UseNormalization
    real*8  :: FillingFactor,error,t1,t2
    type(TSolverSettings) :: SolverSettings

    !First we start by loading the additional arguments, which should have been provided with the call. 
    CALL get_command_argument(1, arg) !Argument 1 is the name of the settingsfile
    IF (LEN_TRIM(arg) == 0) then
      print*,'ERROR not enough arguments given'  
    end if
    read(arg,*) SettingFileName
  
    CALL get_command_argument(2, arg) !Argument 2 is the name of the SparseMatrix files.
    IF (LEN_TRIM(arg) == 0) then
      print*,'ERROR not enough arguments given'  
    end if
    read(arg,*) SparseMatrixFiles

    CALL get_command_argument(3, arg) !Argument 3 is the name of the output file.
    IF (LEN_TRIM(arg) == 0) then
      print*,'ERROR not enough arguments given'  
    end if
    read(arg,*) OutFile
    
    print*,'SettingFileName:',SettingFileName
    print*,'SparseMatrixFiles:',SparseMatrixFiles
    print*,'OutFile',OutFile
    call LoadSettings(SettingFileName,set)
    
    !Before we do anything else openMp related, we call initOpenMp()
    call InitOpenMP(set%AffinityMode,set%AffinityPattern,set%StartThread,set%UseNested,set%NCPUs,set%NCPUsLow,set%NCPUOuter)

    !Now lets create a sparse linear system to work on
    call SMatReadSparse(A,rhs, SparseMatrixFiles,Set%BlockSize) !Read the sparsematrix files into memory
    allocate(x(A%norows))
    allocate(rhs_backup(A%noRows))
    allocate(Residual(A%noRows))
    rhs_backup=rhs
    UseRCM=0
    x=0
    t1=omp_get_wtime() !Start a timer
    call SetSolverSettings(SolverSettings,A,UseRCM,set%UseNormalization,set%SolverType,set%MaxIter,set%blocksize,set%FillingFactor) !Set the solver and create the preconditioner factorization
    call SingleRHS_wrapper(A,rhs,x,SolverSettings) !Propagate the actual linear system
    t2=omp_get_wtime() 
    call SolverNormalizeMatrix(solversettings,A,.false.) !Renormalize matrix A (This is only needed because we want to validate that the solution did indeed solve the initial matrix).
    
    call SMatAmultV(A,x,rhs) !Show that the solution we found is indeed a solution to the linear system by multiplying the solution on the matrix
    Residual=rhs-rhs_backup
    Error=sqrt(ddot(A%NoRows,Residual,1,Residual,1))
    print*,'Error',Error
    print*,'Time for solve',t2-t1
    open(78,file=trim(OutFile),position='append')
    write(78,*) t2-t1
    write(78,*) solversettings%output%status
    write(78,*) Error
    write(78,*) omp_get_max_threads()
    write(78,*) set%SolverType
    close(78)
    
    end program Solver_example

