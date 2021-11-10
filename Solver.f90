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
    
!Tue jan, 2015
!This module contains iterative solvers for linear systems

module mSolver
  use mMisc
  use mSMat
  implicit none 

  type TSolverStatus
    !This type is a subtype in TSolverSettings, meant for output from the solve
    logical, private :: Success=.false.    !True if system converged
    integer          :: Status             !Status is the variable that gets reported to the surrounding code, on a succesfull solve, status=nint(iter) but otherwise it can contain error messages or -1 if no convergence was reached
    real*8, private  :: AbsError           !Average absolute error
    real*8, private  :: RelError           !Average relative error
    real*8, private  :: FactorizationTime  !Time spent factorizing
    real*8, private  :: SolverTime         !Time spent solving
    
    !MultiRHS variables    
    real*8, private   :: Iter              !Average number of iterations used    
    integer, private  :: Iter_min          !Minimum number of iterations
    integer, private  :: Iter_max          !Maximum number of iterations
    real*8, private   :: AbsError_min      !Minimum absolute error 
    real*8, private   :: AbsError_max      !Maximum absolute error
    real*8, private   :: RelError_min      !Minimum relative error
    real*8, private   :: RelError_max      !Maximum relative error
  end type TSolverStatus
 
  type TSolverSettings
    !This type is meant to contain all settings needed for the solver to run, in addition it will also contain the output from the solve.
    logical, private :: IsPrepared=.false. !Before solver runs it checks this to see whether solversettings have been prepared.
    logical, private :: IsComplex          !True if system is complex
    integer, private :: UseNormalization   !If True, the matrix is normalized as: D^-1*A*D^-1 (D*x) = D^-1*b, where D=sqrt(diag(A))
    integer, private :: UseRCM             !If True, the matrix is reordered using Reverse Cuthill-Mckee reordering    
    integer, private :: NoRHS              !Number of right hands sides
    integer, private :: PrecondMaxFill     !Max number of elements each row in LUT factorization can contain (excluding the diagonal)
    integer, private :: NoElements         !Actual number of elements used in LUT factorization (not sure this is set)
    integer, private :: MaxIter            !Max number of Iterations
    integer, private :: SolverType         !Which solver to use
    integer, private :: Inform             !Set how much information solver should print, 0 = silent, 1 = basic statistics, 2 = solversettings+basic statistics, 3=full information
    integer, private :: Blocksize          !The minimum blocksize 
    integer, private :: NBlocks            !The number of blocks used in the parallel preconditioner computation
    integer, private :: FillingFactor      !The filling factor used for computing PrecondMaxFill
    real*8, private  :: PrecondTolerance   !DropTolerance of ILUT filling
    real*8, private  :: StepLimit          !When to stop if convergence is not reached
    real*8, private  :: AbsTol             !Convergence criteria 
    real*8, private  :: RelTol             !Convergence criteria 
    real*8, private  :: DiagDom            !Diagonal dominance, an empirical parameter that gives a rough measure for how difficult this system will be to solve, and hence how many elements we need in LUT factorization
    type(TSolverStatus) :: Output

    !RCM ordering indices
    integer,dimension(:),allocatable,private    :: RCMIndices   !Indices returned from RCM reordering, in case it should be reversed at the end
    !Normalization factor
    real*8,dimension(:),allocatable,private    :: Normalization   !Normalization is used if the logical Variable UseNormalization is true
    real*8,dimension(:),allocatable,private    :: ReNormalization !Renormalization is used if the logical variable UseNormalization is true, this brings back the original Matrix afterwards.
    !Block splitting
    integer,dimension(:,:),allocatable, private :: Blocks         !The actual blocks the matrix will be split up in during a parallel run
    !LU0/LUT preconditioner
    integer,dimension(:),allocatable, private   :: PrecondCols    !Used for LUT precondition to store the ColIdx in CSR format    
    integer,dimension(:),allocatable, private   :: PrecondRows     !Used for LUT precondition to store the RowIdx in CSR format
    real*8,dimension(:),allocatable, private    :: PrecondVals     !Used in real LU0/LUT preconditioning to store the values.
    !LU0 - complex values
    complex*16,dimension(:),allocatable, private  :: PrecondcVals  !Used in complex LU0 preconditioning to store the values
    !BLUT preconditioner
    integer,dimension(:,:),allocatable, private   :: BPrecondIdxs  !Used for BLUT precondition to store the ColIdx in block CSR format  
    integer,dimension(:,:),allocatable, private   :: BPrecondDiag  !Used for BLUT precondition to store the RowIdx in block CSR format
    real*8,dimension(:,:),allocatable, private    :: BPrecondVals  !Used in BLUT preconditioning to store the values of each block.
    !settings for PARDISO solver
    
  end type TSolverSettings

  interface SolverAxB
    module procedure SingleRHS_wrapper   
  end interface    
  
  contains


!****************************************************************************************    
   
      
    subroutine SingleRHS_wrapper(iA,iRHS,ioX,SolverSettings)
    !Tue, Dec 2014
    !This routine is a wrapper routine for all the our solvers which handle systems with a single righthandside.
    !IO
    !iA -> input matrix, in CSR format, and of the type TSparseMat
    !iRHS                    -> a system of right hand sides, RHS(:,1) is one right hand side vector
    !ioX                     -> on input this matrix contains the solution guess for each vector, on output it contains the found solutions
    !SolverSettings          -> contains all the relevant settings the solvers need, it is also where stats about the solve is saved. 
    !                           This needs to be set before calling this routine.
    !
    use omp_lib
    use msMat
    implicit none
    type (TSparseMat),intent(inout)    :: iA
    type(TSolverSettings),intent(inout):: SolverSettings 
    real*8,intent(inout),dimension(:)  :: iRHS
    real*8,intent(inout),dimension(:)  :: ioX
    character*200                      :: S
    integer                            :: i,j
    integer                            :: NoNumaNodes,Threads,ostatus,error
    real*8                             :: RunTimeBegin,RunTimeEnd
    real*8,dimension(:),allocatable    :: InvDiag,tmp
    integer,parameter                  :: NoRowsDenseLimit=35000 !with 35000 elements the dense solver takes about 150 seconds to solve the problem
    integer                            :: OldThreads
    logical                            :: dummy
    logical(4)                        :: Success
    allocate(tmp(iA%NoRows))
    dummy = startOpenMP(2,1)
    
    RuntimeBegin =  omp_get_wtime()
    Call checkSolverSettings(SolverSettings)
 
    if(SolverSettings%UseRCM) then
      tmp(:)=iRHS(SolverSettings%RCMIndices)
      iRHS=tmp
      tmp(:)=ioX(SolverSettings%RCMIndices)
      ioX=tmp
    end if
    if(SolverSettings%UseNormalization) then
      do i=1,iA%NoRows
        iRHS(i)=iRHS(i)*SolverSettings%Normalization(i)
        ioX(i)=ioX(i)*SolverSettings%ReNormalization(i) 
      end do       
    end if
    
    SELECT CASE (SolverSettings%SolverType)
    CASE(1,2)
      call SingleRHS(iA,iRHS,ioX,SolverSettings)
    CASE(3)
     !Pardiso (Not implemented here)   
    END SELECT
    
    !If we used Normalization, we now need to renormalize our solution vector x and our input iRHS.
    if(SolverSettings%UseNormalization) then
      do i=1,iA%NoRows
        ioX(i)=ioX(i)*SolverSettings%Normalization(i)
        iRHS(i)=iRHS(i)*SolverSettings%ReNormalization(i)
      end do
    end if
    if(SolverSettings%UseRCM) then !This is not right, it needs to be the inverse here, but should be easily fixed
      do i=1,iA%NoRows  
        tmp(SolverSettings%RCMIndices(i)) = iRHS(i)
      end do      
      iRHS=tmp
      do i=1,iA%NoRows  
        tmp(SolverSettings%RCMIndices(i)) = ioX(i)
      end do      
      ioX=tmp
    end if

    RuntimeEnd=omp_get_wtime()
    SolverSettings%Output%SolverTime=RunTimeEnd-RunTimeBegin
    if (.not.SolverSettings%Output%success) SolverSettings%Output%Status = -SolverSettings%Output%Status
    call ReturnOldThreadNum()
    return
    end subroutine SingleRHS_wrapper
     
!****************************************************************************************    
    

    
    subroutine SingleRHS(iA,iRHS,ioX,SolverSettings)
    !Tue, November 2014
    !Block Preconditioned Conjugate Gradient, with uma.
    !This routine solves linear systems using the Preconditioned Conjugate gradient method, as a preconditioner the method uses either SOR or BLUT depending on SolverSettings 
    !More information on the method can be found in the book "Iterative methods for sparse linear systems".
    !IO
    !iA -> input matrix, in CSR format, and of the type TSparseMat
    !iRHS                    -> a system of right hand sides, RHS(:,1) is one right hand side vector
    !ioX                     -> on input this matrix contains the solution guess for each vector, on output it contains the found solutions
    !SolverSettings          -> contains all the relevant settings the solvers need, it is also where stats about the solve is saved. 
    !                           This needs to be set before calling this routine.
    !
    !
    use omp_lib
    implicit none

    type (TSparseMat),intent(inout)   :: iA
    type(TSolverSettings),intent(inout)  :: SolverSettings 

    real*8,intent(in),dimension(:)       :: iRHS
    real*8,intent(inout),dimension(:)    :: ioX
    
    !Parallel vars
    Integer             :: Threads
    !Local vars
    integer            :: i
    integer,save       :: j
    integer            :: MaxIter,NoRows,iNuma
    real*8             :: rzDotProd,eps
    real*8             :: alfa, beta
    real*8             :: error,norm
    real*8,allocatable :: r(:),z(:),p(:),Ap(:),OldX(:),StartX(:)
    logical            :: KeepLooping 

    
    !SOR vars    
    integer,dimension(:),allocatable           :: DiagIdx
    integer,dimension(:,:),allocatable         :: BRowIdx
    real*8,dimension(:),allocatable            :: InvDiag

        
    Eps=1e-6
    iNuma=1
    NoRows=size(iRHS)
    allocate(r(NoRows))
    allocate(z(NoRows))
    allocate(p(NoRows))
    allocate(Ap(NoRows))
    allocate(OldX(NoRows))
    allocate(StartX(NoRows))
    OldX=ioX         
    StartX=ioX
    !First the block parallisation
    Threads=OMP_GET_MAX_THREADS ()
    
    call SmatAmultV(iA,ioX(:),r(:))
    call VectorAddition(1,iNuma,NoRows,r,iRHS) !r(:)=iRHS(:)-r(:)
    
    SELECT CASE (SolverSettings%SolverType)
    CASE(1)
      allocate(DiagIdx(iA%NoRows))  
      call SMatGetDiagIdx(iA,DiagIdx) 
      allocate(InvDiag(iA%NoRows))  
      call SMatGetInvDiag(iA,InvDiag)
      allocate(BRowIdx(iA%NoRows,2))
      call SMatGetBlockRowIdx(iA,SolverSettings%Blocks,BRowIdx)        
      call Apply_BPrecondition_BSGS(iA%ColIdxs,iA%RowIdxs,iA%Vals,BRowIdx,DiagIdx,InvDiag,r,z,SolverSettings%Blocks)
    CASE(2)
      call Apply_BPrecondition_BLUT(SolverSettings%Blocks, r, z, SolverSettings%BPrecondVals, SolverSettings%BPrecondIdxs, SolverSettings%BPrecondDiag)
    END SELECT
    
    call VectorAddition(0,iNuma,NoRows,p,z) !p(:)=z(:)
    
    do j = 1 , SolverSettings%MaxIter
      call SMatAMultV(iA,p,Ap)
    
      rzDotProd=ddot(NoRows,r,1,z,1)
      alfa=rzDotProd/ddot(NoRows,Ap(:),1,p,1)
      call VectorAddition(4,iNuma,NoRows,ioX,p,alfa,r,Ap) !ioX(:) = ioX(:) + alfa*p(:) , r(:) = r(:) - alfa*Ap(:)
    
      SELECT CASE (SolverSettings%SolverType)
      CASE(1)
       call Apply_BPrecondition_BSGS(iA%ColIdxs,iA%RowIdxs,iA%Vals,BRowIdx,DiagIdx,InvDiag,r,z,SolverSettings%Blocks)
      CASE(2)
        call Apply_BPrecondition_BLUT(SolverSettings%Blocks, r, z, SolverSettings%BPrecondVals, SolverSettings%BPrecondIdxs, SolverSettings%BPrecondDiag)
      END SELECT
      beta = ddot(NoRows,r,1,z,1)/rzDotProd
      call VectorAddition(3,iNuma,NoRows,p,z,beta) !p(:) = z(:) + beta * p(:)
      
      !Did we converge yet?
      error=0d0
      norm =0d0
      Error=sqrt(ddot(NoRows,r,1,r,1))
      Norm=sqrt(ddot(NoRows,ioX,1,ioX,1))
            
      if (((error/norm).le.SolverSettings%RelTol).or.(error).le.SolverSettings%AbsTol) then
        SolverSettings%Output%Success=.TRUE.
        SolverSettings%Output%Status = min(j,SolverSettings%MaxIter)
        exit
      end if
      if (SolverSettings%StepLimit.gt.0) then
          Keeplooping=.False.
          do i=1,NoRows
            if (abs(ioX(i)-OldX(i))/Norm.gt.SolverSettings%StepLimit) then
              Keeplooping=.True.
              exit
            end if
          end do
          if (.NOT.keeplooping) then 
            exit
          else
            OldX(:) = ioX(:)
          end if
      end if
    end do
    SolverSettings%Output%Status=min(j,SolverSettings%MaxIter)
    SolverSettings%Output%AbsError=error
    SolverSettings%Output%RelError=error/norm   
    if (allocated(DiagIdx)) then
      deallocate(DiagIdx)
    end if
    if (allocated(InvDiag)) then
      deallocate(InvDiag)
    end if
    if (allocated(BRowIdx)) then
      deallocate(BRowIdx)
    end if
    if(.NOT.SolverSettings%Output%Success) ioX=StartX
    if(allocated(r))      deallocate(r)
    if(allocated(z))      deallocate(z)
    if(allocated(p))      deallocate(p)
    if(allocated(Ap))     deallocate(Ap)
    if(allocated(OldX))   deallocate(OldX)
    if(allocated(StartX)) deallocate(StartX)
    
    end subroutine SingleRHS

    subroutine VectorAddition(imode,iNuma,iNoRows,ioX,iY,alfa,ioZ,iV)
    !Tue, Dec 2014
    !This routine does different kinds of vector additions in parallel depending on the mode selected.
    !mode = 0 : ioX = iY
    !mode = 1 : ioX = iY - ioX
    !mode = 2 : ioX = ioX + alfa*iY
    !mode = 3 : ioX = alfa*ioX + iY  
    !mode = 4 : ioX(:) = ioX(:) + alfa*iY(:)    
    !           ioZ(:) = ioZ(:) - alfa*iV(:)
    !IO
    !imode -> selects the type of vector addition used
    !iNuma, -> Number of Numa Nodes
    !iNoRows -> Number of Rows
    !ioX,iY,ioZ,iV     -> vectors
    !alfa              -> scalar
    implicit none
    integer,intent(in)                          :: imode,iNuma,iNoRows
    real*8,intent(inout),dimension(:)           :: ioX
    real*8,intent(in),dimension(:)              :: iY
    real*8,intent(in),optional                  :: alfa
    real*8,intent(inout),dimension(:),optional  :: ioZ
    real*8,intent(in),dimension(:),optional     :: iV
    
    integer i
    if (iNuma.gt.1) then
      if (imode.eq.0) then
    !$OMP PARALLEL DEFAULT(NONE) SHARED(iNoRows,ioX,iY) PRIVATE(I)
    !$OMP DO        
          do i=1,iNoRows
            ioX(i) = iY(i)    
          end do
    !$OMP END DO 
    !$OMP END PARALLEL               
      elseif (imode.eq.1) then
    !$OMP PARALLEL DEFAULT(NONE) SHARED(iNoRows,ioX,iY) PRIVATE(I)
    !$OMP DO        
          do i=1,iNoRows
            ioX(i) = iY(i) - ioX(i)    
          end do
    !$OMP END DO 
    !$OMP END PARALLEL               
      elseif(imode.eq.2) then
    !$OMP PARALLEL DEFAULT(NONE) SHARED(iNoRows,ioX,iY,Alfa) PRIVATE(I)
    !$OMP DO        
          do i=1,iNoRows
            ioX(i) = ioX(i) + alfa*iY(i)    
          end do
    !$OMP END DO 
    !$OMP END PARALLEL               
      elseif(imode.eq.3) then
    !$OMP PARALLEL DEFAULT(NONE) SHARED(iNoRows,ioX,iY,Alfa) PRIVATE(I)
    !$OMP DO        
          do i=1,iNoRows
            ioX(i) = alfa*ioX(i) +iY(i)    
          end do
    !$OMP END DO 
    !$OMP END PARALLEL               
      elseif(imode.eq.4) then
    !$OMP PARALLEL DEFAULT(NONE) SHARED(iNoRows,ioX,iY,alfa,ioZ,iV) PRIVATE(I)
    !$OMP DO        
          do i=1,iNoRows
            ioX(i) = ioX(i) + alfa*iY(i) 
            ioZ(i) = ioZ(i) - alfa*iV(i)
          end do
    !$OMP END DO 
    !$OMP END PARALLEL  
      end if
      
    else
        if (imode.eq.0) then
          ioX(:) = iY(:)    
        elseif (imode.eq.1) then
          ioX(:) = iY(:) - ioX(:)    
        elseif(imode.eq.2) then
          ioX(:) = ioX(:) + alfa*iY(:)    
        elseif(imode.eq.3) then
          ioX(:) = alfa*ioX(:) + iY(:)    
       elseif(imode.eq.4) then
          ioX(:) = ioX(:) + alfa*iY(:)    
          ioZ(:) = ioZ(:) - alfa*iV(:)    
       end if   
    end if
    end subroutine VectorAddition
   
    subroutine Compute_Precondition_BLUT_real(iA,SolverSettings)
    !Tue November 2014
    !This routine sets up all the parameters needed for the Block LUT matrix factorization in parallel.
    !After everything has been setup, then the factorization algorithm is called.
    !IO
    !iA -> input matrix, in CSR format, and of the type TSparseMat
    !SolverSettings          -> contains all the relevant settings the solvers need, it is also where stats about the solve is saved. 
    !                           This needs to be set before calling this routine.
    use omp_lib
    implicit none
    type (TSolverSettings),intent(inout)           :: SolverSettings
    type (TSparseMat),     intent(inout)           :: iA
    
    !Local Vars
    integer                                :: i,Threads,NoRows,MaxRows,NoPrecondElements,SplitNoRows,ierr
    integer,dimension(:),allocatable           :: jw   
    integer,dimension(:,:),allocatable     :: BRowIdx
    real*8                                 :: Timer1, Timer2    
    real*8,dimension(:),allocatable        :: w
    Timer1=omp_get_wtime()     
    Threads=SolverSettings%NBlocks
    allocate(BRowIdx(iA%NoRows,2))
    call SMatGetBlockRowIdx(iA,SolverSettings%BlockS,BRowIdx)
    MaxRows = 0
    do i = 1 , threads
        MaxRows = Max(MaxRows,SolverSettings%BlockS(i,2)-SolverSettings%BlockS(i,1)+1)
    end do             
    NoPrecondElements=(2*SolverSettings%PrecondMaxFill+1)*MaxRows
    Allocate(SolverSettings%BPrecondVals(NoPrecondElements,Threads))
    Allocate(SolverSettings%BPrecondIdxs(NoPrecondElements,Threads))
    Allocate(SolverSettings%BPrecondDiag(MaxRows,Threads))
!$OMP PARALLEL IF(StartOpenMP(2)) DEFAULT(SHARED) PRIVATE(i,SplitNoRows,w,jw,ierr,NoPrecondElements)
!$OMP DO SCHEDULE(STATIC)
    do i = 1 , Threads
      SplitNoRows=SolverSettings%BlockS(i,2)-SolverSettings%BlockS(i,1)+1
      Allocate(w(SplitNoRows+1))
      Allocate(jw(2*SplitNoRows))
      NoPrecondElements=(2*SolverSettings%PrecondMaxFill)*SplitNoRows
      
      call BLUT_Factorization_real(SplitNoRows,iA%Vals,iA%ColIdxs,BRowIdx(SolverSettings%BlockS(i,1):SolverSettings%BlockS(i,2),:),&
          SolverSettings%PrecondMaxFill,SolverSettings%PrecondTolerance,SolverSettings%BPrecondVals(:,i),SolverSettings%BPrecondIdxs(:,i),&
          SolverSettings%BPrecondDiag(:,i),NoPrecondElements,w,jw,SolverSettings%BlockS(i,1)-1,ierr)
      !On Numa architechture we may now need to gather LU into one matrix and copy it around, depending whether we use PCG og BPCG routines.    
      deAllocate(w)
      deAllocate(jw)
    end do
!$OMP END DO
!$OMP END PARALLEL  
    Timer2=omp_get_wtime() 
    SolverSettings%Output%FactorizationTime=Timer2-Timer1 
    if (allocated(BRowIdx)) then
      deallocate(BRowIdx)
    end if  
    end subroutine Compute_Precondition_BLUT_real
    
   
!*************************************************************            
    
  subroutine BlockSplitter(iNoRows,iThreads,iParam,oBlocks)
    !Tue, September 2014
    !Splits iNoRows into iThreads Intervals.
    !Each interval will have a length of iParam*n, where n is a natural number.
    !The result will be saved in oBlocks, where oBlocks(i,1) is the start of the i'th interval and oBlocks(i,2) is the end of the interval
    !IO
    !Update 2015 Tue
    !Fixed the case if iNoRows < iThreads
    implicit none
    Integer,intent(in)                :: iNoRows,iThreads
    Integer,intent(in)                :: iParam
    Integer,intent(out),allocatable       :: oBlocks(:,:)
    !Local vars
    integer                           :: i
    integer                           :: BlockSize
    
    BlockSize=max(iNoRows/iThreads,1) 
    BlockSize =mod(IParam-mod(BlockSize,iParam),iParam)+BlockSize
    Allocate(oBlocks(iThreads,2))
    do i=1,iThreads
        oBlocks(i,1) = (i-1)*(BlockSize)+1
        oBlocks(i,2) = min(i*BlockSize,iNoRows)
    end do     
    print*,'Blocks',oBlocks
  end subroutine BlockSplitter  
  
  
subroutine BLUT_Factorization_real(n,iVals,iColIdxs,iBRowIdxs,lfil,droptol,alu,jlu,ju,iwk,w,jw,Offset,ierr)
!Tue November 2014
!This routine is designed to compute an Incomplete Cholesky factorization, of a sparse linear system, limited by  on a system.
!It has been made by modifying an old ILUT preconditioner made by Saad.
!A slightly modfified description of the ILUT preconditioner is found below. 
!
!IO
!
!n                          -> Number of rows in input matrix
!iVals,iColIdxs             -> Is the values and Column indexes of the input matrix stored in CSR format
!iBRowIdxs                  -> Is the cutoff row indices for the input matrix. iBRowIdxs(i,1) is the first element in row i we use, iBRowIdxs(i,2) is the last element in row i we use
!lfil                       -> The fill in parameter, each row of L and each row of U will have a maximum of lfil elements (excluding the diagonal element)
!droptol                    -> The minimum threshold before we drop a term.
!alu,jlu,ju                 -> the output LU matrix stored in modified sparse row format
!iwk                        -> the length of arrays alu,jlu
!w                          -> work array
!jw                         -> work array
!Offset                     -> The offset between the input matrix idx and the LU idx we want to save them in.
!ierr                       -> status of subroutine
!
!----------------------------------------------------------------------*
!                      *** ILUT preconditioner ***                     *
!      incomplete LU factorization with dual truncation mechanism      *
!----------------------------------------------------------------------*
!c     Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996  *
!c----------------------------------------------------------------------*
!c PARAMETERS                                                           
!c-----------                                                           
!c
!c on entry:
!c========== 
!c n       = integer. The row dimension of the matrix A. The matrix 
!c
!c iVals,iColIdxs, = matrix stored in Compressed Sparse Row format.              
!c iBRowIdxs is special!
!c lfil    = integer. The fill-in parameter. Each row of L and each row
!c           of U will have a maximum of lfil elements (excluding the 
!c           diagonal element). lfil must be .ge. 0.
!c           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
!c           EARLIER VERSIONS. 
!c
!c droptol = real*8. Sets the threshold for dropping small terms in the
!c           factorization. See below for details on dropping strategy.
!c
!c  
!c iwk     = integer. The lengths of arrays alu and jlu. If the arrays
!c           are not big enough to store the ILU factorizations, ilut
!c           will stop with an error message. 
!c
!c On return:
!c===========
!c
!c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!c           the L and U factors together. The diagonal (stored in
!c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!c           contains the i-th row of L (excluding the diagonal entry=1)
!c           followed by the i-th row of U.
!c
!c ju      = integer array of length n containing the pointers to
!c           the beginning of each row of U in the matrix alu,jlu.
!c
!c ierr    = integer. Error message with the following meaning.
!c           ierr  = 0    --> successful return.
!c           ierr .gt. 0  --> zero pivot encountered at step number ierr.
!c           ierr  = -1   --> Error. input matrix may be wrong.
!c                            (The elimination process has generated a
!c                            row in L or U whose length is .gt.  n.)
!c           ierr  = -2   --> The matrix L overflows the array al.
!c           ierr  = -3   --> The matrix U overflows the array alu.
!c           ierr  = -4   --> Illegal value for lfil.
!c           ierr  = -5   --> zero row encountered.
!c
!c work arrays:
!c=============
!c jw      = integer work array of length 2*n.
!c w       = real work array of length n+1.
!c  
!c----------------------------------------------------------------------
!c w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] 
!c jw(n+1:2n)  stores nonzero indicators
!c 
!c Notes:
!c ------
!c The diagonal elements of the input matrix must be  nonzero (at least
!c 'structurally'). 
!c
!c----------------------------------------------------------------------* 
!c---- Dual drop strategy works as follows.                             *
!c                                                                      *
!c     1) Theresholding in L and U as set by droptol. Any element whose *
!c        magnitude is less than some tolerance (relative to the abs    *
!c        value of diagonal element in u) is dropped.                   *
!c                                                                      *
!c     2) Keeping only the largest lfil elements in the i-th row of L   * 
!c        and the largest lfil elements in the i-th row of U (excluding *
!c        diagonal elements).                                           *
!c                                                                      *
!c Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
!c keeping  the largest  elements in  each row  of L  and U.   Taking   *
!c droptol .ne.  0 but lfil=n will give  the usual threshold strategy   *
!c (however, fill-in is then mpredictible).                             *
!c----------------------------------------------------------------------*      
      implicit none  
      integer,intent(in)                :: n,lfil,iwk,offset
      real*8,intent(in),dimension(:)    :: iVals
      real*8,intent(inout),dimension(:) :: alu,w
      real*8,intent(in)                 :: droptol
      
      integer,intent(in),dimension(:)   :: iColIdxs
      integer,intent(in),dimension(:,:) :: iBRowIdxs
      integer,intent(inout),dimension(:):: jlu,ju,jw
      integer,intent(inout)             :: ierr

!c     locals
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,length 
      real*8 tnorm, t, abs, s, fact 
      if (lfil .lt. 0) goto 998
!c-----------------------------------------------------------------------
!c     initialize ju0 (points to next element to be added to alu,jlu)
!c     and pointer array.
!c-----------------------------------------------------------------------
      ju0 = n+2
      jlu(1) = ju0
!c
!c     initialize nonzero indicator array. 
!c
      do j=1,n
         jw(n+j)  = 0
      end do
!c-----------------------------------------------------------------------
!c     beginning of main loop.
!c-----------------------------------------------------------------------
do ii = 1, n
         j1 = iBRowIdxs(ii,1)
         j2 = iBRowIdxs(ii,2)
        !j1=iRowIdxs(ii)
        !j2=iRowIdxs(ii+1)-1
         tnorm = 0.0d0
         do k=j1,j2
            tnorm = tnorm+abs(iVals(k))
         end do
         if (tnorm .eq. 0.0) goto 999
         tnorm = tnorm/real(j2-j1+1)
!c     
!c     unpack L-part and U-part of row of A in arrays w 
!c     
         lenu = 1
         lenl = 0
         jw(ii) = ii
         w(ii) = 0.0
         jw(n+ii) = ii

         do j = j1, j2
            k = iColIdxs(j)-offset
            t = iVals(j)
            if (k .lt. ii) then
               lenl = lenl+1
               jw(lenl) = k
               w(lenl) = t
               jw(n+k) = lenl
            else if (k .eq. ii) then
               w(ii) = t
            else
               lenu = lenu+1
               jpos = ii+lenu-1 
               jw(jpos) = k
               w(jpos) = t
               jw(n+k) = jpos
            endif
         end do
         jj = 0
         length = 0 
!c     
!c     eliminate previous rows
!c     
 150     jj = jj+1
         if (jj .gt. lenl) goto 160
!c-----------------------------------------------------------------------
!c     in order to do the elimination in the correct order we must select
!c     the smallest column index among jw(k), k=jj+1, ..., lenl.
!c-----------------------------------------------------------------------
         jrow = jw(jj)
         k = jj
!c     
!c     determine smallest column index
!c     
         do j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            end if
        end do
         
         if (k .ne. jj) then
!c     exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
!c     exchange in jr
            jw(n+jrow) = jj
            jw(n+j) = k
!c     exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
!c
!c     zero out element in row by setting jw(n+jrow) to zero.
!c     
         jw(n+jrow) = 0
!c
!c     get the multiplier for row to be eliminated (jrow).
!c     
         fact = w(jj)*alu(jrow)
         if (abs(fact) .le. droptol) goto 150
!c     
!c     combine current row and row jrow
!c
         do k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
            j = jlu(k)
            jpos = jw(n+j)
            if (j .ge. ii) then
!c     
!c     dealing with upper part.
!c     
               if (jpos .eq. 0) then
!c
!c     this is a fill-in element
!c     
                  lenu = lenu+1
                  if (lenu .gt. n) goto 995
                  i = ii+lenu-1
                  jw(i) = j
                  jw(n+j) = i
                  w(i) = - s
               else
!c
!c     this is not a fill-in element 
!c
                  w(jpos) = w(jpos) - s

               endif
            else
!c     
!c     dealing  with lower part.
!c     
               if (jpos .eq. 0) then
!c
!c     this is a fill-in element
!c     
                  lenl = lenl+1
                  if (lenl .gt. n) goto 995
                  jw(lenl) = j
                  jw(n+j) = lenl
                  w(lenl) = - s
               else
!c     
!c     this is not a fill-in element 
!c     
                  w(jpos) = w(jpos) - s
               endif
            endif
        end do
!c     
!c     store this pivot element -- (from left to right -- no danger of
!c     overlap with the working elements in L (pivots). 
!c     
         length = length+1 
         w(length) = fact
         jw(length)  = jrow
         goto 150
         
 160     continue
!c     
!c     reset double-pointer to zero (U-part)
!c     
         do k=1, lenu
            jw(n+jw(ii+k-1)) = 0
         end do
!c     
!c     update L-matrix
!c     
         lenl = length 
         length = min0(lenl,lfil)
!c     
!c     sort by quick-split
!c
         call qsplit (w,jw,lenl,length)
!c
!c     store L-part
!c 
         do k=1, length 
            if (ju0 .gt. iwk) goto 996
            alu(ju0) =  w(k)
            jlu(ju0) =  jw(k)
            ju0 = ju0+1
        end do
!c     
!c     save pointer to beginning of row ii of U
!c     
         ju(ii) = ju0
!c
!c     update U-matrix -- first apply dropping strategy 
!c
         length = 0
         do k=1, lenu-1
            if (abs(w(ii+k)) .gt. droptol*tnorm) then 
               length = length+1
               w(ii+length) = w(ii+k) 
               jw(ii+length) = jw(ii+k) 
            endif
         enddo
         lenu = length+1
         length = min0(lenu,lfil)

         call qsplit(w(ii+1:ii+lenu), jw(ii+1:ii+lenu), lenu-1,length)
!c
!c     copy
!c 
         t = abs(w(ii))
         if (length + ju0 .gt. iwk) goto 997
         do k=ii+1,ii+length-1 
            jlu(ju0) = jw(k)
            alu(ju0) = w(k)
            t = t + abs(w(k) )
            ju0 = ju0+1
         end do
!c     
!c     store inverse of diagonal element of u
!c     
         if (w(ii) .eq. 0.0) w(ii) = (0.0001 + droptol)*tnorm
     
         alu(ii) = 1.0d0/ w(ii) 
!c     
!c     update pointer to beginning of next row of U.
!c     
         jlu(ii+1) = ju0
!c-----------------------------------------------------------------------
!c     end main loop
!c-----------------------------------------------------------------------
 end do
      ierr = 0
      return
!c
!c     incomprehensible error. Matrix must be wrong.
!c     
 995  ierr = -1
      return
!c     
!c     insufficient storage in L.
!c     
 996  ierr = -2
      return
!c     
!c     insufficient storage in U.
!c     
 997  ierr = -3
      return
!c     
!c     illegal lfil entered.
!c     
 998  ierr = -4
      return
!c     
!c     zero row encountered
!c     
 999  ierr = -5
      return
!c----------------end-of-ilut--------------------------------------------
end subroutine BLUT_Factorization_real

  
subroutine Apply_BPrecondition_BLUT(iBlocks, iRHS, oX, iLUVals, iLUIdxs, iLUDiag) 
    !Tue November 2014
    !This routine does a forward and backward block solve with the Block LUT matrix
    !
    !IO
    !
    !iBlocks                 -> The blocksplitting of the Matrix
    !iRHS                    -> The righthandside we update our solution after
    !oX                      -> The output solution
    !iLUVals,iLUIdxs,iLUDiag -> The Block LU matrix where each block is a separate matrix stored in modfied sparse row format.
	integer,intent(in),dimension(:,:) :: iBlocks
    real*8,intent(in),dimension(:)    :: iRHS
    real*8,intent(out),dimension(:)   ::  oX
    real*8,intent(in),dimension(:,:)  :: iLUVals
    integer,intent(in),dimension(:,:) :: iLUIdxs, iLUDiag    
    !local variables
    integer             :: i,k,j,threads,m,tmp
    
    
    threads=size(iBlocks,1)
!$OMP PARALLEL IF(StartOpenMP(2)) DEFAULT(NONE)  SHARED(iBlocks,threads,iRHS,oX,iLUVals,iLUIdxs,iLUDiag) PRIVATE(k,i,j,m,tmp)
!$OMP DO SCHEDULE(STATIC)
    do k=1,threads
        do i = iBlocks(k,1),iblocks(k,2)
           oX(i) = iRHS(i)
        end do
!
! forward solve (with U^T)
!
       tmp=iBlocks(k,1)-1
       m=0
        do i = iBlocks(k,1),iblocks(k,2)
           m=m+1
           oX(i) = oX(i) * iLUVals(m,k)
           do j=iLUDiag(m,k),iLUIdxs(m+1,k)-1
              oX(tmp+iLUIdxs(j,k)) = oX(tmp+iLUIdxs(j,k)) - iLUVals(j,k)* oX(i)
           end do
        end do
!     
!     backward solve (with L^T)
!    
      do i = iBlocks(k,2),iBlocks(k,1),-1
       m=i-tmp
	   do j=iLUIdxs(m,k),iLUDiag(m,k)-1
              oX(tmp+iLUIdxs(j,k)) = ox(tmp+iLUIdxs(j,k)) - iLUVals(j,k)*oX(i)
       end do
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL      
end subroutine Apply_BPrecondition_BLUT  
  
   subroutine Apply_BPrecondition_BSGS(iColIdxs,iRowIdxs,iVals,iBRowIdx,iDiagIdx,iInvDiag,iRhS,ioX,Blocks)
    !Tue, August 2014
    !This routine does a forward and backward block solve with the Block SOR preconditioner
    !This algorithm is described in 266-267 in the book "Iterative methods for sparse linear systems" in the normal non-block edition.
    !Note that the SGS preconditioner is not calculated beforehand, but rather calculated on the fly during the applying of it.
    !IO
    !ispA     -> Input, Sparse matrix of the linear system
    !iDiagIdx -> Input, a vector containing indexes to the diagonal entries in the CSR matrix ispA.
    !iInvDiag -> Input, a vector containing the inverse values of the diagonal element in ispA.   
    !iRhS     -> Input, the right hand side of the linear system we try to solve.
    !ioX      -> Input, Solution estimate to update iteratively. Output, our new updated solution
    !Blocks   -> Input, contains the blocksplitting of ispA Blocks(i,1) contains the rowIdx of the first row in the i'th block, Blocks(i,2) contains the last.
    implicit none
    integer,intent(in),dimension(:)   :: iColIdxs,iRowIdxs
    real*8,intent(in),dimension(:)    :: iVals
    integer,intent(in),dimension(:,:) :: iBRowIdx
    integer,intent(in),dimension(:)   :: iDiagIdx
    real*8,intent(in),dimension(:)    :: iInvDiag
    real*8,intent(in),dimension(:)    :: iRHS
    real*8,intent(inout),dimension(:) :: ioX
    integer,intent(in)                :: Blocks(:,:)
    !Local vars
    integer             :: i,j,k
    real*8              :: tmp
    real*8              :: Factor
    integer             :: threads
    real*8              :: omega
    threads=size(Blocks,1)
    
!$OMP PARALLEL IF(StartOpenMP(2)) DEFAULT(NONE) SHARED(omega,Blocks,iDiagIdx,iInvDiag,iColIdxs,iRowIdxs,iVals,ioX,irhs,threads,iBRowIdx) PRIVATE(i,k,j)
!c$OMP CRITICAL
!$OMP DO SCHEDULE(STATIC)
    do k=1,threads
      !We start by doing the  forward step
      do i = Blocks(k,1),Blocks(k,2)
        ioX(i) = iRHS(i)
        do j = iBRowIdx(i,1) , iDiagIdx(i)-1
            ioX(i) = ioX(i) - iVals(j)*iInvDiag(iColIdxs(j))*ioX(iColIdxs(j))
        end do
      end do

      !Now comes the backward step    
      do i = Blocks(k,2),Blocks(k,1),-1
        do j = iBRowIdx(i,2) ,iDiagIdx(i)+1,-1
              ioX(i) = ioX(i) - iVals(j)*ioX(iColIdxs(j))
        end do
        ioX(i) = ioX(i)*iInvDiag(i)        
      end do
    end do  
!$OMP END DO
!c$OMP END CRITICAL
!$OMP END PARALLEL
   end subroutine Apply_BPrecondition_BSGS     
   
  
  subroutine    SetSolverSettings(SolverSettings,A,iUseRCM,iUseNormalization,iSolverType,iIter,iBlockSize,iFillingFactor,iStepLimit)
    !Tue Dec, 2014  
    !This routine Sets the solversettings prior to a linear solve.
    !NOTE that this also calculates the precondition matrix, which the caller takes responsibility to free. The freeing should be done by calling the FreeSolverSettings subroutine.
    !any optional value you wish chosen by the subroutine should be set to -1
    !IO
    !SolverSettings          -> This is where we save all the relevant settings about the solve we want to do - note if we have preconditioner this is also saved in here.
    !A                       -> input matrix, in CSR format, and of the type TSparseMat
    !(Optional)   The rest are optional parameters which can be chosen manually, any parameter set to -1 will be chosen automatically.  
    !iUseNormalization       -> (optional) Use Normalization 
    !iSolverType             -> (optional) manually chooses the solvertype, if none is selected or set to -1, a fitting type will be chosen based on the problem.
    !                           1  SGS
    !                           2  BLUT
    !                           3  Pardiso (not implemented here)
    !iIter                   -> (optional) manually set the maximum iterations. (default is DEFAULT_MAXITERATIONS)
    !iBaseBlockSize          -> (optional) The blocksize will be a multiplum of this number for all blocks if specified.
    !                            this is used for blocksplitting to make sure we don't split in the middle of a block.
    !iFillingFactor          -> (optional) Determines the number of elements allowed in the IC factorization. Factor is based on the average number of elements in iA; default value 1.5
    !iStepLimit              -> (optional) sets the relative minimum difference an iterative step should make in the solution 
    !                           This could be important to set when using block solvers where convergence to the solution is not guaranteed
    use omp_lib
    use msMat
    use mRCM
    implicit none  
    type (TSolverSettings),intent(inout)           :: SolverSettings
    type (TSparseMat),     intent(inout)           :: A
    integer, intent(in),optional                   :: iUseRCM
    integer, intent(in),optional                   :: iUseNormalization
    integer, intent(in),optional                   :: iSolverType
    integer, intent(in),optional                   :: iIter 
    integer, intent(in),optional                   :: iBlockSize      
    real*8, intent(in),optional                   :: iFillingFactor
    integer, intent(in),optional                   :: iStepLimit
    character*256                     :: S
    real*8                            :: NParFactor,DampCutOff,StartTolerance,EndTolerance,DampVal
    real*8,allocatable                :: Diag(:)
    logical                           :: UseSlowButProvenSettings=.false.
    logical                           :: IsLayeredInv
    integer                           :: ErrorCode,Param,NModPar,NParSec,i,j,MaxWidth
    integer                           :: OldThreads
    
    real*8                            :: Elements_pr_row,DiagDom,RunTimeBegin,RunTimeEnd
    !Tweakable factors
    integer,parameter                 :: DEFAULT_USERCM = 0
    integer,parameter                 :: DEFAULT_USENORMALIZATION = 1
    integer,parameter                 :: DEFAULT_INFORMATIONLEVEL = 1
    integer,parameter                 :: DEFAULT_NORHS = 1
    integer,parameter                 :: DEFAULT_MAXITERATIONS = 500
    
    real*8,parameter                  :: DEFAULT_STEPLIMIT = 0 
    real*8,parameter                  :: DEFAULT_FILLINGFACTOR=2.5
    real*8,parameter                  :: MAXDIAGDOM=1        !Not made to be adjusted just yet    
    real*8,parameter                  :: SGS_THRESHOLD=0.1  
    
    logical,parameter                 :: DEBUG_MODE = .FALSE.
    !Test
    real*8                            :: DiagDomTest,test
    integer                           :: NumberOfBlocks
    integer                           :: Error
    real*8,dimension(:),allocatable   :: DummyVal
    logical :: Dummy
    
    
    solverSettings%UseRCM=LoadOptionalParam(-1,iUseRCM)
    if(solverSettings%UseRCM.eq.-1) solverSettings%UseRCM=DEFAULT_USERCM
    if (solverSettings%UseRCM.eq.1) then
      if(allocated(solverSettings%RCMIndices)) deallocate(solverSettings%RCMIndices)   
      allocate(solverSettings%RCMIndices(A%NoRows))  
      call RCMSort(A,solverSettings%RCMIndices)
    end if
    
    solverSettings%UseNormalization=LoadOptionalParam(-1,iUseNormalization)
    if(solverSettings%UseNormalization.eq.-1) solverSettings%UseNormalization=DEFAULT_USENORMALIZATION
    if(present(iSolverType)) then  ! remove normalization for direct solver
        if(iSolverType.eq.3) solverSettings%UseNormalization = 0
    endif
    
    SolverSettings%FillingFactor=LoadOptionalParam(-1,iFillingFactor)
    if(SolverSettings%FillingFactor.eq.-1) SolverSettings%FillingFactor=DEFAULT_FILLINGFACTOR
    
    if (solverSettings%UseNormalization) then
      if(allocated(SolverSettings%ReNormalization)) deallocate(SolverSettings%ReNormalization)   
      if(allocated(SolverSettings%Normalization)) deallocate(SolverSettings%Normalization)   
      allocate(SolverSettings%ReNormalization(A%NoRows))  
      allocate(SolverSettings%Normalization(A%NoRows))
      call SMatGetDiag(A,SolverSettings%ReNormalization,.TRUE.)
      SolverSettings%ReNormalization = sqrt(SolverSettings%ReNormalization)
      SolverSettings%Normalization = 1d0/SolverSettings%ReNormalization      
      call SMatDAD(A,SolverSettings%Normalization)
    end if
    
    SolverSettings%IsPrepared = .False.
    SolverSettings%output%Success = .False.
    
    SolverSettings%IsComplex = A%IsComplex
    SolverSettings%BlockSize=LoadOptionalParam(-1,iBlockSize)
    SolverSettings%RelTol=1e-12
    SolverSettings%AbsTol=1e-9
      
    SolverSettings%PrecondTolerance=1e-6

    SolverSettings%Steplimit=LoadOptionalParam(-1,iSteplimit)
    if(SolverSettings%Steplimit.eq.-1) SolverSettings%Steplimit=DEFAULT_STEPLIMIT
      
            
    SolverSettings%MaxIter=LoadOptionalParam(-1,iIter)
    if(SolverSettings%MaxIter.eq.-1) SolverSettings%MaxIter=DEFAULT_MAXITERATIONS     
      
    if(.not.(iSolverType.eq.3)) call SMatGetDiagonalDominant(A,SolverSettings%DiagDom)
      
    !SolverType is determined
    SolverSettings%SolverType=LoadOptionalParam(-1,iSolverType)
    if (SolverSettings%SolverType.le.0)  call SelectSolver(A,SolverSettings,SolverSettings%DiagDom,SGS_threshold)
    
	if(.not.(isolvertype.eq.3)) then
      Elements_pr_row = (1d0*A%NoElements)/(1d0*A%NoRows)
      !Maximum number of precondition elements is chosen
      MaxWidth=SMatFindMaxColDif(A)
      SolverSettings%PrecondMaxFill=NINT(sqrt(SolverSettings%DiagDom)*SolverSettings%FillingFactor*Elements_pr_row)
      SolverSettings%PrecondMaxFill=min(A%NoRows,min(SolverSettings%PrecondMaxFill,2*MaxWidth))
    endif 
	
    !DEBUG_MODE Only relevant in debug mode
    if (DEBUG_MODE) then
        print*,'IsComplex',SolverSettings%IsComplex
        print*,'Blocksize',SolverSettings%BlockSize
        print*,'Diagdom',SolverSettings%DiagDom
        print*,'SolverType',SolverSettings%SolverType
        print*,'SolverSettings%PrecondMaxFill',SolverSettings%PrecondMaxFill
        print*,'Elements_pr_row',Elements_pr_row     
    end if
      
    !reset pointers
    if (allocated(SolverSettings%PrecondVals))  deallocate(SolverSettings%PrecondVals,Stat=ErrorCode) 
    if (allocated(SolverSettings%PrecondCols))  deallocate(SolverSettings%PrecondCols,Stat=ErrorCode) 
    if (allocated(SolverSettings%PrecondRows))  deallocate(SolverSettings%PrecondRows,Stat=ErrorCode) 
    if (allocated(SolverSettings%BPrecondVals)) deallocate(SolverSettings%BPrecondVals,Stat=ErrorCode) 
    if (allocated(SolverSettings%BPrecondIdxs)) deallocate(SolverSettings%BPrecondIdxs,Stat=ErrorCode) 
    if (allocated(SolverSettings%BPrecondDiag)) deallocate(SolverSettings%BPrecondDiag,Stat=ErrorCode) 
    if (allocated(SolverSettings%BlockS))       deallocate(SolverSettings%BlockS,Stat=ErrorCode) 
    if (allocated(SolverSettings%PrecondcVals)) deallocate(SolverSettings%PrecondcVals,Stat=ErrorCode) 
   
    Dummy=StartOpenMP(2)
    NumberOfBlocks=omp_get_max_threads()
    if(.not.(iSolvertype.eq.3)) call BlockSplitter(A%NoRows,NumberOfBlocks,SolverSettings%BlockSize,SolverSettings%BlockS)
    SolverSettings%NBlocks=size(SolverSettings%BlockS,1)
    !Factorization of the preconditioner is done
    RuntimeBegin =  omp_get_wtime()    
    SELECT CASE(SolverSettings%SolverType)    
    CASE(2)
        !BLUT
        call Compute_Precondition_BLUT_real(A,SolverSettings)
    CASE(3)
        !Pardiso (not implemented here)        
    END SELECT
    RuntimeEnd =  omp_get_wtime()    
    SolverSettings%output%FactorizationTime=RuntimeEnd-RunTimeBegin   
    print*,'Factorization Time:',SolverSettings%output%FactorizationTime
    SolverSettings%IsPrepared=.true. 
    return 
  end subroutine SetSolverSettings
 
     subroutine SolverNormalizeMatrix(SolverSettings,ioA,Norm)
! Kristoffer Andersen, Dec 2016
! Normalizes the input matrix ioA using the normalization vectors in 
! solversettings. There is no check on the size of the vectors compared to the matrix
! Norm controls whether normalization or the renormalization vector is used. 
! If SolverSettings%UseNormalization = .false. nothing is performed
! IO
! I - SolverSettings - settings with the normalization vectors.
!   If normalization is not used nothing happens to ioA
! IO - ioA - the matrix which is transformed to D*A*D
! I,optional - Norm - controls whether normalization or renormalization is used. 
!   .true. (default) for normalization
!   .false. for renormalization
!
    implicit none
    type(TSolverSettings),intent(in)     :: SolverSettings 
    type(Tsparsemat),intent(inout)       :: ioA
    logical,optional,intent(in)          :: Norm
    !var
    logical  :: normalize
    if(SolverSettings%USENORMALIZATION) then ! only if normalization is used
        normalize = .true.
        if(present(Norm)) normalize = norm
        
        if(normalize) then
            call SMatDAD(ioA,SolverSettings%Normalization)
        else
            call SMatDAD(ioA,SolverSettings%ReNormalization)
        endif
    endif
    
    RETURN
    END subroutine SolverNormalizeMatrix
  
  
  subroutine SelectSolver(iA,SolverSettings,iDiagDom,iSGS_threshold)
   !Tue Dec, 2014  
   !This routine selects the most appropiete solvertype based on the diagonal dominance and various thresholds.
   !IO
   !iA                      -> input matrix, in CSR format, and of the type TSparseMat
   !SolverSettings          -> This is where we save all the relevant settings about the solve we want to do.
   !iDiagDom                -> Maximum diagonal dominance factor for the matrix.
   !iSGS_threshold          -> threshold between BLUT and SGS preconditioner
    implicit none  
    type (TSparseMat),intent(in)           :: iA
    type(TSolverSettings),intent(inout)    :: SolverSettings
    real*8,intent(in)                      :: iDiagDom
    real*8,intent(in)                      :: iSGS_threshold    
    
        if (SolverSettings%IsComplex) then
              SolverSettings%SolverType = 12 !LU0
        else
          if (iDiagDom.gt.iSGS_threshold) then
            !BLUT
            SolverSettings%SolverType=2
          else
            !SGS
            SolverSettings%SolverType=1            
          end if    
        end if
  end subroutine SelectSolver
  
  
  
    subroutine CheckSolverSettings(iSolverSettings)
    !Tue, dec 2014
    !Checks whether the solversettings have been set and if informationlevel is sufficiently high it also reports the solver settings
    !
    !use mError
    implicit none
    type(TSolverSettings),intent(in)    :: iSolverSettings
    character*200                       :: S
    integer         :: Output
    logical         :: PrintToLog
    
    if(.NOT.iSolverSettings%IsPrepared) goto 901
        
    return
    901 continue
    print*, 'Critical error! Solversettings not set, make sure SetSolverSettings is called before running the solver'
    stop
    end subroutine CheckSolverSettings

 
    logical function isSolverNormalized(SolverSettings)
    ! Kristoffer Andersen, Dec 2016
    ! return the normalization state of the solver
    implicit none
    type(TSolverSettings),intent(in)     :: SolverSettings 
    isSolverNormalized = SolverSettings%USENORMALIZATION
    RETURN
    END FUNCTION isSolverNormalized
    
   
     subroutine FreeSolverSettings(SolverSettings)
    !Tue June 2015
    !This routine frees the memory allocated by SolverSettings.
    !This routine should always be called whenever Solversettings is no longer needed, since caller takes resposibility to free TSolverSettings.    
    !
    implicit none
    type(TSolverSettings),intent(inout)     :: SolverSettings 
    integer                                 :: ErrorCode
	real*8                                  :: dummyVal(1)
    integer                                 :: dummyInt(1)
    integer                                 :: dummy

    !reset pointers
    if (allocated(SolverSettings%PrecondVals))  deallocate(SolverSettings%PrecondVals,Stat=ErrorCode) 
    if (allocated(SolverSettings%PrecondCols))  deallocate(SolverSettings%PrecondCols,Stat=ErrorCode) 
    if (allocated(SolverSettings%PrecondRows))  deallocate(SolverSettings%PrecondRows,Stat=ErrorCode) 
    if (allocated(SolverSettings%BPrecondVals)) deallocate(SolverSettings%BPrecondVals,Stat=ErrorCode) 
    if (allocated(SolverSettings%BPrecondIdxs)) deallocate(SolverSettings%BPrecondIdxs,Stat=ErrorCode) 
    if (allocated(SolverSettings%BPrecondDiag)) deallocate(SolverSettings%BPrecondDiag,Stat=ErrorCode) 
    if (allocated(SolverSettings%BlockS))       deallocate(SolverSettings%BlockS,Stat=ErrorCode) 
    if (allocated(SolverSettings%PrecondcVals)) deallocate(SolverSettings%PrecondcVals,Stat=ErrorCode) 
    if (allocated(SolverSettings%Renormalization)) deallocate(SolverSettings%Renormalization,Stat=ErrorCode) 
    if (allocated(SolverSettings%normalization)) deallocate(SolverSettings%normalization,Stat=ErrorCode) 
        
     end subroutine FreeSolverSettings
 
  
        subroutine qsplit(a,ind,n,ncut)
        implicit none
        real*8 a(n)
        integer ind(n), n, ncut,j,mid
!c-----------------------------------------------------------------------
!c     does a quick-sort split of a real array.
!c     on input a(1:n). is a real array
!c     on output a(1:n) is permuted such that its elements satisfy:
!c
!c     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
!c     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
!c
!c     ind(1:n) is an integer array which permuted in the same way as a(*).
!c-----------------------------------------------------------------------
        real*8 tmp, abskey
        integer itmp, first, last
!c-----
        first = 1
        last = n
        if (ncut .lt. first .or. ncut .gt. last) return
!c
!c     outer loop -- while mid .ne. ncut do
!c
 1      mid = first
        abskey = abs(a(mid))
        do 2 j=first+1, last
           if (abs(a(j)) .gt. abskey) then
              mid = mid+1
!c     interchange
              tmp = a(mid)
              itmp = ind(mid)
              a(mid) = a(j)
              ind(mid) = ind(j)
              a(j)  = tmp
              ind(j) = itmp
           endif
 2      continue
!c
!c     interchange
!c
        tmp = a(mid)
        a(mid) = a(first)
        a(first)  = tmp
!c
        itmp = ind(mid)
        ind(mid) = ind(first)
        ind(first) = itmp
!c
!c     test for while loop
!c
        if (mid .eq. ncut) return
        if (mid .gt. ncut) then
           last = mid-1
        else
           first = mid+1
        endif
        goto 1
!c----------------end-of-qsplit------------------------------------------
!c-----------------------------------------------------------------------
        end

     
    end module mSolver      
