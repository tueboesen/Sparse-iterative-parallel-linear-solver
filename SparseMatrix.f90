module mSMat
  use mParallel, only: StartOpenMP, ReturnOldThreadNum
  implicit none 
!  logical   :: UseDenseAlgebra,VoxelInversion
  integer,private                :: FInitNoElements=1000000     !Initial allocation size of the matrices
  logical                        :: FCheckDoubleEntries=.false. !Debug var. When .true. the code checks for double entries 
  logical                        :: FShowMessages=.false.       !when adding new values to the matrix and checks that the column indices of the 
  logical,private                :: ColIdxCorrected=.false.     !Internal variable ensuring that column index correction is only applied once. VERY IMPORTANT
                                                                !TO TOGGLE THIS VALUE TRUE IN THE BEGINNING OF ROUTINES THAT ALTER THE INDEX AND FALSE IN THE END!!!!!  
  logical,private                :: UseMKLSparseBLAS=.true.     !If true the intel math kernel library is used for high performance sparse BLAS. If false, full source routines from SparseKit are used.                                                                 
  logical                        :: SMatLogParameters=.true. 
!$OMP THREADPRIVATE(ColIdxCorrected)
                                                                
  type TSparseMat
  !Public variables(General matrix information)                                                               
    character*(256)              :: Label
    integer                      :: NoElements               !Number of elements in the matrix 
    integer                      :: NoRows                   !Number of rows of the matrix
    integer                      :: NoCols                   !Number of columns of the matrix 
    
    logical                      :: IsComplex=.false.        !True if the matrix is complex                                                            

    real*8 , dimension(:),pointer    :: Vals=>null()             !Values of the matrix elements if Real
    Complex*16, dimension(:),pointer :: cVals=>null()            !Values of the matrix elements if Complex 
    integer,  dimension(:),pointer   :: RowIdxs =>null()         !Row indices of the elements
    integer,  dimension(:),pointer   :: ColIdxs=>null()          !Column indices of the elements    
  end type TSparseMat
  
  
    contains
    subroutine SMatCreate(ioSparseMat,iNoRows,iNoCols,iLabel,iIsComplex,iAllocSize)
!Casper Kirkegaard, November 2009.
!This routine creates a TSparseMat structure, ie. initializes
!parameters allocates memory etc.
!
!Updated 2014 Tue, Now it takes complex numbers as well
!
! CK, october 2015.    
!  Removed initialization of value and column index vectors. It takes a long time when creating small matrices and is not necessary.
! KRA, Dec 2016. Added optional allocation size. Has to be >= 1
!IO
!   oSparseMat  -> Sparse matrix to create
!   iNoRows     -> Number of rows in the matrix
!   iNoCols     -> Number of cols in the matrix 
!   iLabel      -> Label of the matrix. Setting this to eg. 'G' or 'B' is very useful
!                  when debugging or optimizing 
!   iIsComplex  -> Default=.false., if true it means that the matrix is complex
!   iAllocSize  -> Default=FInitNoElements, is the allocation size. The default value is a 
!                  model parameter. Has to >=1    
    implicit none
    type (TSparseMat),intent(inout) :: ioSparseMat
    integer,intent(in)              :: iNoRows
    integer,intent(in)              :: iNoCols
    character*(*),intent(in)        :: iLabel          
    logical,intent(in),optional     :: iIsComplex    
    integer,intent(in),optional     :: iAllocSize
    integer                         :: FInitNoElements=1000000     !Initial allocation size of the matrices
!   var
    integer :: NoAlloc
!Implementation
    !If the matrix already exists, destroy it and recreate.
    if (associated(ioSparseMat%Vals)) then
      call SMatDestroy(ioSparseMat)
    end if
    !First initialize the vars of the input matrix
    ioSparseMat%IsComplex=.false.
    NoAlloc = FInitNoElements

    if (present(iIsComplex)) ioSparseMat%IsComplex=iIsComplex
    if (present(iAllocSize)) then
        if(iAllocSize>0) NoAlloc = iAllocSize
    endif
    
    ioSparseMat%NoRows=iNoRows
    ioSparseMat%NoCols=iNoCols
    ioSparseMat%Label=iLabel
    ioSparseMat%NoElements=0
    !Now allocate an initial size for the arrays holding the matrix
    allocate(ioSparseMat%ColIdxs(NoAlloc))
    if (ioSparseMat%IsComplex) then
      allocate(ioSparseMat%cVals(NoAlloc))
    else
      allocate(ioSparseMat%Vals(NoAlloc))
    end if    
    allocate(ioSparseMat%RowIdxs(iNoRows+1)) 
    ioSparseMat%RowIdxs(:)=0
    end subroutine SMatCreate

    
    subroutine SMatWriteSparse(iA,iV,iNDigits,iName,iLabel)
    !Gianluca Fiandaca, April 2013, AU
    !This routine writes the sparse matrix iA in 3 different files, one containing
    !the matrix row indices, one containing the matrix column indices and one containing 
    !the matrix values. The names of the output files are built concatenating the file name 
    !iName and the words "RowIdxs.txt", "ColIdxs.txt" & "Values.txt" 
    !IO
    !  iA          -> Sparse matrix to be written
    !  iV          -> RHS vector
    !  iName       -> the name of the 3 output files 
    !  iLabel      -> Label to be written before the matrix
    !  iNDigits    -> Number of digits after comma in the output format (in scientific notation) 
    implicit none
    type (TSparseMat),intent(inout)        :: iA
    real*8,dimension(:)                    :: iV  
    integer          ,intent(in), optional :: iNDigits
    character*(*)    ,intent(in), optional :: iName
    character*(*)    ,intent(in), optional :: iLabel
    
    integer                         :: i,j,k,Error
    integer                         :: NIntegerDigits,NScientificDigits
    integer,parameter               :: F=99
    real*8                          :: Val
    character*256                   :: RowName,ColName,ValName,S,DefName,RHSName
    integer(kind=1)                 :: CSR
 
    !digits to be used in the formats
    if(present(iNDigits)) then
      NScientificDigits=iNDigits  
    else
      NScientificDigits=4  
    end if
    NIntegerDigits=nint(log10(1.d0*iA%NoElements))+1
    !names and legend
    if(present(iName)) then
      DefName=trim(iName)//'_Definition.txt'
      RHSName=trim(iName)//'_RHS.txt'
      RowName=trim(iName)//'_RowIdxs.txt'
      ColName=trim(iName)//'_ColIdxs.txt'
      ValName=trim(iName)//'_Values.txt'
    else
      DefName='Definition.txt'  
      RHSName='RHS.txt'
      RowName='RowIdxs.txt'
      ColName='ColIdxs.txt'
      ValName='Values.txt'
    end if
    
    !Definition file
    open(unit=F,file=DefName,status='unknown')
    if(present(iLabel)) then
        write(F,'(A)') trim(iLabel)
    else 
        write(F,'(A)')
    end if
    write(F,30) iA%NoRows,iA%NoCols,iA%NoElements,CSR
    close(F)
    
    
    !RHS
    open(unit=F,file=RhsName,status='unknown')
    if(present(iLabel)) then
        write(F,'(A)') trim(iLabel)
    else 
        write(F,'(A)')
    end if
    do i=1,iA%NoCols
      write(F,10) iV(i)
    end do
    close(F)
    
    
    !Row indices
    open(unit=F,file=RowName,status='unknown')
        if(present(iLabel)) then
            write(F,'(A)') trim(iLabel)
        else 
            write(F,'(A)')
        end if
        do i=1,iA%NoRows+1
          write(F,20) iA%RowIdxs(i)
        end do
    close(F)
        
        
    
    !Col indices
    open(unit=F,file=ColName,status='unknown')
    if(present(iLabel)) then
        write(F,'(A)') trim(iLabel)
    else 
        write(F,'(A)')
    end if
    do i=1,iA%NoElements
      write(F,20) iA%ColIdxs(i)
    end do
    close(F)
    
    !Values
    open(unit=F,file=ValName,status='unknown')
    if(present(iLabel)) then
        write(F,'(A)') trim(iLabel)
    else 
        write(F,'(A)')
    end if
    do i=1,iA%NoElements
      write(F,10) iA%Vals(i)
    end do
    close(F)

    return
30  format(4(1xi<NIntegerDigits>))   
20  format(3(1xi<NIntegerDigits>))   
10  format(1x1pe<7+NScientificDigits>.<NScientificDigits>)    
    end subroutine SMatWriteSparse   
    
    
    subroutine SMatReadSparse(oA,oV, iFileName,BlockSize)
    !Gianluca Fiandaca, April 2013, AU
    !This routine reads a sparse matrix iA in 3 different files, one containing
    !the matrix row indices, one containing the matrix column indices and one containing 
    !the matrix values. The names of the output files are built concatenating the file name 
    !iName and the words "RowIdxs.txt", "ColIdxs.txt" & "Values.txt" 
    !IO
    !  iA          -> Sparse matrix to be written
    !  iName       -> the name of the 3 output files 
    !  iLabel      -> Label to be written before the matrix
    !  BlockSize   -> The size of each block in the matrix (Used for blocksplitting)    
    !------------------------------------
    !Updated Tue & Casper, September 2014.
    !The reader can now read both CSR and dense format
    !Furthermore the reader now assumes the file starts with a label line at the very top
    !------------------------------------
    use mArrays
    implicit none
    type (TSparseMat),intent(out)          :: oA
    real*8,dimension(:),pointer            :: oV
    character*(*)    ,intent(in), optional :: iFileName
    integer, intent(out)                   :: Blocksize 
    integer i,j,k,Error
    integer,parameter               :: F=99
    real*8                          :: Val
    character*256                   :: RowFileName,ColFileName,ValFileName,S,DefFileName,RHSFileName
    integer :: NoCols,NoRows,NoElements,IsCSR,Pos
    logical :: ErrorCode
    character*256 :: str

    !names and legend
    if(present(iFileName)) then
      DefFileName=trim(iFileName)//'_definition.txt'  
      RHSFileName=trim(iFileName)//'_RHS.txt'
      RowFileName=trim(iFileName)//'_RowIdxs.txt'
      ColFileName=trim(iFileName)//'_ColIdxs.txt'
      ValFileName=trim(iFileName)//'_Values.txt'
    else
      DefFileName='Definition.txt'
      RHSFileName='RHS.txt'
      RowFileName='RowIdxs.txt'
      ColFileName='ColIdxs.txt'
      ValFileName='Values.txt'
    end if
    
    open(unit=F,file=DefFileName,status='old')
    do i=1,4
      read(F,'(A)') str ! read line
      str = adjustl(str) ! remove left blanks
      pos = scan(str,'!') 
      if(pos.gt.0) then
        str = str(1:pos-1)  ! extract the parameter value
      end if
      select case(i)    
      case(1)  
        read(str,*) NoRows
      case(2)  
        read(str,*) NoCols
      case(3)  
        read(str,*) NoElements
      case(4)  
        read(str,*) Blocksize
      end select
    end do
    close(F)
    
    call SMatCreate(oA,NoRows,NoCols,'from file')
    call mArraysReallocateReals(oA%Vals,NoElements)
    call mArraysReallocateInts(oA%ColIdxs,NoElements)
    oA%NoRows=NoRows
    oA%NoCols=NoCols
    oA%NoElements=NoElements
       
    !Row indices
    open(unit=F,file=RowFileName,status='unknown')
    read(F,'(A)') S
    do i=1,NoRows+1
      read(F,'(A)') S
      read(S,*) oA%RowIdxs(i)
    end do
    close(F)
    
    !Col indices
    open(unit=F,file=ColFileName,status='unknown')
    read(F,'(A)') S
    do i=1,NoElements
      read(F,'(A)') S
      read(S,*) oA%ColIdxs(i)
    end do
    close(F)
    
    !Vals indices
    open(unit=F,file=ValFileName,status='unknown')
    read(F,'(A)') S
    do i=1,NoElements
      read(F,'(A)') S
      read(S,*) oA%Vals(i)
    end do
    close(F)
    
    !RHS
    allocate(oV(NoCols))
    open(unit=F,file=RHSFileName,status='unknown')
    read(F,'(A)') S
    do i=1,NoCols
      read(F,'(A)') S
      read(S,*) oV(i)
    end do
    close(F)
    return

    end subroutine SMatReadSparse    

   
       subroutine SMatDestroy(ioSparseMat)
!Casper Kirkegaard, November 2009.
!This routine destroys a TSparseMat structures, ie. 
!deallocates memory etc.
!
!Updated 2014 Tue, Now it takes complex numbers as well
!
!Updated Dec 2016 Kristoffer Andersen. Now also resets IsSorted and IsNativeCSR
!
!IO
!   ioSparseMat  -> Sparse matrix to destroy    
    implicit none
    type (TSparseMat),intent(inout) :: ioSparseMat
!Implementation
    !Deallocate the coordinate form arrays(Always allocated)
    if (associated(ioSparseMat%Vals)) then
      deallocate(ioSparseMat%Vals)
      deallocate(ioSparseMat%ColIdxs)
      deallocate(ioSparseMat%RowIdxs)
    end if    
    if (associated(ioSparseMat%cVals)) then
      deallocate(ioSparseMat%cVals)
      deallocate(ioSparseMat%ColIdxs)
      deallocate(ioSparseMat%RowIdxs)
    end if  
    nullify(ioSparseMat%Vals)   
    nullify(ioSparseMat%CVals)   
    nullify(ioSparseMat%ColIdxs)
    nullify(ioSparseMat%RowIdxs)
    !Reset vars   
    ioSparseMat%NoRows=0
    ioSparseMat%NoCols=0
    ioSparseMat%NoElements=0
    ioSparseMat%IsComplex=.false.
  end subroutine SMatDestroy

    subroutine SMatGetDiag(iA,oDiag,iEnsurePositiveDefinite)
    !Tue June 2015
    !This function should be called through its interface SMatGetDiag
    !This function returns the value of the diagonal elements in a matrix saved in CSR format
    !IO
    !iA -> input sparse matrix in CSR format
    !oDiag -> output oDiag(i) contains the inverse value of the i'th diagonal element in the CSR matrix given in iA
    !iEnsurePositiveDefinite -> (default = false) can be used to check if the diagonal is positive - NOT FULL POSITIVE DEFINITENESS
    !
    !CK, october 2015
    !Changed oDiag to no longer be allocatable - it not possible to have an output array alloctable. 
    !
    ! KRA, Dec 2016
    ! Added Positive diagonal check and check of the size of oDiag compared to iA%NoRows
    !
    implicit none    
    type (TSparseMat),intent(inout)  :: iA
    real*8,intent(out)   :: oDiag(:)
    logical,intent(in),optional      :: iEnsurePositiveDefinite
    integer                          :: i,k
    integer                          :: oError
    real*8                           :: tmp
    character*256                    :: S
    logical                          :: Error = .False.
    integer                          :: OldThreads
    logical                          :: CheckPositiveDefinite
    !allocate(oDiag(iA%NoRows))
    if(size(oDiag,1).lt.iA%NoRows) goto 903
    if(present(iEnsurePositiveDefinite)) then
      CheckPositiveDefinite=iEnsurePositiveDefinite
    else 
      CheckPositiveDefinite=.False.
    end if
    
    
    if(CheckPositiveDefinite) then
      !$OMP PARALLEL IF(StartOpenMP(2,1)) DEFAULT(SHARED) PRIVATE(i,tmp,oError)
      !$OMP DO SCHEDULE(STATIC)
      do i=1,iA%NoRows
        tmp=SMatGetValue(iA,i,i,oError)
        if (tmp.gt.0d0) then
          oDiag(i)=tmp
        else 
          Error=.true.
          k=i
        end if
      end do
      !$OMP END DO
      !$OMP END PARALLEL  
    else
      !$OMP PARALLEL IF(StartOpenMP(2,1)) DEFAULT(SHARED) PRIVATE(i,tmp,oError)
      !$OMP DO SCHEDULE(STATIC)
      do i=1,iA%NoRows
        oDiag(i)=SMatGetValue(iA,i,i,oError)
      end do
      !$OMP END DO
      !$OMP END PARALLEL  
    end if        
    
    call ReturnOldThreadNum()
    if (Error) goto 901
    return
    901 continue
    print*, 'Critical failure in SMatGetDiag_real, matrix-diagonal is ',SMatGetValue(iA,k,k,oError),' at ',k
903 continue
    print*, 'Output array allocation too small in SMatGetDiag_real'
    End subroutine SMatGetDiag

     
    real*8 function SMatGetValue(iA,iRowNo,iColNo,oError)
!Casper Kirkegaard, November 2009.
!This routine locates and returns the value of a given row,col entry. 
!NOTE THAT THE MATRIX IS FINALIZED WHEN CALLING THIS PROCEDURE!!!  
!IO
!  iA     -> Sparse matrix to search for value
!  iRowNo -> Index of the row to search for
!  iColNo -> Index of the column to search for
!  oError -> -1 if value does not exist, ie. entry equal to zero
!            Positive value if the value entry exists.   
!IO parameters
    implicit none    
    type (TSparseMat),intent(inout)  :: iA
    integer,intent(in)               :: iColNo
    integer,intent(in)               :: iRowNo
    integer,intent(out)              :: oError
!Variable declarations  
    integer                          :: CorColNo     
!Implementation  
    !Apply column index correction if necessary. Also make sure to switch back the ColIdxCorrected status flag
    !at the end of this routine if this is the routine applying the correction.
    CorColNo=iColNo    
    !Return the value
    call SMatLocateCSRIdx(iA,iRowNo,CorColNo,oError)
    if (oError.eq.-1) then
      SMatGetValue=0d0
    else
      SMatGetValue=iA%Vals(oError)  
    end if  

    !Switch back the col index correction lock if necessary.
end function SMatGetValue

    subroutine SMatDAD(A,iD)
      type (TSparseMat),intent(inout)  :: A
      real*8,dimension(:),intent(in)   :: iD 
      integer :: i,j
      do i=1,A%NoRows
        do j=A%RowIdxs(i),A%RowIdxs(i+1)-1
          A%Vals(j)=A%Vals(j)*iD(i)*iD(A%colIdxs(j))
        end do
      end do
      
    end subroutine SMatDAD

       subroutine SMatLocateCSRIdx(iA,iRowNo,iColNo,oIndex,iFindNearest)
!Casper Kirkegaard, November 2009.
!This routine locates and returns the CSR index of a given row,col entry.  
!-------------------------------------
!Updated September 2014, by Tue
!Added the optional parameter iFindNearest.
!------------------------------------
!IO
!  iA     -> Sparse matrix to search
!  iRowNo -> Index of the row to search for
!  iColNo -> Index of the column to search for
!  oIndex -> Index in iA%Vals of the entry iRow,iColNo
!            Set to -1 if the entry does not exits in the matrix 
!  iFindNearest -> optional parameter, which tells whether to return the nearest match if iIndex does not exist.
!                  default is 0, which returns -1 if index is not found
!                  if set to 1, it returns the nearest index searching forward
!                  if set to -1, it returns the nearest index searching backwards
!IO parameters
    implicit none    
    type (TSparseMat),intent(inout)  :: iA
    integer,intent(in)               :: iColNo
    integer,intent(in)               :: iRowNo
    integer,intent(out)              :: oIndex
    integer,intent(in),optional      :: iFindNearest
!Variable declarations    
    integer                          :: CorColNo
    integer                          :: StartIndex,EndIndex 
    integer                          :: FindNearest
!Implementation  
    if (present(iFindNearest)) then
      FindNearest = iFindNearest
    else
      FindNearest = 0
    end if
    !Apply column index correction if necessary. Also make to switch back the ColIdxCorrected status flag
    !at the end of this routine if this is the routine applying the correction.
    CorColNo=iColNo

    !Do an n*log(n) search by exploiting that the column indices of each
    !row are sorted by increasing value. Chop up the entire interval in halfs
    !and determine in which of the intervals the index lies. Continue until the value
    !is found or the width of the search interval is 1 and the values has not been found.
    StartIndex=iA%RowIdxs(iRowNo)
    EndIndex=iA%RowIdxs(iRowNo+1)-1
    call SMatLocateIdx(iA%ColIdxs,StartIndex,EndIndex,CorColNo,oIndex,FindNearest)    
    !Switch back the col index correction lock if necessary.
   end subroutine SMatLocateCSRIdx

        subroutine SMatLocateIdx(iArray,iBegin,iEnd,iIndex,oIndex,iFindNearest)
!Casper Kirkegaard, November 2009.
!This routine locates the position of a given value in a sorted integer array.  
!-------------------------------------
!Updated September 2014, by Tue
!Added the optional parameter iFindNearest.
!------------------------------------
!IO
!  iArray -> Array to search through
!  iBegin -> Index in iArray to start from
!  iEnd   -> Index in iArray to end at
!  iIndex -> Integer to locate index in iArray for
!  oIndex -> index in iarray of the value iIndex
!            Set to -1 if the entry does not exits in the  
!  iFindNearest -> optional parameter, which tells whether to return the nearest match if iIndex does not exist.
!                  default is 0, which returns -1 if index is not found
!                  if set to 1, it returns the nearest index searching forward
!                  if set to -1, it returns the nearest index searching backwards
!IO parameters
    implicit none
    integer, dimension(:),intent(in) :: iArray
    integer,intent(in)               :: iBegin
    integer,intent(in)               :: iEnd
    integer,intent(in)               :: iIndex
    integer,intent(out)              :: oIndex
    integer,intent(in),optional      :: iFindNearest
!Variable declarations    
    integer                          :: i,j   
    integer                          :: StartIndex
    integer                          :: MidIndex
    integer                          :: EndIndex   
    real*8                           :: Dummy 
    integer                          :: BruteIndex
    logical                          :: Searching
    integer                          :: Count
    integer                          :: SearchForNearest
!Implementation  
    !BruteIndex=-1
    !do i=iA%RowIdxs(iRowNo),iA%RowIdxs(iRowNo+1)-1
    !  if (iA%ColIdxs(i).eq.iColNo) then
    !    BruteIndex=i
    !  end if 
    !end do
    !Do a log(n) search by exploiting that the column indices of each
    !row are sorted by increasing value. Chop up the entire interval in halfs
    !and determine in which of the intervals the index lies. Continue until the value
    !is found or the width of the search interval is 1 and the values has not been found.
    if (present(iFindNearest)) then 
        SearchForNearest = iFindNearest
    else
        SearchForNearest = 0
    end if
    Searching=.true.
    oIndex=-1
    StartIndex=iBegin
    EndIndex=iEnd
    !Don't search if the row is empty ...
    if (EndIndex.lt.StartIndex) Searching=.False.
    Count=1
    do while(Searching)
      Count=Count+1
      if (EndIndex-StartIndex.lt.2) then
        Searching=.false.
        do i=StartIndex,EndIndex
          if (iArray(i).eq.iIndex) then
            oIndex=i
            exit
          end if   
        end do
      end if
      Dummy=(StartIndex+EndIndex)/2
      MidIndex=AnInt(Dummy)  
  
      if (iArray(MidIndex).eq.iIndex) then
        Searching=.false.
        oIndex=MidIndex 
      else
        if (iArray(MidIndex).gt.iIndex) then
          EndIndex=MidIndex
        else
          StartIndex=MidIndex
        end if  
      end if      
    end do
    !If the exact index was not found, do we look for nearest index?
    if ((oIndex.eq.-1).and.(SearchForNearest.ne.0)) then
      if (SearchForNearest.eq.1) then
        if (iArray(StartIndex).gt.iIndex) then
          oIndex=StartIndex
        else
          oIndex=EndIndex
        end if
      else if(SearchForNearest.eq.-1) then
        if (iArray(EndIndex).lt.iIndex) then
          oIndex=EndIndex
        else
          oIndex=StartIndex
        end if
      end if
    end if
    end subroutine SMatLocateIdx
   
   subroutine SMatGetDiagonalDominant(iA,oDiagDom)
   !Tue Dec, 2014
   !Determines the 
   !
   !IO
   !iA -> matrix in CSR format
   !oDiagDom -> on output the diagonal dominance of the matrix, if DiagDom<1 it is strictly diagonal dominant   
    implicit none    
    type (TSparseMat),intent(inout)     :: iA   
    real*8,intent(out)                  :: oDiagDom
    !Local vars
    integer                           :: i, j,k
    integer,allocatable,dimension(:)  :: DiagIdxs
    real*8                            :: Diag,RowVal,MinVal,Maxval    
    character*256                     :: S
    logical                           :: Error =.False.
    integer                            :: OldThreads
    allocate(DiagIdxs(iA%NoRows))
    Call SMatGetDiagIdx(iA,DiagIdxs)
    oDiagDom=0
!$OMP PARALLEL IF(StartOpenMP(2,1)) PRIVATE(I,J,Diag,RowVal) SHARED(DiagIdxs,iA,k,error) DEFAULT(NONE)  REDUCTION(+:oDiagDom)    
!$OMP DO SCHEDULE(STATIC)
      do i=1,iA%NoRows
        Diag=(abs(iA%Vals(DiagIdxs(i))))
        if(Diag.eq.0) then 
            Error=.true.
            k=i
        end if        
        RowVal=0d0
        do j=iA%RowIdxs(i),iA%RowIdxs(i+1)-1
            RowVal=RowVal+abs(iA%Vals(j))
        end do
        if (Error) then
          !Set oDiagDom to something very large if the diagelement is 0.
          oDiagDom=50d0
        else
          oDiagDom=oDiagDom+(RowVal-Diag)/Diag
        end if  
    end do
!$OMP END DO
!$OMP END PARALLEL
   oDiagDom=oDiagDom/((iA%NoRows-2)*0.5d0)
    call ReturnOldThreadNum()
    deallocate(DiagIdxs)
    return
   end subroutine SMatGetDiagonalDominant

   
     subroutine SMatAmultV(A,iV,oV)
!Casper Kirkegaard, November 2009.
!This routine calculates the product of the input matrix A with a vector iV
!IO
! A  -> Sparse A matrix
! iV -> Dense vector to multiply with
! oV -> Product A*iV 
!IO parameters
    implicit none
    type (TSparseMat),intent(inout)  :: A
    real*8,dimension(:),intent(in)   :: iV           
    real*8,dimension(:),intent(out)  :: oV
!Variable declarations    
    integer,Dimension(:),allocatable :: TempVals 
    integer                          :: i,j
    integer                          :: AFirst,ALast
    integer                          :: RowStart
    integer                          :: RowEnd    
    real*8                           :: Sum   
    integer                          :: OldThreads
    logical                          :: Dummy
!Implementation  
    dummy = StartOpenMP(2,1)
    call amux (A%NoRows, iV, oV, A%Vals,A%ColIdxs,A%RowIdxs) 
    call ReturnOldThreadNum()
     end subroutine SMatAmultV
     
    subroutine amux (n, x, y, a,ja,ia) 
      real*8  x(*), y(*), a(*) 
      integer n, ja(*), ia(*)
!c-----------------------------------------------------------------------
!c         A times a vector
!c----------------------------------------------------------------------- 
!c multiplies a matrix by a vector using the dot product form
!c Matrix A is stored in compressed sparse row storage.
!c
!c on entry:
!c----------
!c n     = row dimension of A
!c x     = real array of length equal to the column dimension of
!c         the A matrix.
!c a, ja,
!c    ia = input matrix in compressed sparse row format.
!c
!c on return:
!c-----------
!c y     = real array of length n, containing the product y=Ax
!c
!c-----------------------------------------------------------------------
!c local variables
!c
      real*8 t
      integer i, k
!c-----------------------------------------------------------------------
      !call mkl_dcsrgemv('n', n, a, ia, ja, x, y)
      !return
!$OMP PARALLEL SHARED(n, x, y, a,ja,ia) PRIVATE(t,i,k)
!$OMP DO
      do i = 1,n
        t = 0
        do k=ia(i), ia(i+1)-1 
          t = t + a(k)*x(ja(k))
        end do
        y(i) = t
      end do  
!$OMP END DO
!$OMP END PARALLEL      
!c---------end-of-amux---------------------------------------------------
!c-----------------------------------------------------------------------
    end
    
          subroutine smatgetUpperTriangular (ioA)
!   Kristoffer Andersen, Dec 2016
!   extracts the upper trangular part of ioA. 
      type(tsparsemat),intent(inout) :: ioA
      !real*8 a(*), ao(*), t
      !integer ia(*), ja(*), iao(*), jao(*)
!----------------------------------------------------------------------- 
! On entry
!-----------
! nrow  = dimension of the matrix a.
! a, ja, 
!    ia = matrix stored in compressed row sparse format
!
! nzmax = length of arrays ao,  and jao. 
!
! On return:
!----------- 
! ao, jao, 
!     iao = lower part of input matrix (a,ja,ia) stored in compressed sparse 
!          row format format.
!  
! ierr   = integer error indicator. 
!          ierr .eq. 0  means normal return
!          ierr .eq. i  means that the code has stopped when processing
!          row number i, because there is not enough space in ao, jao
!          (according to the value of nzmax) 
!
!----------------------------------------------------------------------- 
      !var
      integer :: ierr,ko,i,kold,kdiag,k,RowIDold
      !prog
      ko = 0
      do  i=1, ioA%NoRows
         RowIdOld = ioA%RowIdxs(i)
         ioA%RowIdxs(i) = ko+1
         do k = RowIdOld, ioA%RowIdxs(i+1) -1
            if (ioA%ColIdxs(k)  .ge. i) then ! upper trangle
                ko = ko+1 ! next new index
                ioA%Vals(ko) = ioA%Vals(k) ! store value
                ioA%ColIdxs(ko) = ioA%ColIdxs(k)
            endif
         enddo
     
      enddo
      
      ioA%RowIDxs(ioA%NoRows+1) = ko+1
      ioA%NoElements = ko
      return
      end subroutine smatgetUpperTriangular 
      

    
     subroutine SMatGetInvDiag(iA,oInvDiag)
    !Tue September 2014
    !This function should be called through its interface SMatGetInvDiag
    !This function returns the inverse value of the diagonal elements in a matrix saved in CSR format
    !IO
    !iA -> input sparse matrix in CSR format
    !oInvDiag -> output oInvDiag(i) contains the inverse value of the i'th diagonal element in the CSR matrix given in iA
    !
    !CK, october 2015
    !Removed allocatable property on oInvDiag. Cannot allocate and pass an allocatable outside a subroutine.
    implicit none    
    type (TSparseMat),intent(inout)  :: iA
    real*8,intent(out)               :: oInvDiag(:)
    !var
    integer                          :: i,k
    integer                          :: oError
    real*8                           :: tmp
    character*256                    :: S
    logical                          :: Error
    integer                          :: OldThreads
    Error = .False.
    !allocate(oInvDiag(iA%NoRows))
    !$OMP PARALLEL IF(StartOpenMP(2,1)) DEFAULT(SHARED) PRIVATE(i,tmp,oError)
    !$OMP DO SCHEDULE(STATIC)
    do i=1,iA%NoRows
        tmp=SMatGetValue(iA,i,i,oError)
        if (abs(tmp).gt.0d0) then
            oInvDiag(i)=1/tmp
        else 
!!$OMP CRITICAL (Errorfound)  ! might be neccesary
            Error=.true.
            k=i
!!$OMP END CRITICAL (Errorfound)
        end if
    end do
    !$OMP END DO
    !$OMP END PARALLEL  
    call ReturnOldThreadNum()
    return
    End subroutine SMatGetInvDiag
    
    double precision function ddot(n,dx,incx,dy,incy)
    double precision dx(*),dy(*),dtemp
    integer i,incx,incy,ix,iy,m,mp1,n
    ddot = 0.0d0
    dtemp = 0.0d0
    if(n<=0)return
    if(incx==1.AND.incy==1)go to 20
    ix = 1
    iy = 1
    if(incx<0)ix = (-n+1)*incx + 1
    if(incy<0)iy = (-n+1)*incy + 1
    do 10 i = 1,n
      dtemp = dtemp + dx(ix)*dy(iy)
      ix = ix + incx
      iy = iy + incy
 10 continue
    ddot = dtemp
    return

 20 m = mod(n,5)
    if( m == 0 ) go to 40
    do 30 i = 1,m
      dtemp = dtemp + dx(i)*dy(i)
 30 continue
    if( n < 5 ) go to 60
 40 mp1 = m + 1
    do 50 i = mp1,n,5
      dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
 50 continue
 60 ddot = dtemp
    return
    end
    
     subroutine SMatGetDiagIdx(iA,oDiagIdxs)
    !Tue September 2014
    !This function should be called through its interface SMatGetDiagIdx
    !This function returns the indices for the diagonal elements in a matrix saved in CSR format
    !CK, october 2015
    !Changed oDiagIdxs to no longer be allocatable - it not possible to have an output array alloctable. 
    !IO
    !iA -> input sparse matrix in CSR format
    !SMatGetDiagIdx_type -> output SMatGetDiagIdx_type(i) contains the indices for the i'th diagonal element in the CSR matrix given in iA
    implicit none    
    type (TSparseMat),intent(inout)  :: iA
    integer,intent(out)              :: oDiagIdxs(:)
    integer                          :: i
    integer                            :: OldThreads
    
    !$OMP PARALLEL IF(StartOpenMP(2,1)) DEFAULT(SHARED) PRIVATE(i)
    !$OMP DO SCHEDULE(STATIC)
    do i=1,iA%NoRows
      call SMatLocateCSRIdx(iA,i,i,oDiagIdxs(i))
    end do
    !$OMP END DO
    !$OMP END PARALLEL  
    call ReturnOldThreadNum()
    End subroutine SMatGetDiagIdx
    
      subroutine SMatGetBlockRowIdx(iA,iBlocks,oBRowIdxs)
    !Tue September 2014
    !This routine finds the first and last element in each row which are within the blocks given in iBlocks.
    !In other words, this routine is used to sort out the elements which are not used in a blocksplitting of the matrix.
    !CK, October 2015
    !Removed allocation of browidxs - the allocation hasto be performed from outside not to cause data corruption.
    !IO
    !iA -> input sparse matrix in CSR format
    !iBlocks -> The blocksplitting we wish to use on iA.
    !oBRowIdxs -> oBRowIdxs(i,1) is the first indices to first element in row i within the block, oBRowIdxs(i,2) is the last element.
    implicit none    
    type (TSparseMat),intent(inout)  :: iA
    integer,dimension(:,:),intent(in):: iBlocks
    integer,            intent(out)  :: oBRowIdxs(:,:)
    integer                          :: i,k
    integer                          :: threads
    integer                            :: OldThreads
    !allocate(oBRowIdxs(iA%NoRows,2))
    threads = size(iBlocks,1)    
    do k=1,threads
    !$OMP PARALLEL IF(StartOpenMP(2,1)) DEFAULT(SHARED) PRIVATE(i)
    !$OMP DO SCHEDULE(STATIC)
      do i = iBlocks(k,1),iBlocks(k,2)
        call SMatLocateCSRIdx(iA,i,iBlocks(k,1),oBRowIdxs(i,1),1)
        call SMatLocateCSRIdx(iA,i,iBlocks(k,2),oBRowIdxs(i,2),-1)
      end do
    !$OMP END DO
    !$OMP END PARALLEL  
    end do
    call ReturnOldThreadNum()
      End subroutine SMatGetBlockRowIdx
    
      
    function SMatFindMaxColDif(iA)
    !Tue September 2014
    !This function finds the maximum distance two elements can be apart in the same row in the matrix
    !IO
    !iA -> input sparse matrix in CSR format
    !FindMaxColDif -> output FindMaxColDif contains the maximum difference between the first and last element in each row across all rows.
    implicit none    
    type (TSparseMat),intent(in)     :: iA
    integer                          :: SMatFindMaxColDif
    integer                          :: i,j
    SMatFindMaxColDif=0
    do i=1,iA%NoRows
        SMatFindMaxColDif=max(SMatFindMaxColDif,(iA%ColIdxs(iA%RowIdxs(i+1)-1)-iA%ColIdxs(iA%RowIdxs(i))+1))
    end do
    end function SMatFindMaxColDif
   
	
   
end module mSMat      

