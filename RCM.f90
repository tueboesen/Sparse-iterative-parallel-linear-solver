module mRCM
! By Kristoffer Andersen, Dec 2016
! Updated by Tue Boesen, August 2017
! These routines can be used to reorder a matrix.     
! The only dependecies are to mSmat, mError, mArrays and mQuicksort
!
!   The main routines are>
!       RCMSort - performs a reverse Cuthill-McKee ordering of a matrix
!       SymPerm - performs a symmetric permutation of a matrix: PAP^t
!       HalfPermOnlyTranspose - by using permutation operations this 
!               routine finds the transpose of a matrix    
    
    use mSMat
    implicit none
    
    contains
    subroutine RCMSort(ioA,oPerm)
    ! Kristoffer Andersen, Dec 2016
    ! this routine finds the reverse Cuthill-McKee ordering of matrix ioA
    ! and makes the symmetric permutation A = PAP^t
    !
    ! IO
    ! IO - ioA - on input the unpermuted matrix
    !            on output the permuted matrix
    ! O - oPerm - the permuatation vector - this can ie be used for the vectors used together with the matrix
    !
    implicit none
    type(Tsparsemat),intent(inout)  :: ioA
    integer,intent(out)             :: oPerm(:)
    ! prog
    call RCMOrdering(ioA,oPerm)
    call symperm(ioA,oPerm)
    end subroutine
    
    
    subroutine RCMOrdering(iA,oRCMind,opStartIndex)
    ! Kristoffer Andersen, Dec 2016
    ! this routine finds the reverse Cuthill-McKee ordering for a 
    ! symmetric sparse matrix iA.
    ! IO
    ! IO - iA - the unpermuted matrix
    ! O - RCMind - the permutation vector - this can ie be used for the vectors used together with the matrix
    ! I (optional) - opStartIndex - optional start index (default = 1)
    !
    use mArrays
    use mQuickSort
    implicit none
    type(Tsparsemat),intent(inout)  :: iA
    integer,intent(out)             :: oRCMind(:)
    integer,intent(in),optional     :: opStartIndex
    !var
    logical, allocatable        :: inRCMind(:)
    integer :: i,j,k,ni,ni1,nneighbours,jn,ad,itot,inArray,MinEle,n,MaxEle
    integer,allocatable :: Neighbours(:),Adj(:),R(:),ele(:)
    !prog
    
    allocate(inRCMind(iA%NoRows))
    allocate(Ele(iA%NoRows))
    allocate(R(iA%NoRows))
 
    inRCMInd(:) = .false. ! non sorted
    
    MaxEle=-1
    Ele=-1
    do j=1,iA%NoRows
       n=iA%RowIdxs(j+1)-iA%RowIdxs(j)
       MaxEle=max(MaxEle,n)
       Ele(j)=n
    end do
    
    allocate(Neighbours(MaxEle))
    allocate(Adj(MaxEle))    
        
    
    R=-1
    itot = 1
    inArray = 0
    do 
        if(itot.gt.iA%NoRows) exit ! termination
        if(inArray.ge.iA%NoRows ) exit
        
        if(R(itot).eq.-1) then
          MinEle=999999
          inArray=inArray+1
          do j=1,iA%NoRows
            if (inRCMInd(j)) then
              cycle    
            end if            
            n=iA%RowIdxs(j+1)-iA%RowIdxs(j)  
            if (n.lt.MinEle) then
              MinEle=n
              R(itot)=j
            end if
          end do
        end if
        i = R(itot)  !present index
        inRCMInd(i) = .true. !Mark it 
        
         
        
        ! get neighbours
        ni = iA%Rowidxs(I)
        ni1 = iA%Rowidxs(I+1)
        Nneighbours = ni1-ni
        Neighbours(1:Nneighbours) = ia%Colidxs(ni:(ni1-1))
        ! remove already present neighbours - this will also remove the self reference
        jn = 1
        do j = 1,Nneighbours
            if(.not.inRCMInd(Neighbours(J))) then 
                Neighbours(Jn) = Neighbours(J)
                jn = jn+1
            endif
        enddo
        Nneighbours = jn-1
        
        Adj(1:Nneighbours)=Ele(Neighbours(1:Nneighbours))        
        ! sorting
        call IntQSortC(adj(1:Nneighbours),Neighbours(1:Nneighbours))  ! we need a real type routine for this. 
        
        !Add new nodes to the queue
        do j = 1,Nneighbours
            inArray = inArray+1
            R(inArray) = Neighbours(j) ! add neighbours
            inRCMInd(Neighbours(J)) = .true. ! mark as included
        enddo
        itot = itot+1        
    enddo
    
    ! finally reverse the order
    do i=1,(iA%NoRows)
        oRCMind(i) = R(iA%NoRows-i+1)
    enddo
    
    return
    ! error messages
    end subroutine RCMOrdering
    
    
    subroutine halfperm(iA,oB,iPerm,opColCounts)
    ! by Kristoffer Andersen, Dec 2016
    ! this routine performs the operation B = (P*A)^t, where P is er 
    ! permutation vector of A
    ! The out is sorted independent on the status of the input - nice feature!
    !
    ! IO
    ! I iA - unpermuted matrix
    ! O oB - permuted and sorted matrix
    ! I iPerm - permutation vector
    ! I opColCounts (optional) - used to speed up repeated calls of the routine. 
    implicit none
    !io
    type(Tsparsemat),intent(in) :: iA
    type(Tsparsemat),intent(out) :: oB
    integer,intent(in) :: iPerm(:)
    integer,intent(in),optional :: opColCounts(:)
    !var
    integer :: pi,i,j,jp,p,pm1,pm2,jpt
    !prog
    
    call SMatCreate(oB,ia%NoCols,iA%NoRows,'(PA)^T',.false.,iA%NoElements)
    
    ! create column structure in permuted matrix - the column count is unchanged by row permutations
    if(present(opColCounts)) then ! when halfperm is called twice in symperm, one already knows this
        oB%RowIdxs(1:iA%NoRows) = opColCounts(1:iA%NoRows)
    else
        oB%RowIdxs(:) = 0
        do i=1,iA%NoRows
            do j=iA%RowIDxs(i),iA%RowIdxs(i+1)-1
                oB%RowIdxs(iA%ColIDxs(j)) = oB%RowIdxs(iA%ColIDxs(j))+1  ! find coumn counts
            enddo
        enddo
    endif
    ! create column pointers
    pm2 = 0
    pm1 = 1
    do i=3,iA%NoCols+3
        p = pm1+ob%RowIdxs(i-2) 
        ob%RowIdxs(i-2) = pm2
        pm2 =pm1
        pm1 = p
    enddo
    
    ! fill matrix
    do i=1,iA%NoRows
        pi = iPerm(i) ! get the permuted row
        do jp=iA%RowIDxs(pi),iA%RowIdxs(pi+1)-1
            j = iA%ColIdxs(jp)
            jpt = oB%RowIdxs(j+1)
            ob%ColIdxs(jpt) = i
            oB%Vals(jpt) = iA%Vals(jp)
            oB%RowIdxs(j+1) = jpt+1
        enddo
    enddo
    oB%RowIdxs(1) = 1
    oB%NoElements = iA%NoElements
    
    return
    ! error messages

    end subroutine halfperm
    
    subroutine halfpermOnlyTranspose(iA,oB)
    ! by Kristoffer Andersen, Dec 2016
    ! this routine performs the operation B = A^t
    ! The out is sorted independent on the status of the input - nice feature!
    !
    ! IO
    ! I iA - unpermuted matrix
    ! O oB - permuted and sorted matrix
    !
    implicit none
    !io
    type(Tsparsemat),intent(in) :: iA
    type(Tsparsemat),intent(out) :: oB
    !var
    integer :: pi,i,j,jp,p,pm1,pm2,jpt
    !prog
    call SMatCreate(oB,ia%NoCols,iA%NoRows,'(A)^T',.false.,iA%NoElements)
    
    ! create column structure in permuted matrix - the column count is unchanged by row permutations
    oB%RowIdxs(:) = 0
    do i=1,iA%NoRows
        do j=iA%RowIDxs(i),iA%RowIdxs(i+1)-1
            oB%RowIdxs(iA%ColIDxs(j)) = oB%RowIdxs(iA%ColIDxs(j))+1  ! find coumn counts
        enddo
    enddo
    
    ! create column pointers
    pm2 = 0
    pm1 = 1
    do i=3,iA%NoCols+3
        p = pm1+ob%RowIdxs(i-2) 
        ob%RowIdxs(i-2) = pm2
        pm2 =pm1
        pm1 = p
    enddo
    
    ! fill matrix
    do i=1,iA%NoRows
        pi = i ! get row
        do jp=iA%RowIDxs(pi),iA%RowIdxs(pi+1)-1
            j = iA%ColIdxs(jp)
            jpt = oB%RowIdxs(j+1)
            ob%ColIdxs(jpt) = i
            oB%Vals(jpt) = iA%Vals(jp)
            oB%RowIdxs(j+1) = jpt+1
        enddo
    enddo
    oB%RowIdxs(1) = 1
    oB%NoElements = iA%NoElements
    
    return
    ! error messages

    end subroutine halfpermOnlyTranspose
    
    subroutine symperm(ioA,iPerm)
    ! by Kristoffer Andersen, Dec 2016
    ! performs the halfperm twice in order to get B = PAP^t = (P(PA)^t)^t
    !
    ! IO
    ! IO ioA - in input: unpermuted matrix
    !          on output: permuted and sorted matrix
    ! I iPerm- permutation vector
    implicit none
    !io
    type(Tsparsemat),intent(inout) :: ioA
    integer,intent(in) :: iPerm(:)
    !var
    type(tsparsemat) :: temp
    integer,allocatable :: ColCounts(:)
    integer :: i,pi,N
    !prog
    call halfperm(ioA,temp,iPerm)
    N = ioA%NoRows
    allocate(ColCounts(N))
    do i = 1,N
        pi = iPerm(i)
        ColCounts(i) = ioA%RowIdxs(pi+1)-ioA%RowIdxs(pi)
    enddo
    call halfperm(temp,ioA,iPerm,ColCounts)
    return
    end subroutine symperm
    
    subroutine PermuteVector(ioV,iPerm)
    ! by Kristoffer Andersen, Dec 2016
    ! this routine permutes a vector according to iPerm
    !
    ! IO
    ! IO ioV - in input: unpermuted vector
    !          on output: permuted vector
    ! I iPerm- permutation vector
    real*8,intent(inout) :: ioV(:)
    integer,intent(in) :: iPerm(:)
    ! var
    real*8,allocatable :: temp(:)
    integer :: i,N
    !prog
    N = ubound(ioV,1)
    allocate(temp(N))
    temp(1:N) = ioV(1:N)
    do i=1,N
        ioV(iPerm(i)) = temp(i)
    enddo
    return
    end subroutine PermuteVector
    
end module mRCM
    
    
    