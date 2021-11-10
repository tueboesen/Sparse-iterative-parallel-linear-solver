! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd
!
!Updated 2016 Tue
!The RealQsort now uses real*8 instead of real.
!
! Dec 2016 KRA, added integer version
!

module mQuickSort

implicit none
public :: RealQsortC
private :: RealPartition
public :: IntQsortC
private :: IntPartition

interface  IntPartition
  module procedure IntPartition_real,IntPartition_complex,IntPartition_int
end interface
interface IntQsortC
  module procedure IntQsortC_complex, IntQsortC_real, IntQsortC_int
end interface

contains
recursive subroutine RealQsortC(A)
  real*8, intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call RealPartition(A, iq)
     call RealQsortC(A(:iq-1))
     call RealQsortC(A(iq:))
  endif
end subroutine RealQsortC

subroutine RealPartition(A, marker)
  real*8, intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  real*8 :: temp
  real*8 :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine RealPartition

recursive subroutine IntQsortC_int(SortInts,Secondarys)
  integer, intent(inout), dimension(:) :: SortInts
  integer, intent(inout), dimension(:) :: Secondarys
  integer                              :: iq

  if(size(SortInts) > 1) then
     call IntPartition(SortInts,Secondarys, iq)
     call IntQsortC(SortInts(:iq-1),Secondarys(:iq-1))
     call IntQsortC(SortInts(iq:),Secondarys(iq:))
  endif
end subroutine IntQsortC_int
recursive subroutine IntQsortC_real(SortInts,Secondarys)
  integer, intent(inout), dimension(:) :: SortInts
  real*8, intent(inout), dimension(:) :: Secondarys
  integer                              :: iq

  if(size(SortInts) > 1) then
     call IntPartition(SortInts,Secondarys, iq)
     call IntQsortC(SortInts(:iq-1),Secondarys(:iq-1))
     call IntQsortC(SortInts(iq:),Secondarys(iq:))
  endif
end subroutine IntQsortC_real


recursive subroutine IntQsortC_complex(SortInts,Secondarys)
  integer, intent(inout), dimension(:) :: SortInts
  complex*16, intent(inout), dimension(:) :: Secondarys
  integer                              :: iq

  if(size(SortInts) > 1) then
     call IntPartition(SortInts,Secondarys, iq)
     call IntQsortC(SortInts(:iq-1),Secondarys(:iq-1))
     call IntQsortC(SortInts(iq:),Secondarys(iq:))
  endif
end subroutine IntQsortC_complex

subroutine IntPartition_int(SortInts, Secondarys,marker)
  integer, intent(inout), dimension(:) :: SortInts
  integer, intent(inout), dimension(:)  :: Secondarys
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
  real*8  :: rtemp
  integer :: x      ! pivot point
  x = SortInts(1)
  i= 0
  j= size(SortInts) + 1
  do
     j = j-1
     do
        if (SortInts(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (SortInts(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = SortInts(i)
        SortInts(i) = SortInts(j)
        SortInts(j) = temp
        rtemp=Secondarys(i)
        Secondarys(i)=Secondarys(j)
        Secondarys(j)=rtemp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine IntPartition_int
subroutine IntPartition_real(SortInts, Secondarys,marker)
  integer, intent(inout), dimension(:) :: SortInts
  real*8, intent(inout), dimension(:)  :: Secondarys
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
  real*8  :: rtemp
  integer :: x      ! pivot point
  x = SortInts(1)
  i= 0
  j= size(SortInts) + 1
  do
     j = j-1
     do
        if (SortInts(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (SortInts(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = SortInts(i)
        SortInts(i) = SortInts(j)
        SortInts(j) = temp
        rtemp=Secondarys(i)
        Secondarys(i)=Secondarys(j)
        Secondarys(j)=rtemp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine IntPartition_real


subroutine IntPartition_complex(SortInts, Secondarys,marker)
  integer, intent(inout), dimension(:) :: SortInts
  complex*16, intent(inout), dimension(:)  :: Secondarys
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
  complex*16  :: rtemp
  integer :: x      ! pivot point
  x = SortInts(1)
  i= 0
  j= size(SortInts) + 1
  do
     j = j-1
     do
        if (SortInts(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (SortInts(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = SortInts(i)
        SortInts(i) = SortInts(j)
        SortInts(j) = temp
        rtemp=Secondarys(i)
        Secondarys(i)=Secondarys(j)
        Secondarys(j)=rtemp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine IntPartition_complex

end module mQuickSort



