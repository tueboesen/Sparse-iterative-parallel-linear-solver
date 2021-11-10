Module mMisc
  !Tue a module for small usefull routines that doesn't fit in anywhere else, but is used through different modules
  implicit none
  interface LoadOptionalParam
    module procedure LoadOptionalParam_int, LoadOptionalParam_real, LoadOptionalParam_real2
  end interface
  contains

  Function LoadOptionalParam_int(DefaultValue,iParam)
    !Tue Dec, 2014
    !This function is used to load optional parameters
    !IO
    !DefaultValue -> A default value which we load if the optional parameter does not exist.
    !iParam       -> the optional parameter that we want to load
    implicit none
    integer, intent(in)          :: DefaultValue
    integer, intent(in),optional :: iParam
    integer                      :: LoadOptionalParam_int
    if(present(iParam)) then
        LoadOptionalParam_int=iParam
      else 
        LoadOptionalParam_int=DefaultValue
      end if
  end function LoadOptionalParam_int
  Function LoadOptionalParam_real(DefaultValue,iParam)
    !Tue Dec, 2014
    !This function is used to load optional parameters
    !IO
    !DefaultValue -> A default value which we load if the optional parameter does not exist.
    !iParam       -> the optional parameter that we want to load
    implicit none
    real*8, intent(in)          :: DefaultValue
    real*8, intent(in),optional :: iParam
    real*8                      :: LoadOptionalParam_real
    if(present(iParam)) then
        LoadOptionalParam_real=iParam
      else 
        LoadOptionalParam_real=DefaultValue
      end if
  end function LoadOptionalParam_real
  Function LoadOptionalParam_real2(DefaultValue,iParam)
    !Tue Dec, 2014
    !This function is used to load optional parameters
    !IO
    !DefaultValue -> A default value which we load if the optional parameter does not exist.
    !iParam       -> the optional parameter that we want to load
    implicit none
    integer, intent(in)          :: DefaultValue
    real*8, intent(in),optional :: iParam
    real*8                      :: LoadOptionalParam_real2
    if(present(iParam)) then
        LoadOptionalParam_real2=iParam
      else 
        LoadOptionalParam_real2=DefaultValue
      end if
  end function LoadOptionalParam_real2
    end module mMisc

    
    module mArrays
   interface ArraysCheckAndIncAlloc
    module procedure mArraysCheckAndIncAllocReal, mArraysCheckAndIncAllocInt
   end interface
   
  contains  
    subroutine mArraysReallocateInts(iArray,iNewSize)
    implicit none
    integer,dimension(:),pointer,intent(inout) :: iArray
    integer,intent(in)                         :: iNewSize
    integer,dimension(:),pointer               :: TempArray 
    integer                                    :: CurSize 
    integer                                    :: TransferSize 
    integer                                    :: Error
    character*256                              :: S
    !Find the size to transfer to the resized array to be allocated.
    CurSize=Size(iArray)

    if(CurSize==0) then
      !goto 902
      !CS: Below: in case we don't want an error message
      allocate(iArray(iNewSize),Stat=Error)
      if (Error.ne.0)&
        goto 900
      return
    end if  

    if (CurSize.ne.iNewsize) then
      TransferSize=CurSize
      if (TransferSize.gt.iNewSize)&
        TransferSize=iNewSize      
      !Allocate temporary arrays.
      allocate(TempArray(TransferSize),Stat=Error) 
      if (Error.ne.0)&
        goto 900
      !Transfer values to temporary arrays - first the real*8 array holding the values.    
      TempArray(1:TransferSize)=iArray(1:TransferSize)
      !Deallocate old array
      if (associated(iArray))&
        deallocate(iArray)
      !Allocate new array of updated length and transfer the old values
      allocate(iArray(iNewSize),Stat=Error) 
      if (Error.ne.0)&
        goto 900
      if (TransferSize+1.le.iNewSize) then
        iArray(TransferSize+1:iNewSize)=0
      end if
      iArray(1:TransferSize)=TempArray(1:TransferSize)  
      deallocate(TempArray)  
    end if
    return
900 continue
    S='Error re-allocating memory. Please increase the free physical memory or decrease the number of data.'
     go to 999
901 continue
    S='Error de-allocating memory. Please increase the free physical memory or decrease the number of data.'
     go to 999
902 continue
    S='Error in mArraysReallocateInts iArray array has a 0 size at the entry of the subroutine.'
     go to 999
999 print*, S    
    end subroutine mArraysReallocateInts
      
!******************************************************************************      
    subroutine mArraysReallocateReals(iArray,iNewSize)
    implicit none
    real*8,dimension(:),pointer,intent(inout) :: iArray
    integer,intent(in)                        :: iNewSize
    real*8,dimension(:),pointer               :: TempArray 
    integer                                   :: CurSize 
    integer                                   :: TransferSize 
    integer                                   :: Error
    character*256                             :: S
    !Find the size to transfer to the resized array to be allocated.
    CurSize=Size(iArray)

    if(CurSize==0) then
      !goto 902
      !CS: Below: in case we don't want an error message
      allocate(iArray(iNewSize),Stat=Error)
      iArray(:)=0.d0
      if (Error.ne.0)&
        goto 900
      return
    end if  

    if (CurSize.ne.iNewsize) then
      TransferSize=CurSize
      if (TransferSize.gt.iNewSize)&
        TransferSize=iNewSize      
      !Allocate temporary arrays.
      allocate(TempArray(TransferSize),Stat=Error) 
      if (Error.ne.0)&
        goto 900
      !Transfer values to temporary arrays - first the real*8 array holding the values.    
      TempArray(1:TransferSize)=iArray(1:TransferSize)
      !Deallocate old array
      if (associated(iArray)) then
        deallocate(iArray,Stat=Error)
        if (Error.ne.0)&
          goto 901
      end if  
      !Allocate new array of updated length and transfer the old values
      allocate(iArray(iNewSize),Stat=Error) 
      iArray(:)=0.d0
      if (Error.ne.0)&
        goto 900
      iArray(1:TransferSize)=TempArray(1:TransferSize) 
      deallocate(TempArray,Stat=Error) 
      if (Error.ne.0)&
        goto 901  
    end if
    return
900 continue
    S='Error re-allocating memory. Please increase the free physical memory or decrease the number of data.'
     go to 999
901 continue
    S='Error de-allocating memory. Please increase the free physical memory or decrease the number of data.'
     go to 999
902 continue
    S='Error in mArraysReallocateReals iArray array has a 0 size at the entry of the subroutine.'
     go to 999
999 print*, S    
    end subroutine mArraysReallocateReals

        subroutine mArraysCheckAndIncAllocReal(Nnew,ioAR)
! This routine increases the allocation of ioAR if
! Nnew is larger than the allocation of ioAR.
! The Nnew allocation will be 2*Nnew to avoid future allocations
! No data are transferred !    
    integer,intent(in) :: Nnew
    real*8,allocatable,intent(inout) :: ioAR(:)
    integer :: k
    if (allocated(ioAR)) then
        K=ubound(ioAR,dim=1)
        if (k.lt.Nnew) then
            deallocate(ioAR)  
            allocate(ioAR(2*Nnew+1))
        end if    
    else
        allocate(ioAR(2*Nnew))
    end if
    end subroutine mArraysCheckAndIncAllocReal
    
	!dir$ attributes offload : mic :: marrayscheckAnsIncAllocInt
    subroutine mArraysCheckAndIncAllocInt(Nnew,ioAR)
! This routine increases the allocation of ioAR if
! Nnew is larger than the allocation of ioAR.
! The Nnew allocation will be 2*Nnew to avoid future allocations
! No data are transferred !    
    integer,intent(in) :: Nnew
    integer,allocatable,intent(inout) :: ioAR(:)
    integer :: k
    if (allocated(ioAR)) then
        K=ubound(ioAR,dim=1)
        if (k.lt.Nnew) then
            deallocate(ioAR)  
            allocate(ioAR(2*Nnew+1))
        end if    
    else
        allocate(ioAR(2*Nnew))
    end if
    end subroutine mArraysCheckAndIncAllocInt    
end module mArrays    