module mInput
    implicit none 
    !Loads the settings used by the parallelization framework and the iterative solver      
  type Tsettings
    integer :: AffinityMode
    integer :: AffinityPattern
    integer :: NCPUs
    integer :: UseNested
    integer :: StartThread
    integer :: NCPUsLow
    integer :: NCPUOuter
    integer :: NCPUInner
    integer :: NCPUInnerLow
    integer :: SolverType
    integer :: MaxIter
    integer :: BlockSize
    integer :: UseNormalization
    real*8  :: FillingFactor
  end type Tsettings
    contains
    
    subroutine LoadSettings(FName,Settings)
    character*256,intent(in)           :: Fname 
    type(Tsettings),intent(inout)      :: Settings
    integer :: Nsettings,i,j,f,pos    
    character*256 :: str,str2
    F=55
    Nsettings=13
    open(unit=F,file=FName,status='old')
    do i=1,Nsettings  
      read(F,'(A)') str ! read line
      str = adjustl(str) ! remove left blanks
      pos = scan(str,'!') 
      if(pos.gt.0) then
        str = str(1:pos-1)  ! extract the parameter value
      end if
      select case(i)    
      case(2)  
        read(str,*) Settings%AffinityMode 
      case(3)  
        read(str,*) Settings%AffinityPattern
      case(4)  
        read(str,*) Settings%UseNested
      case(5)  
        read(str,*) Settings%StartThread 
      case(6)  
        read(str,*) Settings%NCPUs 
      case(7)  
        read(str,*) Settings%NCPUsLow 
      case(8)  
        read(str,*) Settings%NCPUOuter 
      case(9)  
        read(str,*) Settings%SolverType 
      case(10)  
        read(str,*) Settings%MaxIter 
      case(11)  
        read(str,*) Settings%BlockSize 
      case(12)  
        read(str,*) Settings%FillingFactor 
      case(13)  
        read(str,*) Settings%UseNormalization 
      end select
    end do
    close(F)  
    end subroutine LoadSettings
    
    end module mInput