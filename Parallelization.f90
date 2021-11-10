Module mParallel
  !Tue Boesen 2017  
  !This module is meant to work as a general parallelization framework, and is meant to replace all other parallelization existing in a code.
  !The general structure of the code is as follows:
  !On initialization the subroutine InitOpenMP() should be called, it is imperative that this routine is called before any OpenMP routines are called - even omp_get_wtime()
  !When creating parallel regions anywhere in the code, these regions should in general always be created with an if call to StartOpenMP 
  use omp_lib
  use IFPORT, only: SETENVQQ
  implicit none

  !Module variables - These variables changes dynamically during a run  
  logical,private :: RunningSingleParallelRegion !Used to lock against further parallelization while running a nonnested region.
  integer,private :: OldThreadNumber             !Saves the previous threadnumber so we can return to this number afterwards if desired.
  !Initial settings - These variables get set in InitOpenMP, and then remain static after that
  integer,private :: UseNested                   !Use nested parallelization?   
  integer,private :: NCPUs                       !Total number of threads
  integer,private :: NCPUsLow                    !Total number of threads in a memory bandwidth intense region  
  integer,private :: NCPUOuter                   !Number of CPUs in the outer region of a nested parallelization
  integer,private :: NCPUInner                   !Number of CPUs in the inner region of a nested parallelization
  integer,private :: NCPUInnerLow                !Number of CPUs in the inner region of a nested parallelization with memory bandwith limits
  integer,private :: AffinityMode                !Dedicated or shared server?
  integer,private :: AffinityPattern             !Scattered or compact affinity
  integer,private :: StartThread                 !ThreadID of first affinity binding
  contains

  subroutine InitOpenMP(iAffinityMode,iAffinityPattern,iStartThread,iUseNested, iNCPUs,iNCPUsLow,iNCPUOuter)
    !This is the initialization routine for OpenMP. 
    !It NEEDS to be run before any OpenMP routines are called - even omp_get_wtime().
    !This routine setup the entire parallel framework. 
    !Based on the input variables it determines the number of threads to use, the affinity which they should be used with, the parallelization structure 
    !On output this routine will have set all the following private module variables: NCPUs, NCPUInner, NCPUOuter, NCPUInnerLow, NCPUsLow, AffinityMode, UseNested, StartThread.
    !All these variables will be used by the subroutine StartOpenMP, which should be called everytime a parallel region is potentially about to be created.
    
    !IO
      !Note any variable can be set to -99 for automatic determination. 
      !AffinityMode    -> Specifies different modes thread binding can occur. 
      !  0=shared mode. 
      !                 Threads bind to a socket, but are free to migrate between the cores of a socket.
      !                 This mode is usefull when running on a server with multiple users, however,
      !                 performance can also be lower than the dedicated mode. Shared mode does not use nested parallelization in any way.
      !                 Note if all CPU's are used, AffinityMode will automatically be changed from shared to dedicated
      !  1=Dedicated mode. 
      !                 Threads bind to individual cores, and hence are locked in place. This gives better performance, but only if the cores are not used for anything else.  
      !  2=manual mode. No affinity settings are made inside the program. It is left to the user to set proper environment variables, before starting
      !                 the program.
      !iAffinityPattern -> Sets the affinity pattern for the threadbinding. This is mostly relevant on NUMA systems, or when using nested parallelization.
      !  0=Scatter threads.
      !                 Threads are distributed as far apart as possible.
      !  1=Compact threads.
      !                 Threads are clustered together.  
      !iStartThread    -> The threadID the first thread is bound to. 
      !iUseNested      -> Determines whether to enable nested parallelization or not                   (0=no nested parallelization allowed,1=Nested parallelization enabled)
      !iNCPUs          -> Total number of CPUs to use                                                  (0=all, -1=all-1, -2=all-2, ect.)  
      !iNCPUsLow       -> Total number of CPUs to use in a bandwidth intensive regions                 (0=all, -1=all-1, -2=all-2, ect.)  
    !The following are only relevant if nested parallelization is used.                                
      !iNCPUOuter      -> Number of CPUs to use in outer parallel region.                              (Should be set between 2..and Number of numanodes)                               
    use ifwin
    implicit none
    integer,intent(in) :: iAffinityMode,iAffinityPattern,iStartThread,iUseNested,iNCPUs,iNCPUsLow,iNCPUOuter
    !Local vars
    logical(4)           :: Success
    integer              :: i,j
    integer              :: numaNodeCount,processorPackageCount,processorCoreCount,logicalProcessorCount
    integer,dimension(3) :: processorCacheCount
    real*8               :: rand
    logical              :: Debug = .false. !Only used for debugging, to see how the threads are bound and what was detected.
    character(1024)      :: KmpAff          
    !Vars for manual binding
    integer(HANDLE) Process
    integer(DWORD_PTR) AffinityMask
    integer Core
    integer(BOOL) Retval
    integer(DWORD) iError
    
      call GetMachineArch(numaNodeCount,processorPackageCount,processorCoreCount, logicalProcessorCount,processorCacheCount ) !Get machine architecture
      if(Debug) then
        print*, "GetLogicalProcessorInformation results:"
        print*, "  Number of NUMA nodes: ", numaNodeCount
        print*, "  Number of physical processor packages: ", processorPackageCount
        print*, "  Number of processor cores: ", processorCoreCount
        print*, "  Number of logical processors: ", logicalProcessorCount
        print*, "  Number of processor L1/L2/L3 caches: ",processorCacheCount    
      end if
      
      !We set StartThread
      StartThread = iStartThread
      if (StartThread.le.-1) then
         call Random_seed()
         call Random_Number(rand) !Random number between 0-1
         StartThread=floor(rand*logicalProcessorCount) !scaled to 0-(logicalprocessorcount-1)
      end if   
      StartThread=max(0,StartThread) !Sanity check
      StartThread=min(logicalProcessorCount-1,StartThread) !Sanity check
      
      !We set NCPUs
      NCPUs=iNCPUs
      if (NCPUs.eq.-99) then
        NCPUs=processorCoreCount/2 !We use approximately half the cores if nothing else is specified.   
      end if      
      if (NCpus.le.0) then
        NCPUs = max(processorCoreCount+NCpus,1)
      end if
      NCPUs=min(logicalProcessorCount,NCPUs) !Sanity check
      NCPUs=max(1,NCPUs) !Sanity check
      
      !We set NCPUsLow
      NCPUsLow=iNCPUsLow
      if (NCPUsLow.eq.-99) then
        NCPUsLow=NCPUs*2/3 !We use approximately 2/3 of the threads specified for a bandwith intensive region, if nothing else is specified.   
      end if      
      if (NCpusLow.le.0) then
        NCPUsLow = max(processorCoreCount+NCpusLow,1)
      end if
      NCPUsLow=min(NCPUs,NCPUsLow) !Sanity check
      NCPUsLow=max(1,NCPUsLow) !Sanity check
            
      !We set AffinityMode
      AffinityMode=iAffinityMode
      if (AffinityMode.eq.-99) then
        AffinityMode=0
      end if
      if(NCpus.ge.logicalProcessorCount.AND.AffinityMode.EQ.0) then !We use all cores, no reason to run shared mode in that case.
        AffinityMode=1 
      end if
      if (AffinityMode.lt.0.OR.AffinityMode.gt.2) then !Sanity check
        print*, 'Warning! Affinity mode not set correctly. Setting shared AffinityMode'        
        AffinityMode=0
      end if
      
      !We set UseNested
      UseNested=iUseNested
      if ((AffinityMode.eq.0).OR.(NCPUs.lt.2).OR.(UseNested.ne.1)) then
        UseNested=0
      end if
      
      !We set AffinityPattern
      AffinityPattern=iAffinityPattern
      if((AffinityPattern.ne.0).OR.(AffinityPattern.ne.1)) then
        if(UseNested) then
          AffinityPattern=0 !Scatter
        else
          AffinityPattern=1 !Compact
        end if
      end if
      
      if(UseNested.eq.1) then
        !We set NCPUOuter
        if (NCPUOuter.eq.-99) then
          NCPUOuter=numaNodeCount 
        end if
        NCPUOuter=min(max(NCPUOuter,2),numaNodeCount) !Sanity check
        NCPUOuter=min(NCPUOuter,NCPUs) !Sanity check
        
        !We set NCPUInner
        NCPUInner=NCPUs/NCPUOuter 
 
        !We set NCPUInnerLow
        NCPUInnerLow=NCPUsLow/NCPUOuter 
        
        if(NCPUs.ne.NCPUOuter*NCPUInner) then
          !The original CPU numbers do not match a multiplum of the detected NumaNodes
          NCPUs=NCPUOuter*NCPUInner
          if(iNCPUs.ne.-99) then
            print*,'NCPUs was originally chosen as',iNCPUs
            print*,'In order to preserve even load across the NUMAnodes NCPUs is changed to:', NCPUs
          end if
        end if
        if(NCPUsLow.ne.NCPUOuter*NCPUInnerLow) then
          !The original Low CPU numbers do not match a multiplum of the detected NumaNodes
          NCPUsLow=NCPUOuter*NCPUInnerLow
          if(iNCPUsLow.ne.-99) then
            print*,'NCPUsLow was originally chosen as',iNCPUsLow
            print*,'In order to preserve even load across the NUMAnodes NCPUs is changed to:', NCPUsLow
          end if
        end if
      else
        NCPUOuter = NCPUs
        NCPUInner = 1
      end if
      
      !All parameters determined and set, now we set all the actual environment variables
      if(NCPUs.ge.1.AND.AffinityMode.ne.2) then       
        call GetKmpAffinityStr(KmpAff,NCPUs,numaNodeCount,processorCoreCount, logicalProcessorCount,AffinityMode,AffinityPattern,StartThread)
        Success = SETENVQQ(KmpAff) 
      else
        Process = GetCurrentProcess()
        Core = StartThread 
        AffinityMask = ishft(1_DWORD_PTR, Core)
        Retval = SetProcessAffinityMask(Process,AffinityMask)
        if(Retval == 0) then
          iError = GetLastError()
          write(*,*) 'Failed to set the affinity properly',iError
          !stop
        end if  
      end if
      if(Debug) then
        print*,'omp_bind_string:',KmpAff
        print*,'NCPU',NCPUs
        print*,'NCPUsLow',NCPUsLow
        print*,'NCPUouter',NCPUouter
        print*,'NCPUInner',NCPUinner
        print*,'NCPUInnerLow',NCPUinnerLow
        print*,'UseNested',UseNested
        print*,'AffinityMode',AffinityMode
        print*,'AffinityPattern',AffinityPattern
        print*,'StartThread',StartThread
      end if
      call KMP_SET_STACKSIZE_S(8000000) !In our code this is essential otherwise we get very subtle stack memory corruption in some of our routines
      if (UseNested) then
        call omp_set_nested(1) 
        call omp_set_max_active_levels(2)        
      else 
        call omp_set_nested(0)
      end if
      call omp_set_num_threads(NCPUs)     
      !OpenMP settings are applied when the first parrallel loop is found. We do a dummy loop here to get it done now..
      !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(NCpus) REDUCTION(+:J)
        J=0
      !$OMP DO
        do i=1,NCpus
          J=J+1      
        end do   
      !$OMP END DO
      !$OMP END PARALLEL
      if(debug) then
        print*,'InitOpenMP complete'
      end if            
      return
  end subroutine InitOpenMP

  subroutine GetMachineArch(numaNodeCount,processorPackageCount,processorCoreCount, logicalProcessorCount,processorCacheCount )
    !This routine gets all the current hardware system information.
    !IO
    !numaNodeCount         -> Number of NUMA nodes
    !processorPackageCount -> Number of physical processor packages
    !processorCoreCount    -> Number of processor cores
    !logicalProcessorCount -> Number of logical processors
    !processorCacheCount   -> Number of processor L1/L2/L3 caches
    use, intrinsic :: ISO_C_BINDING
    use kernel32
    implicit none   
    integer,intent(out)              :: numaNodeCount,processorPackageCount,processorCoreCount,logicalProcessorCount
    integer,dimension(3),intent(out) :: processorCacheCount     
    !Local vars
    procedure(GetLogicalProcessorInformation), pointer :: glpi
    type(T_SYSTEM_LOGICAL_PROCESSOR_INFORMATION), allocatable, dimension(:) :: buffer
    integer(DWORD) :: returnLength = 0
    integer(DWORD) :: ret
    integer :: nlpi, lpi_element_length, i      
      processorCacheCount=0
      numaNodeCount = 0
      processorCoreCount = 0
      logicalProcessorCount = 0
      processorPackageCount = 0
      call c_f_procpointer( & !Get kernel information
        transfer( &
        GetProcAddress( &
        GetModuleHandle("kernel32"//C_NULL_CHAR), &
        "GetLogicalProcessorInformation"//C_NULL_CHAR &
        ), &
        C_NULL_FUNPTR &
        ), &
        glpi)
    
      if (.not. associated(glpi)) then
        print*,  "GetLogicalProcessorInformation not supported"
        error stop
      end if
      !We don't know in advance the size of the buffer we need. We'll pick a number, allocate it,
      !and see if that's sufficient.  If not, we'll use the returned size information and reallocate
      !the buffer to the required size.
      allocate (buffer(20))
      lpi_element_length = C_SIZEOF(buffer(1))
      returnLength = C_SIZEOF(buffer)
      ret = glpi(buffer, returnLength)
      if (ret == FALSE) then ! Failed
        if (GetLastError() == ERROR_INSUFFICIENT_BUFFER) then
          deallocate (buffer)
          allocate (buffer(returnLength/lpi_element_length))
          ret = glpi(buffer, returnLength)
          if (ret == FALSE) then
            print*,  "GetLogicalProcessorInformation call failed with error code ", GetLastError()
            error stop
          end if
        else
          print*, "GetLogicalProcessorInformation call failed with error code ", GetLastError()
          error stop
        end if
      end if
      !Now we can iterate through the elements of buffer and see what we can see
      do i=1, returnLength / lpi_element_length ! Number of elements in buffer
        select case (buffer(i)%Relationship)
        case(RelationNumaNode)
          numaNodeCount = numaNodeCount + 1
        case(RelationProcessorCore)
          processorCoreCount = processorCoreCount + 1
          !A Hyperthreaded core supplies more than one logical processor
          logicalProcessorCount = logicalProcessorCount + popcnt(buffer(i)%processorMask)
        case(RelationCache)
          !One cache descriptor for each cache
          if (buffer(i)%Cache%Level > 0 .and. buffer(i)%Cache%Level <= 3) then
            processorCacheCount(buffer(i)%Cache%Level) = processorCacheCount(buffer(i)%Cache%Level) + 1
          else
            print*, "Invalid processor cache level ", buffer(i)%Cache%Level
          end if
        case(RelationProcessorPackage)
          !Logical processors share a physical package (socket)
          processorPackageCount = processorPackageCount + 1
        case default
          print*, "Unrecognized relationship code ", buffer(i)%Relationship
        end select
      end do
      if(allocated(buffer)) deallocate(buffer)
  end subroutine GetMachineArch
    
    
  subroutine GetKmpAffinityStr(Str,NCPUs,numaNodeCount,processorCoreCount, logicalProcessorCount,AffinityMode,IsClose,StartThread)   
    !This routine builds a custom AFFINITY string, that binds each thread.
    !The actual binding depends on whether AffinityMode is 0 or 1. 
    !  AffinityMode=0 
    !                The binding will be loose and each thread will only be bound to a processor package.
    !                Furthermore all threads will be distributed among the NUMA nodes round-robin style (ie. like dealing cards). (IsClose is not respected in this case) 
    !  AffinityMode=1
    !                Each thread will be bound to a logical processor. 
    !IO
    !Str           -> On output the string will contain a KMP_AFFINITY enviroment variable string.
    !NCPUs         -> Number of threads we wish to use
    !numaNodeCount -> Number of numa nodes in the system
    !LogicalProcessorCount -> Number of logical processors in the system (ie. including hyperthreads).
    !processorCoreCount    -> Number of physical processors in the system (ie. excluding hyperthreads).
    !AffinityMode          -> Affinitymode to generate string for. 0->shared server, 1->dedicated server
    !IsClose               -> If 1, then threads will be distributed on the same NUMA node first until full, before moving on to the next, if not then threads will be spread across all NUMA nodes
    !StartThread           -> which threadID to start the affinity binding with
    implicit none
    character*(*),intent(inout) :: Str      
    integer,intent(in)          :: numaNodeCount,logicalProcessorCount,AffinityMode,processorCoreCount,NCPUs,StartThread,IsClose
    !Local vars
    integer                     :: FirstIdx,LastIdx,i,j,ij,ProcPrNode,thread,bound,inc,CoresPrNuma,NumaID,NInc,reminder,k,numa
    character(len=4)            :: TmpStr
    integer,allocatable         :: usedthreads(:)
    logical                     :: ThreadAccepted, Matchfound
    character*256               :: String
      if (affinityMode.eq.0) then !Shared server
        str="KMP_AFFINITY=verbose,granularity=fine,proclist=["   !48 characters in string
        FirstIdx=49
        LastIdx=49
        ProcPrNode=LogicalProcessorCount/numaNodeCount
        numa=StartThread/ProcPrNode+1
        do i=1,numaNodeCount
          if (i.ne.1) then
            Str(FirstIdx:FirstIdx)=','          
            FirstIdx=FirstIdx+1
          end if      
          Str(FirstIdx:FirstIdx)='{'
           FirstIdx=FirstIdx+1
          do j=1,ProcPrNode
            ij=(numa-1)*ProcPrNode+j-1
            if (j.ne.1) then
              Str(FirstIdx:FirstIdx)=','          
              FirstIdx=FirstIdx+1
            end if    
            LastIdx=FirstIdx
            if (ij.gt.9) then
              LastIdx=FirstIdx+1
            end if
            if (ij.gt.99) then
              LastIdx=FirstIdx+2
            end if
            write (TmpStr, '(I3)') ij
            str(FirstIdx:LastIdx)=TmpStr(3-(LastIdx-FirstIdx):3)
            FirstIdx=LastIdx+1
          end do
          Str(FirstIdx:FirstIdx)='}'
          FirstIdx=FirstIdx+1
          numa=numa+1
          if (numa.gt.numaNodeCount) then
            numa=1
          end if
        end do    
        Str(FirstIdx:FirstIdx+10)="],explicit" 
        print *,Str
      else if (affinityMode.eq.1) then !Dedicated server
        inc=LogicalProcessorCount/processorCoreCount    !Takes care of hyperthreading
        CoresPrNuma=LogicalProcessorCount/NumaNodeCount 
        str="OMP_PLACES="   !12 characters in string
        FirstIdx=12
        LastIdx=12
        bound=min(NCPUs,LogicalProcessorCount)
        thread=StartThread
        allocate(usedthreads(bound))
        usedthreads=-1
        do i=1,bound
          if(i.eq.1) then
            str(FirstIdx:LastIdx)="{"  
          else
            LastIdx=FirstIdx+1
            str(FirstIdx:LastIdx)=",{"
          end if
          FirstIdx=LastIdx+1
          LastIdx=FirstIdx
          do k=1,bound          !We keep looping until we find a thread that is acceptable to all our criterias, or until we run out of threads 
            threadaccepted=.true. 
            matchfound=.false.    !First we check whether we already have bound this particular thread
            do j=1,bound
              if(usedthreads(j).eq.thread) then
                matchfound=.true.
                threadaccepted=.false.  !we take another round
                exit
              end if
            end do
            if (matchfound) then  !If we have bound this before we jump to the next thread
              thread=thread+inc    
            end if
            if(thread.ge.LogicalProcessorCount) then   !If we are trying to bound to a thread higher than the number of threads available, we loop around
              thread=mod(thread,LogicalProcessorCount)
              threadaccepted=.false. !we take another round
            end if
            if(threadaccepted) then
              exit !Finally we found a thread that fits our criteria, so we can exit
            end if
          end do
          if (.not.threadaccepted) then
            !exit   !If you really insist on binding more threads than cpu's just enable this Exit.
            goto 11 !Safety check, we could not find a thread to bind to. This can for instance happen if you try to bind more threads than you have on your system, (hyperthreads do not count)
          else
            usedthreads(i)=thread 
          end if
          !Finally we bind the thread
          if(thread.gt.9) then
            LastIdx=FirstIdx+1
          end if
          if (thread.gt.99) then
            LastIdx=FirstIdx+2
          end if
          write (TmpStr, '(I3)') thread
          str(FirstIdx:LastIdx)=TmpStr(3-(LastIdx-FirstIdx):3)
          FirstIdx=LastIdx+1
          LastIdx=FirstIdx
          str(FirstIdx:LastIdx)="}"
          FirstIdx=LastIdx+1
          !Next up we increase the thread number, since this is distributed, we increase it by jumping up one numa node.
          if(IsClose) then
            thread=thread+inc         !we increase to next non-hyperthread
          else
            thread=thread+CoresPrNuma !we increase to next numa node                 
          end if
        end do  
      end if
      return
  11  print*,'No suitable thread to bind the process to was found, make sure your hardware settings are correct. If they are report this issue.'
  end subroutine GetKmpAffinityStr
    
  function StartOpenMP(iType,iSaveOldThreadNumber) result(UseOpenMP)
    use mMisc
    !This function will return a logical parameter which tells whether further parallelization should occur.
    !Furthermore it will set the number of threads depending on: iType, the level of parallelization set in InitOpenMP and the level of parallelization already enabled.
    !
    !IO
    !iType -> !Determine what kind of parallel region we wish to make
              !iType = 1 : (Default) Adaptive parallelization. 
                  !If no previous parallelization is detected it set the number of threads to NCPUs and locks against further parallelization. 
                  !If 1 layer of previous parallelization is detected it sets the number of threads to NCPUInner.
                  !If 2 or more layers of parallelization is detected it sets the number of threads to 1 and returns UseOpenMP=false.
              !iType = 2 : Adaptive parallelization for memorybandwidth intensive parallelization. Locks against further parallelization. 
                  !If no previous parallelization is detected it set the number of threads to NCPUsLow and locks against further parallelization.  
                  !If 1 layer of previous parallelization is detected it sets the number of threads to NCPUInnerLow.
                  !If 2 or more layers of parallelization is detected it sets the number of threads to 1 and returns UseOpenMP=false.
              !iType = 3 : NUMA spreading. Only Locks against further parallelization if previous parallelization is detected.
                  !If no previous parallelization is detected it set the number of threads to NCPUOuter. 
                  !If 1 layer of previous parallelization is detected it set the number of threads to NCPUInner. 
                  !If 2 or more layers of parallelization is detected it sets the number of threads to 1 and returns UseOpenMP=false.
              !iType = 4 : NUMA spreading for memorybandwidth intensive parallelization. Only Locks against further parallelization if previous parallelization is detected.
                  !If no previous parallelization is detected it spawns NCPUOuter threads  in the parallelization spread across the different NUMA-nodes. 
                  !If 1 layer of previous parallelization is detected it set the number of threads to NCPUInnerLow. 
                  !If 2 or more layers of parallelization is detected it sets the number of threads to 1 and returns UseOpenMP=false.
              !iType = 5 : XeonPhi Not implemented
              !Note only Use iType=3/4 if you are sure further parallelization will occur before the CPU intensive tasks 
    !iSaveOldThreadNumber -> if set to 1 it saves the old threadnumber such that it can be restored after the parallel region is done, if desired. In order to recall it call, the subroutine ReturnOldThreadNumber()
    integer,intent(in)             :: iType 
    integer,intent(in),optional    :: iSaveOldThreadNumber 
    !Local vars
    logical  :: UseOpenMP 
    integer  :: i,saveOldThreadNumber,CPUs   
      !Do we wish to save the old threadnumber before changing it?
      SaveOldThreadNumber=LoadOptionalParam(0,iSaveOldThreadNumber)
      if(SaveOldThreadNumber.eq.1) then
          OldThreadNumber = omp_get_max_threads()
      end if
      i=omp_get_active_level() !Get previous level of parallelization
      if (i.eq.0) then !No previous parallelization detected
        RunningSingleParallelRegion = .false. !We unlock any potential previous parallelization locks
        UseOpenMP=.true.
        SELECT CASE (iType) 
        CASE(1)
          RunningSingleParallelRegion = .true.
          CPUs=NCPUs
          call omp_set_num_threads(NCPUs)
        CASE(2)
          RunningSingleParallelRegion = .true.
          call omp_set_num_threads(NCPUsLow)
        CASE(3,4)
          call omp_set_num_threads(NCPUOuter)
          RunningSingleParallelRegion=.false.
        END SELECT
      elseif(i.eq.1.AND.UseNested.AND.(.not.RunningSingleParallelRegion)) then !Outer parallelization detected, but room for inner parallelization
        UseOpenMP=.true.
        SELECT CASE (iType) 
        CASE(1,3)
          call omp_set_num_threads(NCPUInner)
        CASE(2,4)
          call omp_set_num_threads(NCPUInnerLow)     
        END SELECT
      else !Full parallelization already detected
        UseOpenMP=.false.
        call omp_set_num_threads(1)
      end if
  end function StartOpenMP

  subroutine ReturnOldThreadNum()
    !This returns the thread number to the value saved in the module variable OldThreadNumber
      call omp_set_num_threads(OldThreadNumber)
  end subroutine ReturnOldThreadNum

end module mParallel
    
    
