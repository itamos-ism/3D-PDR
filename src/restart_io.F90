! Strict iteration-boundary checkpoints: B commits a prefix of append-only A.
! Validate both files completely before applying any cell state.
!
! A stores each frozen cell once; B.tmp stores active cells and global solver
! state, then RENAME commits it as B. B records the byte end/count of its A
! prefix, so an interrupted append never invalidates the previous checkpoint.
! Ignore any A tail on load and overwrite it only on the next save. Do not
! promote B.tmp: without a committed B the user must explicitly start afresh.
!
! Native format version 1 requires the same executable, parameter/input and
! chemistry/coolant files, level dimensions and OpenMP thread count. Old,
! incomplete, corrupt or incompatible checkpoints stop with an error, never
! silently cold-start. CRC64 detects accidental corruption, not malicious edits.
! CLOSE/RENAME protect against process interruption, not sudden power loss.
module restart_checkpoint
  use iso_fortran_env, only : restart_int64 => int64, int8
  use ieee_arithmetic, only : ieee_is_finite
  use maincode_module
  implicit none
  private
  public :: checkpoint_setup, checkpoint_load, checkpoint_save, checkpoint_finished
  integer, parameter :: format_version=1
  character(len=8), parameter :: magic='PDRST001'
  integer(restart_int64) :: committed_end=1, pair_id=0, signature=0
  integer :: archived=0
  logical :: signature_ready=.false., checkpoint_finished=.false.
  integer(restart_int64) :: crc_table(0:255)
  logical :: crc_ready=.false.
contains
  subroutine prepare_crc
    integer :: j,bit
    integer(restart_int64) :: crc
    if (crc_ready) return
    do j=0,255
      crc=int(j,restart_int64)
      do bit=1,8
        if (btest(crc,0)) then
          crc=ieor(shiftr(crc,1),int(z'C96C5795D7870F42',restart_int64))
        else
          crc=shiftr(crc,1)
        endif
      enddo
      crc_table(j)=crc
    enddo
    crc_ready=.true.
  end subroutine

  function hash_buffer(data,length,seed) result(hash)
    real(dp), intent(in) :: data(:)
    integer(restart_int64), intent(in) :: length,seed
    integer(restart_int64) :: hash,bits
    integer :: j,byte,index
    call prepare_crc
    call require(length==8_restart_int64*size(data),'invalid checksum buffer length')
    hash=seed
    do j=1,size(data)
      bits=transfer(data(j),bits)
      do byte=1,8
        index=int(iand(ieor(hash,bits),255_restart_int64))
        hash=ieor(shiftr(hash,8),crc_table(index))
        bits=shiftr(bits,8)
      enddo
    enddo
  end function

  integer function hash_file(path,length,hash) result(rc)
    character(*), intent(in) :: path
    integer(restart_int64), intent(in) :: length
    integer(restart_int64), intent(out) :: hash
    integer(int8) :: bytes(65536)
    integer(restart_int64) :: remaining,total
    integer :: unit,n,j,index,ios
    call prepare_crc
    rc=1
    hash=0
    open(newunit=unit,file=path,access='stream',form='unformatted',status='old',action='read',iostat=ios)
    if (ios/=0) return
    inquire(unit=unit,size=total)
    remaining=total
    if (length>=0) remaining=length
    if (remaining>total) then
      close(unit)
      return
    endif
    do while (remaining>0)
      n=int(min(remaining,int(size(bytes),restart_int64)))
      read(unit,iostat=ios) bytes(:n)
      if (ios/=0) then
        close(unit)
        return
      endif
      do j=1,n
        index=int(iand(ieor(hash,int(bytes(j),restart_int64)),255_restart_int64))
        hash=ieor(shiftr(hash,8),crc_table(index))
      enddo
      remaining=remaining-n
    enddo
    close(unit,iostat=rc)
  end function

  integer function atomic_replace(temporary,target) result(rc)
    character(*), intent(in) :: temporary,target
    ! Both paths are in the same directory/filesystem. RENAME is the commit
    ! point; CLOSE flushes process buffers but does not guarantee power-loss durability.
    call rename(temporary,target,rc)
  end function

  subroutine require(ok,message)
    logical, intent(in) :: ok
    character(*), intent(in) :: message
    if (ok) return
    write(*,'(A)') ' [RESTART] ERROR: '//message
    error stop 1
  end subroutine

  subroutine checkpoint_setup
#if defined RESTART && defined THERMALBALANCE
    logical :: have_a,have_tmp
    inquire(file=restart_file_B,exist=restart)
    inquire(file=restart_file_A,exist=have_a)
    inquire(file=trim(restart_file_B)//'.tmp',exist=have_tmp)
    call require(restart.or..not.(have_a.or.have_tmp), &
         'no committed B checkpoint; preserve/remove orphan A/tmp before a new run')
#endif
  end subroutine

  subroutine model_signature
#if defined RESTART && defined THERMALBALANCE
    character(len=4096) :: executable
    character(len=16) :: network
    integer :: j,stat
    if (signature_ready) return
    ! Native binary format: require the exact executable and input data.
    call get_command_argument(0,executable,status=stat)
    call require(stat==0,'cannot identify executable')
    signature=0
    call add_file(trim(executable))
    call add_file(trim(paramFile))
    call add_file(trim(input))
#ifdef REDUCED
    network='reduced'
#elif MEDIUM
    network='medium'
#elif FULL
    network='full'
#else
    network='mynetwork'
#endif
    call add_file('chemfiles/species_'//trim(network)//'.d')
    call add_file('chemfiles/rates_'//trim(network)//'.d')
    do j=1,coo
      call add_file(trim(coolfile(j)))
    enddo
    signature_ready=.true.
#endif
  end subroutine

  subroutine add_file(path)
    character(*), intent(in) :: path
    integer(restart_int64) :: h
    call require(hash_file(path,-1_restart_int64,h)==0,'cannot hash '//path)
    signature=ieor(ishftc(signature,7),h)
  end subroutine

  integer function cell_size() result(n)
    integer :: j,nl
    n=5+nspec+coo+12+5*nrays
#ifdef THERMALBALANCE
    n=n+7
#ifdef ILLINOIS
    n=n+3
#endif
#ifdef REBRACKET
    n=n+1
#endif
#endif
#ifdef CHEMANALYSIS
    n=n+nreac
#endif
    do j=1,coo
      nl=coolant(j)%cnlev
      n=n+2+3*nl+nl*nl
#ifdef NGACCEL
      n=n+3*nl
#endif
#ifdef NGRELAX
      n=n+1
#endif
    enddo
  end function

  subroutine cell_state(point,data,restore)
    integer, intent(in) :: point
    real(dp), intent(inout) :: data(:)
    logical, intent(in) :: restore
    integer :: pos,j,col
    pos=0
    if (.not.restore) data=0
    call scalar(pdr(point)%Tgas)
    call scalar(pdr(point)%nTgas)
    call scalar(pdr(point)%totalcooling)
    call flag(pdr(point)%levelconverged)
    call letter(pdr(point)%previouschange)
#ifdef THERMALBALANCE
    call flag(pdr(point)%fullyconverged)
    call flag(pdr(point)%doleveltmin)
    call flag(pdr(point)%dobinarychop)
    call scalar(pdr(point)%Tlow)
    call scalar(pdr(point)%Thigh)
    call scalar(pdr(point)%Fmean)
    call scalar(pdr(point)%Fratio)
#ifdef ILLINOIS
    call scalar(pdr(point)%Flow)
    call scalar(pdr(point)%Fhigh)
    call letter(pdr(point)%ill_last)
#endif
#ifdef REBRACKET
    call number(pdr(point)%nrebracket)
#endif
#endif
    call vector(pdr(point)%abundance)
    call vector(pdr(point)%cooling)
    call vector(pdr(point)%heating)
    call vector(pdr(point)%column_NH2)
    call vector(pdr(point)%column_NHD)
    call vector(pdr(point)%column_NCO)
    call vector(pdr(point)%column_NC)
    call vector(pdr(point)%column_NS)
#ifdef CHEMANALYSIS
    call vector(temp_rate(:,point))
#endif
    do j=1,coo
      call flag(pdr(point)%coolant(j)%isconverged)
      call scalar(pdr(point)%coolant(j)%coolprev)
      call vector(pdr(point)%coolant(j)%pop)
      call vector(pdr(point)%coolant(j)%solution)
      call vector(pdr(point)%coolant(j)%relativechange)
      do col=1,coolant(j)%cnlev
        call vector(pdr(point)%coolant(j)%line(:,col))
      enddo
#ifdef NGACCEL
      do col=1,3
        call vector(pdr(point)%coolant(j)%pophist(:,col))
      enddo
#endif
#ifdef NGRELAX
      call number(pdr(point)%coolant(j)%noscil)
#endif
    enddo
    call require(pos==size(data),'internal checkpoint cell layout mismatch')
  contains
    subroutine scalar(value)
      real(dp), intent(inout) :: value
      pos=pos+1
      if (restore) then
        value=data(pos)
      else
        data(pos)=value
      endif
    end subroutine
    subroutine vector(value)
      real(dp), intent(inout) :: value(:)
      if (restore) then
        value=data(pos+1:pos+size(value))
      else
        data(pos+1:pos+size(value))=value
      endif
      pos=pos+size(value)
    end subroutine
    subroutine number(value)
      integer, intent(inout) :: value
      real(dp) :: v
      v=real(value,dp)
      call scalar(v)
      if (restore) value=nint(v)
    end subroutine
    subroutine flag(value)
      logical, intent(inout) :: value
      real(dp) :: v
      v=0
      if (value) v=1
      call scalar(v)
      if (restore) value=(v==1)
    end subroutine
    subroutine letter(value)
      character, intent(inout) :: value
      integer :: v
      v=iachar(value)
      call number(v)
      if (restore) value=achar(v)
    end subroutine
  end subroutine

  subroutine global_state(data,restore,finished)
    real(dp), intent(inout) :: data(:)
    logical, intent(in) :: restore,finished
    integer :: j
    if (restore) then
      iteration=nint(data(1))
      levpop_iteration=nint(data(2))
#ifdef THERMALBALANCE
      first_time=data(3)==1
      level_conv=data(4)==1
#endif
#ifdef NGACCEL
      ng_nhist=nint(data(5))
#endif
      checkpoint_finished=data(6)==1
      thermal_percentage=data(7)
      levpop_percentage=data(8)
      do j=1,coo
        coolant(j)%percentage=real(data(16+j))
      enddo
    else
      data=0
      data(1)=iteration
      data(2)=levpop_iteration
#ifdef THERMALBALANCE
      data(3)=merge(1,0,first_time)
      data(4)=merge(1,0,level_conv)
#endif
#ifdef NGACCEL
      data(5)=ng_nhist
#endif
      data(6)=merge(1,0,finished)
      data(7)=thermal_percentage
      data(8)=levpop_percentage
      do j=1,coo
        data(16+j)=coolant(j)%percentage
      enddo
    endif
  end subroutine

  subroutine put_cell(unit,point,data)
    integer, intent(in) :: unit,point
    real(dp), intent(inout) :: data(:)
    integer :: ios
    integer(restart_int64) :: h
    call cell_state(point,data,.false.)
    call require(all(ieee_is_finite(data)),'non-finite cell state cannot be checkpointed')
    h=hash_buffer(data,8_restart_int64*size(data),int(point,restart_int64))
    write(unit,iostat=ios) point,data,h
    call require(ios==0,'writing cell record')
  end subroutine

  subroutine checkpoint_save(finished)
    logical, intent(in) :: finished
#if defined RESTART && defined THERMALBALANCE
    real(dp), allocatable :: data(:),globals(:)
    integer :: a,b,ios,j,new_count,threads
    integer(restart_int64) :: next_end,h
    character(len=:), allocatable :: temporary
    call model_signature
    allocate(data(cell_size()),globals(16+coo))
    if (committed_end==1) then
      call system_clock(pair_id)
      pair_id=ieor(pair_id,signature)
      open(newunit=a,file=restart_file_A,access='stream',form='unformatted',status='replace',iostat=ios)
      call require(ios==0,'creating A checkpoint')
      write(a,iostat=ios) magic,format_version,pair_id,signature
      call require(ios==0,'writing A header')
    else
      open(newunit=a,file=restart_file_A,access='stream',form='unformatted',status='old',iostat=ios)
      call require(ios==0,'opening A checkpoint')
      ! Position at the committed prefix; subsequent ENDFILE discards only
      ! uncommitted tail records, never data referenced by the previous B.
      write(a,pos=committed_end,iostat=ios)
      call require(ios==0,'positioning A checkpoint')
    endif
    new_count=archived
    do j=1,pdr_ptot
      if (.not.pdr(j)%fullyconverged.or.pdr(j)%restconverged) cycle
      call put_cell(a,j,data)
      new_count=new_count+1
    enddo
    inquire(unit=a,pos=next_end)
    endfile(a,iostat=ios)
    call require(ios==0,'truncating uncommitted A tail')
    close(a,iostat=ios)
    call require(ios==0,'closing A checkpoint')

    temporary=trim(restart_file_B)//'.tmp'
    open(newunit=b,file=temporary,access='stream',form='unformatted',status='replace',iostat=ios)
    call require(ios==0,'creating temporary B checkpoint')
    threads=1
#ifdef OPENMP
    threads=CPUs
#endif
    write(b,iostat=ios) magic,format_version,pair_id,signature,pdr_ptot,coo,nspec,nreac,nrays, &
         lmax_levels,threads,size(data),coolant(:)%cnlev,new_count,pdr_ptot-new_count,next_end
    call require(ios==0,'writing B header')
    call global_state(globals,.false.,finished)
    write(b,iostat=ios) globals
    call require(ios==0,'writing global state')
    do j=1,pdr_ptot
      if (.not.pdr(j)%fullyconverged) call put_cell(b,j,data)
    enddo
    close(b,iostat=ios)
    call require(ios==0,'closing temporary B checkpoint')
    call require(hash_file(temporary,-1_restart_int64,h)==0,'hashing B checkpoint')
    open(newunit=b,file=temporary,access='stream',form='unformatted',status='old',position='append',iostat=ios)
    call require(ios==0,'opening B for completion marker')
    write(b,iostat=ios) h,magic
    call require(ios==0,'writing B completion marker')
    close(b,iostat=ios)
    call require(ios==0,'closing completed B checkpoint')
    call require(atomic_replace(temporary,trim(restart_file_B))==0, &
         'atomically committing B checkpoint')
    committed_end=next_end
    archived=new_count
    pdr(:)%restconverged=pdr(:)%fullyconverged
    write(*,'(A,I0,A,I0,A,I0)') ' [RESTART] committed iteration ',iteration,': archived=',archived, &
         ', active=',pdr_ptot-archived
#endif
  end subroutine

  subroutine checkpoint_load
#if defined RESTART && defined THERMALBALANCE
    integer :: a,b,ios,v,j,id,pass,na,nb,threads,meta(8),nl(coo)
    integer(restart_int64) :: token,sig,endpoint,bsize,h,stored,apos,bpos,pos,asize
    real(dp), allocatable :: data(:),globals(:)
    logical, allocatable :: seen(:)
    character(len=8) :: tag
    call model_signature
    open(newunit=b,file=restart_file_B,access='stream',form='unformatted',status='old',iostat=ios)
    call require(ios==0,'opening committed B checkpoint')
    inquire(unit=b,size=bsize)
    call require(bsize>=16,'truncated B checkpoint')
    read(b,pos=bsize-15,iostat=ios) stored,tag
    call require(ios==0.and.tag==magic,'missing B completion marker (old/truncated format)')
    call require(hash_file(trim(restart_file_B),bsize-16,h)==0,'hashing committed B')
    call require(h==stored,'B checkpoint checksum mismatch')
    read(b,pos=1,iostat=ios) tag,v,token,sig,meta,nl,na,nb,endpoint
    call require(ios==0,'truncated B header')
    call require(tag==magic.and.v==format_version,'unsupported checkpoint format')
    call require(sig==signature,'executable, parameters or input data differ from checkpoint')
    threads=1
#ifdef OPENMP
    threads=CPUs
#endif
    call require(all(meta==[pdr_ptot,coo,nspec,nreac,nrays,lmax_levels,threads,cell_size()]), &
         'checkpoint dimensions, LEVMAX or thread count differ')
    call require(all(nl==coolant(:)%cnlev),'coolant level counts differ')
    call require(na>=0.and.na<=pdr_ptot.and.nb==pdr_ptot-na,'invalid cell counts')
    allocate(data(cell_size()),globals(16+coo),seen(pdr_ptot))
    read(b,iostat=ios) globals
    call require(ios==0,'truncated global state')
    call require(all(ieee_is_finite(globals)),'non-finite global state')
    call require(globals(1)>=1.and.globals(1)<=itertot,'invalid iteration counter')
    inquire(unit=b,pos=bpos)
    open(newunit=a,file=restart_file_A,access='stream',form='unformatted',status='old',iostat=ios)
    call require(ios==0,'missing A checkpoint')
    read(a,iostat=ios) tag,v,pair_id,signature
    call require(ios==0,'truncated A header')
    call require(tag==magic.and.v==format_version.and.pair_id==token.and.signature==sig, &
         'A/B checkpoint pair mismatch')
    inquire(unit=a,pos=apos,size=asize)
    call require(endpoint==apos+int(na,restart_int64)*(12+8_restart_int64*size(data)), &
         'invalid committed A prefix')
    call require(asize>=endpoint-1,'truncated committed A prefix')
    call require(bsize==bpos-1+int(nb,restart_int64)*(12+8_restart_int64*size(data))+16, &
         'invalid B payload length')
    ! Validate all IDs, lengths, checksums and coverage before applying state.
    do pass=1,2
      seen=.false.
      pos=apos
      do j=1,na
        call get_cell(a,pos,.true.)
      enddo
      pos=bpos
      do j=1,nb
        call get_cell(b,pos,.false.)
      enddo
      call require(all(seen),'checkpoint is missing cells')
    enddo
    close(a,iostat=ios)
    call require(ios==0,'closing A after restore')
    close(b,iostat=ios)
    call require(ios==0,'closing B after restore')
    call global_state(globals,.true.,.false.)
    committed_end=endpoint
    archived=na
    write(*,'(A,I0)') ' [RESTART] validated complete checkpoint at iteration ',iteration
    if (asize>endpoint-1) write(*,'(A)') ' [RESTART] ignoring uncommitted A tail'
  contains
    subroutine get_cell(unit,position,frozen)
      integer, intent(in) :: unit
      integer(restart_int64), intent(inout) :: position
      logical, intent(in) :: frozen
      read(unit,pos=position,iostat=ios) id,data,stored
      call require(ios==0,'truncated cell payload')
      call require(id>=1.and.id<=pdr_ptot,'cell index out of range')
      call require(.not.seen(id),'duplicate cell index')
      seen(id)=.true.
      h=hash_buffer(data,8_restart_int64*size(data),int(id,restart_int64))
      call require(h==stored,'cell checksum mismatch')
      call require(all(ieee_is_finite(data)),'non-finite cell payload')
      call require(data(6)==merge(1.0_dp,0.0_dp,frozen),'cell in wrong checkpoint partition')
      if (pass==2) then
        call cell_state(id,data,.true.)
        pdr(id)%restconverged=frozen
      endif
      position=position+12+8_restart_int64*size(data)
    end subroutine
#endif
  end subroutine

end module restart_checkpoint
