module m_IOAndVisu
#ifdef OPENMP
  use omp_lib
#endif
  use maincode_module
  use m_Mesh
  implicit none
  private

  integer,parameter,public:: mytype_save = KIND(0.0D0)
  public:: InitVisu, dump_visu

contains

  !******************************************************************
  ! InitVisu
  !******************************************************************
  subroutine InitVisu
    implicit none

    ! locals
    character(128)::XdmfFile
    integer::nUnitFile,ierror,nflds,ifld,iprec,i,j,k
 
    ! write XDMF file
    write(xdmfFile,"(A)") "VisuFor3DPDR.xmf"
    open(newunit=nUnitFile, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror/=0) print*,"Cannot open xmf file"
    ! XDMF/XMF Title
    write(nUnitFile,'(A)') '<?xml version="1.0" ?>'
    write(nUnitFile,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(nUnitFile,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(nUnitFile,'(A)') '<Domain>'

    ! grid
    iprec=mytype_save
    write(nUnitFile,'(A,3I7,A)')'    <Topology name="TOPO" TopologyType="3DRectMesh" Dimensions="',nzc,nyc,nxc,'"/>'
    write(nUnitFile,'(A)')'    <Geometry name="GEO" GeometryType="VXVYVZ">'
    ! x-grid
    write(nUnitFile,'(A,I1,A,I5,A)') '        <DataItem Format="XML" DataType="Float" Precision="',&
    &iprec,'" Endian="Native" Dimensions="',nxc,'">'
    write(nUnitFile,'(A)',advance='no') '        '
    do i=1,nxc
      write(nUnitFile,'(E14.7)',advance='no') (i-1)*dx+dx*0.5_RK
    enddo
    write(nUnitFile,'(A)')' '; write(nUnitFile,'(A)')'        </DataItem>'
    ! y-grid
    write(nUnitFile,'(A,I1,A,I5,A)') '        <DataItem Format="XML" DataType="Float" Precision="',&
    &iprec,'" Endian="Native" Dimensions="',nyc,'">'
    write(nUnitFile,'(A)',advance='no') '        '
    do j=1,nyc
      write(nUnitFile,'(E14.7)',advance='no') (j-1)*dy+dy*0.5_RK
    enddo
    write(nUnitFile,'(A)')' '; write(nUnitFile,'(A)')'        </DataItem>'
    ! z-grid
    write(nUnitFile,'(A,I1,A,I5,A)') '        <DataItem Format="XML" DataType="Float" Precision="',&
    &iprec,'" Endian="Native" Dimensions="',nzc,'">'
    write(nUnitFile,'(A)',advance='no') '        '
    do k=1,nzc
      write(nUnitFile,'(E14.7)',advance='no') (k-1)*dz+dz*0.5_RK
    enddo
    write(nUnitFile,'(A)')' '; write(nUnitFile,'(A)')'        </DataItem>'
    write(nUnitFile,'(A)')'    </Geometry>'

    ! Time series
    nflds = 1
    write(nUnitFile,'(A)')'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    write(nUnitFile,'(A)')'        <Time TimeType="List">'
    write(nUnitFile,'(A,I6,A)')'        <DataItem Format="XML" NumberType="Int" Dimensions="',nflds,'">' 
    write(nUnitFile,'(A)',advance='no')'        '
    ifld = 0
    write(nUnitFile,'(I10)',advance='no') ifld
    write(nUnitFile,'(A)')'        </DataItem>'
    write(nUnitFile,'(A)')' '; write(nUnitFile,'(A)')'       </Time>'

    ! attribute
    write(nUnitFile,'(A,I10.10,A)')'        <Grid Name="T',ifld,'" GridType="Uniform">'
    write(nUnitFile,'(A)')'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
    write(nUnitFile,'(A)')'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
    ! add virables
    call Write_XDMF_One(nUnitFile,ifld,'Density')
    call Write_XDMF_One(nUnitFile,ifld,'Tgas')
    write(nUnitFile,'(A)')'        </Grid>'

    write(nUnitFile,'(A)')'    </Grid>'
    write(nUnitFile,'(A)')'</Domain>'
    write(nUnitFile,'(A)')'</Xdmf>'
    close(nUnitFile,IOSTAT=ierror)
  end subroutine InitVisu

  !******************************************************************
  ! Write_XDMF_One
  !******************************************************************
  subroutine Write_XDMF_One(nUnitFile, ifld,chAttribute)
    implicit none
    integer,intent(in)::nUnitFile,ifld
    character(*),intent(in)::chAttribute

    ! locals
    character(128)::chFile
    integer::iprec=mytype_save

    write(chFile,'(A,A,I10.10)')"VisuFor3DPDR","_"//trim(adjustl(chAttribute))//"_",ifld
    write(nUnitFile,'(A)')'            <Attribute Name="'//trim(chAttribute)//'" Center="Node">'
    write(nUnitFile,'(A,I1,A,3I7,A)')'                <DataItem Format="Binary" DataType="Float" Precision="',&
    &iprec,'" Endian="Native" Dimensions="',nzc,nyc,nxc,'">'
    write(nUnitFile,'(A)')'                    '//trim(chFile)
    write(nUnitFile,'(A)')'                </DataItem>'
    write(nUnitFile,'(A)')'            </Attribute>'
  end subroutine Write_XDMF_One

  !******************************************************************
  ! dump_visu
  !******************************************************************
  subroutine dump_visu
    implicit none

    ! locals
    character(128)::chFile
    integer::I,J,K,p,cI,ntime,nUnit,ierror
    real(RK),allocatable,dimension(:) :: GridCell

    ntime = 0

    ! write Density
    allocate(GridCell(nxc*nyc*nzc))
    GridCell=0.D0
#ifdef XYZ
!$omp parallel do
    do p=1,pdr_ptot
      ! GridCell(p)=log10(pdr(p)%rho)
      GridCell(p)=pdr(p)%rho
    enddo
!$omp end parallel do
#else
!$omp parallel do private (p,cI)
    do K = 1,nzc
      do J = 1,nyc
        do I = 1,nxc
          p = I + (J-1)*nxc + (K-1)*(nxc*nyc)
          cI = K + (J-1)*nzc + (I-1)*(nzc*nyc)
          GridCell(p) = pdr(cI)%rho
        enddo
      enddo
    enddo
!$omp end parallel do
#endif

    write(chFile,"(A,I10.10)") "VisuFor3DPDR_Density_",ntime
    open(newunit=nUnit, file=chFile, status='replace', form='unformatted', access='stream')
      write(nUnit) GridCell
    close(nUnit)

    ! write Tgas
    GridCell=0.D0
#ifdef XYZ
!$omp parallel do
    do p=1,pdr_ptot
      GridCell(p)=pdr(p)%Tgas
    enddo
!$omp end parallel do
#else
!$omp parallel do private (p,cI)
    do K = 1,nzc
      do J = 1,nyc
        do I = 1,nxc
          p = I + (J-1)*nxc + (K-1)*(nxc*nyc)
          cI = K + (J-1)*nzc + (I-1)*(nzc*nyc)
          GridCell(p) = pdr(cI)%Tgas
        enddo
      enddo
    enddo
!$omp end parallel do
#endif

    write(chFile,"(A,I10.10)") "VisuFor3DPDR_Tgas_",ntime
    open(newunit=nUnit, file=chFile, status='replace', form='unformatted', access='stream')
      write(nUnit) GridCell
    close(nUnit)
    deallocate(GridCell)

!     ! write deltau
!     allocate(GridCell(nxc,nyc,nzc),GridCell_temp(nxc,nyc,nzc))
!     GridCell=0.D0; GridCell_temp=0.D0
! !$omp parallel do private(II,JJ,KK)
!     do k=1,nznp
!       do j=1,nynp
!         do i=1,nxnp
!         II=i+IID*nxnp
!         JJ=j+JID*nynp
!         KK=k+KID*nznp
!         GridCell_temp(II,JJ,KK)=Deltau(i,j,k)
!         enddo
!       enddo
!     enddo
! !$omp end parallel do
!     call MPI_REDUCE(GridCell_temp,GridCell,nxc*nyc*nzc,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
!     ! call MPI_ALLREDUCE(GridCell_temp,GridCell,nxc*nyc*nzc,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)

!     if(nrank==0) then
!       write(chFile,"(A,I10.10)") "VisuForVolume_deltau_",ntime
!       open(newunit=nUnit, file=chFile, status='replace', form='unformatted', access='stream')
!         write(nUnit) GridCell
!       close(nUnit)
!     endif
!     deallocate(GridCell,GridCell_temp)

!     ! write J0
!     allocate(GridCell(nxc,nyc,nzc),GridCell_temp(nxc,nyc,nzc))
!     GridCell=0.D0; GridCell_temp=0.D0
! !$omp parallel do private(II,JJ,KK)
!     do k=1,nznp
!       do j=1,nynp
!         do i=1,nxnp
!         II=i+IID*nxnp
!         JJ=j+JID*nynp
!         KK=k+KID*nznp
! #ifndef SPHERE
!         GridCell_temp(II,JJ,KK)=log10(spop(1,i,j,k))
! #else
!         GridCell_temp(II,JJ,KK)=spop(1,i,j,k)
! #endif
!         enddo
!       enddo
!     enddo
! !$omp end parallel do
!     call MPI_REDUCE(GridCell_temp,GridCell,nxc*nyc*nzc,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
!     ! call MPI_ALLREDUCE(GridCell_temp,GridCell,nxc*nyc*nzc,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)

!     if(nrank==0) then
!       open(newunit=nUnit,file='spop1.dat')
!       do k = 1,nzc
!         do j = 1,nyc
!           do i = 1,nxc
!             xc=(real(i,kind=RK)-0.5D0)*ndSize
!             yc=(real(j,kind=RK)-0.5D0)*ndSize
!             zc=(real(k,kind=RK)-0.5D0)*ndSize
!             rr=radius-sqrt((xc-xp)**2.D0+(yc-yp)**2.D0+(zc-zp)**2.D0)
!             if(k==nzc/2.and.j==nyc/2.and.i.ge.nxc/2) then           
!               write(nUnit,'(1ES11.3)') GridCell(i,j,k)
!             endif
!           enddo
!         enddo
!       enddo
!       close(nUnit)

!       write(chFile,"(A,I10.10)") "VisuForVolume_J0_",ntime
!       open(newunit=nUnit, file=chFile, status='replace', form='unformatted', access='stream')
!         write(nUnit) GridCell
!       close(nUnit)
!     endif
!     deallocate(GridCell,GridCell_temp)

!     ! write J1
!     allocate(GridCell(nxc,nyc,nzc),GridCell_temp(nxc,nyc,nzc))
!     GridCell=0.D0; GridCell_temp=0.D0
! !$omp parallel do private(II,JJ,KK)
!     do k=1,nznp
!       do j=1,nynp
!         do i=1,nxnp
!         II=i+IID*nxnp
!         JJ=j+JID*nynp
!         KK=k+KID*nznp
! #ifndef SPHERE
!         GridCell_temp(II,JJ,KK)=log10(spop(2,i,j,k))
! #else
!         GridCell_temp(II,JJ,KK)=spop(2,i,j,k)
! #endif
!         enddo
!       enddo
!     enddo
! !$omp end parallel do
!     call MPI_REDUCE(GridCell_temp,GridCell,nxc*nyc*nzc,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
!     ! call MPI_ALLREDUCE(GridCell_temp,GridCell,nxc*nyc*nzc,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)

!     if(nrank==0) then
!       open(newunit=nUnit,file='spop2.dat')
!       do k = 1,nzc
!         do j = 1,nyc
!           do i = 1,nxc
!             xc=(real(i,kind=RK)-0.5D0)*ndSize
!             yc=(real(j,kind=RK)-0.5D0)*ndSize
!             zc=(real(k,kind=RK)-0.5D0)*ndSize
!             rr=radius-sqrt((xc-xp)**2.D0+(yc-yp)**2.D0+(zc-zp)**2.D0)
!             if(k==nzc/2.and.j==nyc/2.and.i.ge.nxc/2) then           
!               write(nUnit,'(1ES11.3)') GridCell(i,j,k)
!             endif
!           enddo
!         enddo
!       enddo
!       close(nUnit)

!       write(chFile,"(A,I10.10)") "VisuForVolume_J1_",ntime
!       open(newunit=nUnit, file=chFile, status='replace', form='unformatted', access='stream')
!         write(nUnit) GridCell
!       close(nUnit)
!     endif
!     deallocate(GridCell,GridCell_temp)

!     ! write optical depth and escaprprobability
!     allocate(GridCellPix(nxc,nyc,nzc,0:tr_nPix-1),GridCellPix_temp(nxc,nyc,nzc,0:tr_nPix-1))
!     allocate(GridCell(nxc,nyc,nzc),GridCell_temp(nxc,nyc,nzc))
!     GridCellPix=0.D0; GridCell_temp=0.D0
!     GridCell=0.D0; GridCell_temp=0.D0
! !$omp parallel do private(II,JJ,KK)
!     do k=1,nznp
!       do j=1,nynp
!         do i=1,nxnp
!         II=i+IID*nxnp
!         JJ=j+JID*nynp
!         KK=k+KID*nznp
!         GridCell_temp(II,JJ,KK)=beta(i,j,k)
!           do il=0,tr_nPix-1
!             GridCellPix_temp(II,JJ,KK,il)=OpticalDepth(i,j,k,il)
!           enddo
!         enddo
!       enddo
!     enddo
! !$omp end parallel do
!     call MPI_REDUCE(GridCell_temp,GridCell,nxc*nyc*nzc,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
!     call MPI_REDUCE(GridCellPix_temp,GridCellPix,nxc*nyc*nzc*tr_nPix,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
!     ! call MPI_ALLREDUCE(GridCell_temp,GridCell,nxc*nyc*nzc,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
!     ! call MPI_ALLREDUCE(GridCellPix_temp,GridCellPix,nxc*nyc*nzc*tr_nPix,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)

!     ! find ray index
!     max_dotprod = 0.0
!     xnode = 1.D0
!     ynode = 0.D0
!     znode = 0.D0
!     do ipix = 0,tr_nPix-1
!       call pix2ang_nest(tr_nSide, ipix, thfpix, phfpix)
!       x = sin(thfpix)*cos(phfpix)
!       y = sin(thfpix)*sin(phfpix)
!       z = cos(thfpix)
!       dotprod = x*xnode + y*ynode + z*znode
!       if (dotprod > max_dotprod) then
!         max_dotprod = dotprod
!         iclosest = ipix
!       endif
!     enddo
!     call pix2ang_nest(tr_nSide, iclosest, thfpix, phfpix)
!     if(nrank==0) print*,'direction',iclosest,thfpix,phfpix

!     if(nrank==0) then
!       open(newunit=nUnit,file='escape.dat')
!       do k = 1,nzc
!         do j = 1,nyc
!           do i = 1,nxc
!             xc=(real(i,kind=RK)-0.5D0)*ndSize
!             yc=(real(j,kind=RK)-0.5D0)*ndSize
!             zc=(real(k,kind=RK)-0.5D0)*ndSize
!             rr=radius-sqrt((xc-xp)**2.D0+(yc-yp)**2.D0+(zc-zp)**2.D0)
!             if(k==nzc/2.and.j==nyc/2.and.i.ge.nxc/2) then           
!               write(nUnit,*) rr,GridCellPix(i,j,k,iclosest),GridCell(i,j,k)
!             endif
!           enddo
!         enddo
!       enddo
!       close(nUnit)

!       write(chFile,"(A,I10.10)") "VisuForVolume_beta_",ntime
!       open(newunit=nUnit, file=chFile, status='replace', form='unformatted', access='stream')
!         write(nUnit) GridCell
!       close(nUnit)
!     endif
!     deallocate(GridCellPix,GridCellPix_temp)
!     deallocate(GridCell,GridCell_temp)

  end subroutine dump_visu

end module m_IOAndVisu
