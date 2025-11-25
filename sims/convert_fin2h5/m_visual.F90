module m_visual
  use m_paramters
  implicit none
  private

  integer,parameter :: mytype_save = KIND(0.0D0)
  public:: initvisual

contains

  !******************************************************************
  ! InitVisu
  !******************************************************************
  subroutine initvisual
    implicit none

    ! locals
    character(128)::XdmfFile, out_file_h5
    integer::nUnitFile,ierror
 
    ! write XDMF file
    write(xdmfFile,"(A)") "VisuFor3DPDR.xmf"
    open(newunit=nUnitFile, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror/=0) print*,"Cannot open xmf file"
    ! XDMF/XMF Title
    write(nUnitFile,'(A)') '<?xml version="1.0" ?>'
    write(nUnitFile,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(nUnitFile,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.0">'
    write(nUnitFile,'(A)') '  <Domain>'

    ! grid
    write(nUnitFile,'(A)')'    <Grid name="Mesh" GridType="Uniform">'
    write(nUnitFile,'(A,3I7,A)')'      <Topology TopologyType="3DCoRectMesh" Dimensions="',nxc,nyc,nzc,'"/>'
    write(nUnitFile,'(A)')'      <Geometry Type="ORIGIN_DXDYDZ">'
    write(nUnitFile,'(A,3I2,A)') '        <DataItem Format="XML" Dimensions="3">',0,0,0,'</DataItem>'
    write(nUnitFile,'(A,3ES15.7,A)') '        <DataItem Format="XML" Dimensions="3">',dx,dy,dz,'</DataItem>'
    write(nUnitFile,'(A)')'      </Geometry>'
    write(nUnitFile,'(A)')'      <Time Value="0"/>'

    ! add virables
    out_file_h5 = trim(adjustl(outname))//".pdr.h5"
    call Write_XDMF_One(nUnitFile,out_file_h5,'tgas')
    call Write_XDMF_One(nUnitFile,out_file_h5,'rho')
    call Write_XDMF_One(nUnitFile,out_file_h5,'uv')
    call Write_XDMF_One(nUnitFile,out_file_h5,'abundance033')
    call Write_XDMF_One(nUnitFile,out_file_h5,'av011')

    out_file_h5 = trim(adjustl(outname))//".cool.h5"
    call Write_XDMF_One(nUnitFile,out_file_h5,'cool004')
    call Write_XDMF_One(nUnitFile,out_file_h5,'totalcool')

    out_file_h5 = trim(adjustl(outname))//".heat.h5"
    call Write_XDMF_One(nUnitFile,out_file_h5,'heat012')

    write(nUnitFile,'(A)')'    </Grid>'
    write(nUnitFile,'(A)')'  </Domain>'
    write(nUnitFile,'(A)')'</Xdmf>'
    close(nUnitFile,IOSTAT=ierror)
  end subroutine initvisual

  !******************************************************************
  ! Write_XDMF_One
  !******************************************************************
  subroutine Write_XDMF_One(nUnitFile, outname, chAttribute)
    implicit none
    integer,intent(in)::nUnitFile
    character(*),intent(in)::outname, chAttribute

    ! locals
    character(128)::chFile
    integer::iprec=mytype_save

    write(chFile,'(3A)')trim(outname)//':/'//trim(chAttribute)
    write(nUnitFile,'(A)')'      <Attribute Name="'//trim(chAttribute)//'" AttributeType="Scalar" Center="Node">'
    write(nUnitFile,'(A,I1,A,3I7,A)')'        <DataItem Format="HDF" DataType="Float" Precision="',&
    &iprec,'" Dimensions="',nxc,nyc,nzc,'">'
    write(nUnitFile,'(A)')'          '//trim(chFile)
    write(nUnitFile,'(A)')'        </DataItem>'
    write(nUnitFile,'(A)')'      </Attribute>'
  end subroutine Write_XDMF_One

end module m_visual
