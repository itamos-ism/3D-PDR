module m_writeoutputs
  use hdf5
  use m_paramters
  implicit none
  
  integer :: h5err
  integer :: coords(3)
  character(len=256) :: fin_file, out_file_h5
  integer(HID_T) :: file_id, filespace, filespace2, filespace3
  integer(HID_T) :: dset_id_network,dset_id_cooarray,dset_id_allpop,dtype_id,arrtype_id
  integer(HSIZE_T), dimension(3) :: dims
  integer(HSIZE_T), dimension(4) :: dims3 
  real(RK), allocatable :: temp_real(:,:,:)
  contains
    
  subroutine writeoutputs
    implicit none
    
    !locals
    integer :: i, j, k, isp, ipix, ico, ilevel, offset, file_unit
    integer :: access_flags
    real(RK), allocatable :: allpop(:,:,:,:)
    integer(hsize_t), dimension(1) :: dims_network, dims_cooarray
    character(len=128) :: network_name
    character(len=128), allocatable :: coo_array(:)
    character(len=128), DIMENSION(1) :: string_array
    character(len=128) :: file_name, file_index
    integer :: total_levels
    integer(8) :: str_size
    logical :: file_exists
    
    dims = (/ nxc,nyc,nzc /)

    allocate(temp_real(nxc, nyc, nzc))
    
    call h5open_f(h5err)
    
    !-------------------------------------
    !outname FOR ABUNDANCES AND TEMPERATURE
    !-------------------------------------
    fin_file = "../"//trim(adjustl(outname))//".pdr.fin"
    out_file_h5 = trim(adjustl(outname))//".pdr.h5"
    inquire(file=trim(fin_file), exist=file_exists)
    if (file_exists) then
      access_flags = H5F_ACC_TRUNC_F
      call h5fcreate_f(out_file_h5, access_flags, file_id, h5err)

      call write_h5_one_r3d(file_id, dims, pdr%x, "x")
      call write_h5_one_r3d(file_id, dims, pdr%y, "y")
      call write_h5_one_r3d(file_id, dims, pdr%z, "z")
      call write_h5_one_r3d(file_id, dims, pdr%Tgas, "tgas")
      call write_h5_one_r3d(file_id, dims, pdr%Tdust, "tdust")
      call write_h5_one_r3d(file_id, dims, DBLE(pdr%etype), "etype")
      call write_h5_one_r3d(file_id, dims, pdr%rho, "rho")
      call write_h5_one_r3d(file_id, dims, pdr%UVfield, "uv")

      do isp = 1, nspec
        write(file_index, '(I3.3)') isp
        file_name = "abundance"//trim(adjustl(file_index))
        do k=1,nzc
        do j=1,nyc
        do i=1,nxc
          temp_real(i,j,k) = pdr(i,j,k)%abundance(isp)
        enddo
        enddo
        enddo  
        call write_h5_one_r3d(file_id, dims, temp_real, file_name)
      enddo

      do ipix = 0, nrays-1
        write(file_index, '(I3.3)') ipix
        file_name = "av"//trim(adjustl(file_index))
        do k=1,nzc
        do j=1,nyc
        do i=1,nxc
          temp_real(i,j,k) = pdr(i,j,k)%AV(ipix)
        enddo
        enddo
        enddo
        call write_h5_one_r3d(file_id, dims, temp_real, file_name)
      enddo

      call h5fclose_f(file_id, h5err)
    endif
    !-----------------------------------------
    !END outname FOR ABUNDANCES AND TEMPERATURE
    !-----------------------------------------
    
    !-------------------------
    !outname FOR LOCAL CR VALUE
    !-------------------------
    fin_file = "../"//trim(adjustl(outname))//".cr.fin"
    out_file_h5 = trim(adjustl(outname))//".cr.h5"
    inquire(file=trim(fin_file), exist=file_exists)
    if (file_exists) then
      access_flags = H5F_ACC_TRUNC_F
      call h5fcreate_f(out_file_h5, access_flags, file_id, h5err)

      call write_h5_one_r3d(file_id, dims, pdr%zetalocal, "zetalocal")

      call h5fclose_f(file_id, h5err)
    endif
    !-----------------------------
    !END outname FOR LOCAL CR VALUE
    !-----------------------------
    
    !---------------------------
    !outname FOR COOLING FUNCTION
    !---------------------------
    fin_file = "../"//trim(adjustl(outname))//trim(adjustl(".cool"))//".fin"
    out_file_h5 = trim(adjustl(outname))//".cool.h5"
    inquire(file=trim(fin_file), exist=file_exists)
    if (file_exists) then
      access_flags = H5F_ACC_TRUNC_F
      call h5fcreate_f(out_file_h5, access_flags, file_id, h5err)

      do ico = 1, coo
        write(file_index, '(I3.3)') ico
        file_name = "cool"//trim(adjustl(file_index))
        do k=1,nzc
        do j=1,nyc
        do i=1,nxc
          temp_real(i,j,k) = pdr(i,j,k)%cooling(ico)
        enddo
        enddo
        enddo  
        call write_h5_one_r3d(file_id, dims, temp_real, file_name)
      enddo

      call write_h5_one_r3d(file_id, dims, pdr%totalcooling, "totalcool")

      call h5fclose_f(file_id, h5err)
    endif
    !-------------------------------
    !END outname FOR COOLING FUNCTION
    !-------------------------------
    
    !---------------------------
    !outname FOR HEATING FUNCTION
    !---------------------------
    fin_file = "../"//trim(adjustl(outname))//trim(adjustl(".heat"))//".fin"
    out_file_h5 = trim(adjustl(outname))//".heat.h5"
    inquire(file=trim(fin_file), exist=file_exists)
    if (file_exists) then
      access_flags = H5F_ACC_TRUNC_F
      call h5fcreate_f(out_file_h5, access_flags, file_id, h5err)

      do ico = 1, 12
        write(file_index, '(I3.3)') ico
        file_name = "heat"//trim(adjustl(file_index))
        do k=1,nzc
        do j=1,nyc
        do i=1,nxc
          temp_real(i,j,k) = pdr(i,j,k)%heating(ico)
        enddo
        enddo
        enddo  
        call write_h5_one_r3d(file_id, dims, temp_real, file_name)
      enddo

      call h5fclose_f(file_id, h5err)
    endif
    !-------------------------------
    !END outname FOR HEATING FUNCTION
    !-------------------------------
    
    !----------------------------
    !outname FOR RTtool
    !----------------------------
    fin_file = "../"//trim(adjustl(outname))//trim(adjustl(".RTspop"))//".fin"
    out_file_h5 = trim(adjustl(outname))//trim(adjustl(".RTspop"))//".h5"
    inquire(file=trim(fin_file), exist=file_exists)
    if (file_exists) then
      access_flags = H5F_ACC_TRUNC_F
      call h5fcreate_f(out_file_h5, access_flags, file_id, h5err)
      
      allocate(coo_array(coo))
      network_name = coolfile(1)
      string_array(1)=adjustl(network_name)
      do ico=1,coo
        coo_array(ico) = adjustl(coolfile(ico + 1))
      enddo

      total_levels = 0
      do k = 1, coo
        total_levels = total_levels + cnlev(k)
      end do
      
      allocate(allpop(nxc, nyc, nzc, total_levels))

      do k=1,nzc
      do j=1,nyc
      do i=1,nxc
        offset = 0 
        do ico = 1, coo
          do ilevel=1,cnlev(ico)
            if (pdr(i,j,k)%coolant(ico)%pop(ilevel).lt.1d-99) pdr(i,j,k)%coolant(ico)%pop(ilevel) = 0.0D0
          end do
          allpop(i, j, k, offset + 1 : offset + cnlev(ico)) = pdr(i,j,k)%coolant(ico)%pop(1:cnlev(ico))
          offset = offset + cnlev(ico)
        end do
      enddo
      enddo
      enddo
      
      dims_network(1) = 1
      dims_cooarray(1) = coo
      dims3= (/nxc,nyc,nzc, total_levels/)

      call h5screate_simple_f(4, dims3, filespace, h5err)
      call h5screate_simple_f(1, dims_network, filespace2, h5err)
      call h5screate_simple_f(1, dims_cooarray, filespace3, h5err)

      call h5tcopy_f(H5T_C_S1, dtype_id, h5err)
      call h5tcopy_f(H5T_C_S1, arrtype_id, h5err)
      str_size = 128
      call h5tset_size_f(dtype_id, str_size, h5err)
      call h5tset_size_f(arrtype_id, str_size, h5err)

      call h5dcreate_f(file_id, "network", dtype_id, filespace2, dset_id_network, h5err)
      call h5dcreate_f(file_id, "coolfiles", arrtype_id, filespace3, dset_id_cooarray, h5err)
      call h5dcreate_f(file_id, "allpop", H5T_IEEE_F64LE, filespace, dset_id_allpop, h5err)

      call h5dwrite_f(dset_id_network, dtype_id, string_array, dims_network, h5err)
      call h5dwrite_f(dset_id_cooarray, arrtype_id, coo_array, dims_cooarray, h5err)
      call h5dwrite_f(dset_id_allpop, H5T_IEEE_F64LE, allpop, dims3, h5err)
      
      call h5sclose_f(filespace, h5err)
      call h5sclose_f(filespace2, h5err)
      call h5sclose_f(filespace3, h5err)
      call h5dclose_f(dset_id_network, h5err)
      call h5dclose_f(dset_id_cooarray, h5err)
      call h5dclose_f(dset_id_allpop, h5err)
      call h5fclose_f(file_id, h5err)

      deallocate(allpop)
    endif
    !--------------------------------
    !END outname FOR LEVEL POPULATIONS
    !--------------------------------

    !---------------------------
    !outname FOR VELOCITY
    !---------------------------
    fin_file = "../"//trim(adjustl(outname))//trim(adjustl(".velocity"))//".dat"
    out_file_h5 = trim(adjustl(outname))//".velocity.h5"
    inquire(file=trim(fin_file), exist=file_exists)
    if (file_exists) then
      access_flags = H5F_ACC_TRUNC_F
      call h5fcreate_f(out_file_h5, access_flags, file_id, h5err)

      call write_h5_one_r3d(file_id, dims, pdr%vx, "vx")
      call write_h5_one_r3d(file_id, dims, pdr%vy, "vy")
      call write_h5_one_r3d(file_id, dims, pdr%vz, "vz")

      call h5fclose_f(file_id, h5err)
    endif
    !--------------------------------
    !END outname FOR VELOCITY
    !--------------------------------

    deallocate(temp_real)
    call h5close_f(h5err)
  end subroutine writeoutputs

  subroutine write_h5_one_r3d(file_id, dims, temp_real, name)
    implicit none
    integer(HID_T), intent(in) :: file_id
    integer(HSIZE_T), dimension(3), intent(in) :: dims
    real(RK), intent(in) :: temp_real(:,:,:)
    character(len=*), intent(in) :: name

    !locals
    integer(HID_T) :: filespace, dset_id
    integer :: h5err

    call h5screate_simple_f(3, dims, filespace, h5err)
    call h5dcreate_f(file_id, trim(name), H5T_IEEE_F64LE, filespace, dset_id, h5err)

    call h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp_real, dims, h5err)

    if (h5err /= 0) print *, "h5dwrite_f failed with error:", name, h5err

    call h5sclose_f(filespace, h5err)
    call h5dclose_f(dset_id, h5err)
  end subroutine write_h5_one_r3d
    
end module m_writeoutputs