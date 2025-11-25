subroutine writeoutputs
!for ASCII outputs
use hdf5
use healpix_types
use maincode_module
use global_module

integer :: error
integer(hid_t) :: file_id,coolfile_id,heatfile_id, dset_id, dspace_id, x_id, xset_id,av6_id, av6set_id, y_id, yset_id, z_id, zset_id
integer(hid_t) :: tgas_id, tgasset_id, tdust_id, tdustset_id, etype_id, etypeset_id, rho_id, rhoset_id, uvfield_id, uvfieldset_id
integer(hid_t) :: abundance_id, abundanceset_id, av_id, avset_id,cool_id, coolset_id,totalcooling_id,totalcoolingset_id,dtype_id
integer(hid_t) :: heat_id, heatset_id, line_id, lineset_id, linefile_id, spopfile_id, rtspopfile_id, network_id, networkset_id
integer(hid_t) :: cooarray_id,  cooarrayset_id,arrtype_id
integer(hsize_t), dimension(2) :: dims_abundance
integer(hsize_t), dimension(2) :: dims_av
integer(hsize_t), dimension(2) :: dims_cool
integer(hsize_t), dimension(2) :: dims_heat
integer(hsize_t), dimension(2) :: dims_x
integer(hsize_t), dimension(2) :: dims_line
integer(hsize_t), dimension(1) :: dims_network
integer(hsize_t), dimension(1) :: dims_cooarray
real(kind=dp), dimension(pdr_ptot, nspec) :: abundance_data
real(kind=dp), dimension(pdr_ptot, nrays) :: av_data
real(kind=dp), dimension(pdr_ptot, coo) :: cool_data
real(kind=dp), dimension(pdr_ptot, coo) :: coolant_data
real(kind=dp), dimension(pdr_ptot, 12) :: heat_data
real(kind=dp), dimension(pdr_ptot, 1) :: x_data, av6_data, y_data, z_data, tgas_data, tdust_data,rho_data,uvfield_data
real(kind=dp), dimension(pdr_ptot, 1) :: totalcooling_data
real(kind=dp), allocatable :: line_data(:,:) 
integer(kind=i4b) , dimension(pdr_ptot, nspec) :: etype_data
real(kind=dp), allocatable :: write_buffer(:)
character(len=128), DIMENSION(1) :: string_array
integer(8) :: str_size

character(len=20), DIMENSION(coo) :: coo_array
integer :: num_levels, ilevel, total_levels
integer(8) :: var_size

! real, dimension(3, 3) :: data = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 3])
character(len=256) :: out_file_h5
character(len=128) :: network_name
integer :: access_flags,offset

dims_network(1) = 1

dims_abundance(1) = pdr_ptot
dims_abundance(2) = nspec

dims_av(1) = pdr_ptot
dims_av(2) = nrays

dims_cool(1) = pdr_ptot
dims_cool(2) = coo

dims_heat(1) = pdr_ptot
dims_heat(2) = 12

dims_x(1) = pdr_ptot
dims_x(2) = 1

dims_cooarray(1) = 1


#ifdef PYWRAP
out_file = trim(adjustl(output))//".pywrap"
out_file2 = trim(adjustl(out_file))//"]"
write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
open(unit=21,file=out_file,status='replace')

write(21,'("x Av Tgas Tdust nH fuv")',advance='no')
do i=1,nspec
  write(21,'(A6)',advance='no') trim(species(i))
enddo
write(21,'(" dummy1 photoelectric dummy2 Cionization H2formation H2photodis FUVpumping")',advance='no')
write(21,'(" cosmicray turbulent chemical gasgrain totheat ")',advance='no')
do i=1,coo
   write(21,'(A10)',advance='no') trim(adjustl(coolant(i)%cname))//"cool"
enddo
write(21,'(" totcool")',advance='no')
write(21,*) ''

do p=1,pdr_ptot
     write(21,'(300ES15.7)') pdr(p)%x, pdr(p)%AV(6), pdr(p)%Tgas,pdr(p)%Tdust,pdr(p)%rho,pdr(p)%UVfield,&
        &pdr(p)%abundance,pdr(p)%heating,pdr(p)%cooling,pdr(p)%totalcooling
enddo

close(21)
#endif

!-------------------------------------
!OUTPUT FOR CHEMICAL ANALYSIS
!-------------------------------------
 
!   do p=1,pdr_ptot
!      call analyse_chemistry(p, end_time, pdr(p)%rho, pdr(p)%Tgas, &
!        &12, pdr(p)%AV(6), nspec, species,pdr(p)%abundance(1:nspec),nreac, reactant, &
!        & product, temp_rate(:,p))
!   enddo


!-------------------------------------
!OUTPUT FOR PDR PARAMETERS
!-------------------------------------
   out_file = trim(adjustl(output))//".params"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
   open(unit=21,file=out_file,status='replace')
   write(21,'(3ES11.3)') Gext(1), zeta*1.3d-17, metallicity
   close(21)

  !  Prepare output data
  do p=1,pdr_ptot
    x_data(p, :) = pdr(p)%x
    av6_data(p, :) = pdr(p)%AV(6)
    y_data(p, :) = pdr(p)%y
    z_data(p, :) = pdr(p)%z
    tgas_data(p, :) = pdr(p)%Tgas
    tdust_data(p, :) = pdr(p)%Tdust
    etype_data(p, :) = pdr(p)%etype
    rho_data(p, :) = pdr(p)%rho
    uvfield_data(p, :) = pdr(p)%UVfield
    abundance_data(p, :) = pdr(p)%abundance
    av_data(p, :) = pdr(p)%AV
    cool_data(p, :) = pdr(p)%cooling
    totalcooling_data(p, :) = pdr(p)%totalcooling
    heat_data(p, :) = pdr(p)%heating

  enddo

   call h5open_f(error)
!-------------------------------------
!OUTPUT FOR ABUNDANCES AND TEMPERATURE
!-------------------------------------
  !  out_file = trim(adjustl(output))//".pdr.fin"
  !  out_file2 = trim(adjustl(out_file))//"]"
  !  write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
  !  open(unit=21,file=out_file,status='replace')
   out_file_h5 = trim(adjustl(output))//".pdr.h5"
   out_file2 = trim(adjustl(out_file_h5))//"]"
   write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
   access_flags = H5F_ACC_TRUNC_F 
   call h5fcreate_f(out_file_h5, access_flags, file_id, error)
   if (error /= 0) stop "Error creating HDF5 file"
   
   call h5screate_simple_f(2, dims_x, x_id, error)
   call h5screate_simple_f(2, dims_x, tgas_id, error)
   call h5screate_simple_f(2, dims_x, tdust_id, error)
   call h5screate_simple_f(2, dims_x, etype_id, error)  
   call h5screate_simple_f(2, dims_x, rho_id, error)
   call h5screate_simple_f(2, dims_x, uvfield_id, error)
   call h5screate_simple_f(2, dims_abundance, abundance_id, error)
   
  !  call h5screate_simple_f(2, dims_simple, dspace_id, error)
   if (error /= 0) stop "Error creating data space"
  !  call h5dcreate_f(file_id, "abundance", H5T_IEEE_F64LE, dspace_id, dset_id, error)
   call h5dcreate_f(file_id, "x", H5T_IEEE_F64LE, x_id, xset_id, error)
   call h5dcreate_f(file_id, "tgas", H5T_IEEE_F64LE, tgas_id, tgasset_id, error)
   call h5dcreate_f(file_id, "tdust", H5T_IEEE_F64LE, tdust_id, tdustset_id, error)
   call h5dcreate_f(file_id, "etype", H5T_STD_I32LE, etype_id, etypeset_id, error)
   call h5dcreate_f(file_id, "rho", H5T_IEEE_F64LE, rho_id, rhoset_id, error)
   call h5dcreate_f(file_id, "uvfield", H5T_IEEE_F64LE, uvfield_id, uvfieldset_id, error)
   call h5dcreate_f(file_id, "abundance", H5T_IEEE_F64LE, abundance_id, abundanceset_id, error)
   
   if (error /= 0) stop "Error creating dataset"

  ! do p=1,pdr_ptot
#ifdef ONEDIMENSIONAL
    !  write(21,'(I7,4ES15.7,I5,300ES15.7)') p,pdr(p)%x, pdr(p)%AV(6), pdr(p)%Tgas,pdr(p)%Tdust,pdr(p)%etype,&
    !  &pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance
     call h5screate_simple_f(2, dims_x, av6_id, error)
     call h5dcreate_f(file_id, "av6", H5T_IEEE_F64LE, av6_id, av6set_id, error)

     call h5dwrite_f(xset_id, H5T_IEEE_F64LE, x_data, dims_x, error)
     call h5dwrite_f(av6set_id, H5T_IEEE_F64LE, av6_data, dims_x, error)
     call h5dwrite_f(tgasset_id, H5T_IEEE_F64LE, tgas_data, dims_x, error)
     call h5dwrite_f(tdustset_id, H5T_IEEE_F64LE, tdust_data, dims_x, error)
     call h5dwrite_f(etypeset_id, H5T_STD_I32LE, etype_data, dims_x, error)
     call h5dwrite_f(rhoset_id, H5T_IEEE_F64LE, rho_data, dims_x, error)
     call h5dwrite_f(uvfieldset_id, H5T_IEEE_F64LE, uvfield_data, dims_x, error)
     call h5dwrite_f(abundanceset_id, H5T_IEEE_F64LE, abundance_data, dims_abundance, error)

     call h5dclose_f(av6set_id, error)
#else
    !  write(21,'(I7,5ES15.7,I5,300ES15.7)') p,pdr(p)%x, pdr(p)%y, pdr(p)%z, pdr(p)%Tgas,pdr(p)%Tdust,&
    !  &pdr(p)%etype,pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance,pdr(p)%AV
     call h5screate_simple_f(2, dims_x, y_id, error)
     call h5dcreate_f(file_id, "y", H5T_IEEE_F64LE, y_id, yset_id, error)
     call h5screate_simple_f(2, dims_x, z_id, error)
     call h5dcreate_f(file_id, "z", H5T_IEEE_F64LE, z_id, zset_id, error)
     call h5screate_simple_f(2, dims_av, av_id, error)
     call h5dcreate_f(file_id, "av", H5T_IEEE_F64LE, av_id, avset_id, error)

     call h5dwrite_f(xset_id, H5T_IEEE_F64LE, x_data, dims_x, error)
     call h5dwrite_f(yset_id, H5T_IEEE_F64LE, y_data, dims_x, error)
     call h5dwrite_f(zset_id, H5T_IEEE_F64LE, z_data, dims_x, error)
     call h5dwrite_f(tgasset_id, H5T_IEEE_F64LE, tgas_data, dims_x, error)
     call h5dwrite_f(tdustset_id, H5T_IEEE_F64LE, tdust_data, dims_x, error)
     call h5dwrite_f(etypeset_id, H5T_STD_I32LE, etype_data, dims_x, error)
     call h5dwrite_f(rhoset_id, H5T_IEEE_F64LE, rho_data, dims_x, error)
     call h5dwrite_f(uvfieldset_id, H5T_IEEE_F64LE, uvfield_data, dims_x, error)
     call h5dwrite_f(abundanceset_id, H5T_IEEE_F64LE, abundance_data, dims_abundance, error)
     call h5dwrite_f(avset_id, H5T_IEEE_F64LE, av_data, dims_av, error)

     call h5dclose_f(yset_id, error)
     call h5dclose_f(zset_id, error)
     call h5dclose_f(avset_id, error)
#endif
  ! enddo
  ! call h5dclose_f(dset_id, error)
  
  call h5dclose_f(xset_id, error)
  call h5dclose_f(tgasset_id, error)
  call h5dclose_f(tdustset_id, error)
  call h5dclose_f(etypeset_id, error)
  call h5dclose_f(rhoset_id, error)
  call h5dclose_f(uvfieldset_id, error)
  call h5dclose_f(abundanceset_id, error)

  if (error /= 0) stop "Error closing dataset"

  call h5fclose_f(file_id, error)
  if (error /= 0) stop "Error closing HDF5 file"

  print *, "Data written to HDF5 file successfully!"

! close(21)
!-----------------------------------------
!END OUTPUT FOR ABUNDANCES AND TEMPERATURE
!-----------------------------------------


!---------------------------
!OUTPUT FOR COOLING FUNCTION
!---------------------------
  !  out_file = trim(adjustl(output))//trim(adjustl(".cool"))//".fin"
  !  out_file2 = trim(adjustl(out_file))//"]"
  !  write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
  !  open(unit=13,file=out_file,status='replace')
  out_file_h5 = trim(adjustl(output))//trim(adjustl(".cool"))//".h5"
  out_file2 = trim(adjustl(out_file_h5))//"]"
  write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
  access_flags = H5F_ACC_TRUNC_F 
  call h5fcreate_f(out_file_h5, access_flags, coolfile_id, error)
  if (error /= 0) stop "Error creating HDF5 file"

   
  call h5screate_simple_f(2, dims_x, x_id, error)
  call h5screate_simple_f(2, dims_cool, cool_id, error)
  call h5screate_simple_f(2, dims_x, totalcooling_id, error)
  call h5dcreate_f(coolfile_id, "x", H5T_IEEE_F64LE, x_id, xset_id, error)
  call h5dcreate_f(coolfile_id, "cooling", H5T_IEEE_F64LE, cool_id, coolset_id, error)
  call h5dcreate_f(coolfile_id, "totalcooling", H5T_IEEE_F64LE, totalcooling_id, totalcoolingset_id, error)
  !  do p=1,pdr_ptot
#ifdef ONEDIMENSIONAL
      ! write(13,'(I7,200ES15.7)') p, pdr(p)%x, pdr(p)%AV(6), pdr(p)%cooling(:),pdr(p)%totalcooling
    call h5screate_simple_f(2, dims_x, av6_id, error)
    call h5dcreate_f(coolfile_id, "av6", H5T_IEEE_F64LE, av6_id, av6set_id, error)

    call h5dwrite_f(xset_id, H5T_IEEE_F64LE, x_data, dims_x, error)
    call h5dwrite_f(av6set_id, H5T_IEEE_F64LE, av6_data, dims_x, error)
    call h5dwrite_f(coolset_id, H5T_IEEE_F64LE, cool_data, dims_cool, error)
    call h5dwrite_f(totalcoolingset_id, H5T_IEEE_F64LE, totalcooling_data, dims_x, error)

    call h5dclose_f(av6set_id, error)
     
#else
      ! write(13,'(I7,200ES15.7)') p, pdr(p)%x, pdr(p)%y, pdr(p)%z, pdr(p)%cooling(:),pdr(p)%totalcooling, pdr(p)%AV(:)
    call h5screate_simple_f(2, dims_x, y_id, error)
    call h5screate_simple_f(2, dims_x, z_id, error)
    call h5screate_simple_f(2, dims_av, av_id, error)
  
    call h5dcreate_f(coolfile_id, "y", H5T_IEEE_F64LE, y_id, yset_id, error)
    call h5dcreate_f(coolfile_id, "z", H5T_IEEE_F64LE, z_id, zset_id, error)
    call h5dcreate_f(coolfile_id, "av", H5T_IEEE_F64LE, av_id, avset_id, error)

    call h5dwrite_f(xset_id, H5T_IEEE_F64LE, x_data, dims_x, error)
    call h5dwrite_f(yset_id, H5T_IEEE_F64LE, y_data, dims_x, error)
    call h5dwrite_f(zset_id, H5T_IEEE_F64LE, z_data, dims_x, error)
    call h5dwrite_f(coolset_id, H5T_IEEE_F64LE, cool_data, dims_cool, error)
    call h5dwrite_f(totalcoolingset_id, H5T_IEEE_F64LE, totalcooling_data, dims_x, error)
    call h5dwrite_f(avset_id, H5T_IEEE_F64LE, av_data, dims_av, error)

    call h5dclose_f(yset_id, error)
    call h5dclose_f(zset_id, error)
    call h5dclose_f(avset_id, error)
#endif
  !  enddo
    call h5dclose_f(xset_id, error)
    call h5dclose_f(coolset_id, error)
    call h5dclose_f(totalcoolingset_id, error)
    call h5fclose_f(coolfile_id, error)
!-------------------------------
!END OUTPUT FOR COOLING FUNCTION
!-------------------------------


!---------------------------
!OUTPUT FOR HEATING FUNCTION
!---------------------------
  !  out_file = trim(adjustl(output))//trim(adjustl(".heat"))//".fin"
  !  out_file2 = trim(adjustl(out_file))//"]"
  !  write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
  !  open(unit=14,file=out_file,status='replace')

  out_file_h5 = trim(adjustl(output))//trim(adjustl(".heat"))//".h5"
  out_file2 = trim(adjustl(out_file_h5))//"]"
  write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
  access_flags = H5F_ACC_TRUNC_F 
  call h5fcreate_f(out_file_h5, access_flags, heatfile_id, error)
  if (error /= 0) stop "Error creating HDF5 file"

  call h5screate_simple_f(2, dims_x, x_id, error)
  call h5screate_simple_f(2, dims_heat, heat_id, error)
  call h5dcreate_f(heatfile_id, "x", H5T_IEEE_F64LE, x_id, xset_id, error)
  call h5dcreate_f(heatfile_id, "heating", H5T_IEEE_F64LE, heat_id, heatset_id, error)
  !  do p=1,pdr_ptot
#ifdef ONEDIMENSIONAL
      ! write(14,'(I7,200ES15.7)') p, pdr(p)%x, pdr(p)%AV(6), pdr(p)%heating
    call h5screate_simple_f(2, dims_x, av6_id, error)
    call h5dcreate_f(heatfile_id, "av6", H5T_IEEE_F64LE, av6_id, av6set_id, error)
    call h5dwrite_f(xset_id, H5T_IEEE_F64LE, x_data, dims_x, error)
    call h5dwrite_f(av6set_id, H5T_IEEE_F64LE, av6_data, dims_x, error)
    call h5dwrite_f(heatset_id, H5T_IEEE_F64LE, heat_data, dims_heat, error)
    call h5dclose_f(av6set_id, error)

#else
      ! write(14,'(I7,200ES15.7)') p, pdr(p)%x, pdr(p)%y, pdr(p)%z, pdr(p)%heating, pdr(p)%AV
    call h5screate_simple_f(2, dims_x, y_id, error)
    call h5screate_simple_f(2, dims_x, z_id, error)
    call h5screate_simple_f(2, dims_av, av_id, error)
    call h5dcreate_f(heatfile_id, "y", H5T_IEEE_F64LE, y_id, yset_id, error)
    call h5dcreate_f(heatfile_id, "z", H5T_IEEE_F64LE, z_id, zset_id, error)
    call h5dcreate_f(heatfile_id, "av", H5T_IEEE_F64LE, av_id, avset_id, error)
    call h5dwrite_f(xset_id, H5T_IEEE_F64LE, x_data, dims_x, error)
    call h5dwrite_f(yset_id, H5T_IEEE_F64LE, y_data, dims_x, error)
    call h5dwrite_f(zset_id, H5T_IEEE_F64LE, z_data, dims_x, error)
    call h5dwrite_f(heatset_id, H5T_IEEE_F64LE, heat_data, dims_heat, error)
    call h5dwrite_f(avset_id, H5T_IEEE_F64LE, av_data, dims_av, error)
    call h5dclose_f(yset_id, error)
    call h5dclose_f(zset_id, error)
    call h5dclose_f(avset_id, error)
#endif
  !  enddo
  !  close(14)
  call h5dclose_f(xset_id, error)
  call h5dclose_f(heatset_id, error)
  call h5fclose_f(heatfile_id, error)
!-------------------------------
!END OUTPUT FOR HEATING FUNCTION
!-------------------------------


#ifdef ONEDIMENSIONAL
!-----------------------
!OUTPUT FOR EMISSIVITIES
!-----------------------
 do k=1,coo
  !  out_file = trim(adjustl(output))//"."//trim(adjustl(coolant(k)%cname))//trim(adjustl(".line"))//".fin"
  !  out_file2 = trim(adjustl(out_file))//"]"
  !  write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
  !  open(unit=16,file=out_file,status='replace')
   
   out_file_h5 = trim(adjustl(output))//"."//trim(adjustl(coolant(k)%cname))//trim(adjustl(".line"))//".h5"
   out_file2 = trim(adjustl(out_file_h5))//"]"
   write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
   access_flags = H5F_ACC_TRUNC_F 
   call h5fcreate_f(out_file_h5, access_flags, linefile_id, error)
   dims_line(1) = pdr_ptot
   dims_line(2) = coolant(k)%cnlev - 1
   if (error /= 0) stop "Error creating HDF5 file"
   call h5screate_simple_f(2, dims_x, x_id, error)
   call h5screate_simple_f(2, dims_x, av6_id, error)
   call h5screate_simple_f(2, dims_line, line_id, error)
   
   call h5dcreate_f(linefile_id, "x", H5T_IEEE_F64LE, x_id, xset_id, error)
   call h5dcreate_f(linefile_id, "av6", H5T_IEEE_F64LE, av6_id, av6set_id, error)
   call h5dcreate_f(linefile_id, "line", H5T_IEEE_F64LE, line_id, lineset_id, error)
   
   num_levels = coolant(k)%cnlev - 1
   allocate(line_data(pdr_ptot, num_levels))
   allocate(write_buffer(num_levels))
   do p = 1, pdr_ptot
    do ilevel = 1, num_levels
     IF (pdr(p)%coolant(k)%line(ilevel+1, ilevel) .LT. 1D-99) THEN
       pdr(p)%coolant(k)%line(ilevel+1, ilevel) = 0.0D0
     endif
       write_buffer(ilevel) = pdr(p)%coolant(k)%line(ilevel+1, ilevel)
    end do
    line_data(p, :) = write_buffer
   end do
   
   call h5dwrite_f(xset_id, H5T_IEEE_F64LE, x_data, dims_x, error)
   call h5dwrite_f(av6set_id, H5T_IEEE_F64LE, av6_data, dims_x, error) 
   call h5dwrite_f(lineset_id, H5T_IEEE_F64LE, line_data, dims_line, error)

   call h5dclose_f(xset_id, error)
   call h5dclose_f(av6set_id, error)
   call h5dclose_f(lineset_id, error)
   deallocate(write_buffer)
   deallocate(line_data)
   !experimental [1D models only]
  !  do p = 1, pdr_ptot
  !   num_levels = coolant(k)%cnlev - 1
  !   allocate(write_buffer(num_levels))
 
  !   do ilevel = 1, num_levels
  !    IF (pdr(p)%coolant(k)%line(ilevel+1, ilevel) .LT. 1D-99) THEN
  !      pdr(p)%coolant(k)%line(ilevel+1, ilevel) = 0.0D0
  !    endif
  !      write_buffer(ilevel) = pdr(p)%coolant(k)%line(ilevel+1, ilevel)
  !   end do
  !   write(16, '(I5, 1X, ES15.7, 1X, ES15.7, 1X)', advance='no') p,pdr(p)%x, pdr(p)%AV(6)
  !   write(16, '(100ES15.7)') write_buffer
  !   deallocate(write_buffer)
    ! write(16, *)
  !    write(16, '(I5, 1X, ES15.7, 1X, ES15.7, 1X)', advance='no') p,pdr(p)%x, pdr(p)%AV(6)
  !    do ilevel=1,coolant(k)%cnlev-1
  !      if (pdr(p)%coolant(k)%line(ilevel+1,ilevel).lt.1d-99) pdr(p)%coolant(k)%line(ilevel+1,ilevel)=0.0D0
  !      write(16, '(100ES15.7)', advance='no') pdr(p)%coolant(k)%line(ilevel+1,ilevel)
  !    end do
  !    write(16, *)
   end do
  !  close(16)
    
    ! close(16)
!  enddo
!---------------------------
!END OUTPUT FOR EMISSIVITIES
!---------------------------

!----------------------------
!OUTPUT FOR LEVEL POPULATIONS
!----------------------------
 do k=1,coo
  !  out_file = trim(adjustl(output))//"."//trim(adjustl(coolant(k)%cname))//trim(adjustl(".spop"))//".fin"
  !  out_file2 = trim(adjustl(out_file))//"]"
  !  write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
  !  open(unit=16,file=out_file,status='replace')
   out_file_h5 = trim(adjustl(output))//"."//trim(adjustl(coolant(k)%cname))//trim(adjustl(".spop"))//".h5"
   out_file2 = trim(adjustl(out_file_h5))//"]"
   write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
   access_flags = H5F_ACC_TRUNC_F 
   call h5fcreate_f(out_file_h5, access_flags, spopfile_id, error)
   dims_line(1) = pdr_ptot
   dims_line(2) = coolant(k)%cnlev
   if (error /= 0) stop "Error creating HDF5 file"
   call h5screate_simple_f(2, dims_x, x_id, error)
   call h5screate_simple_f(2, dims_x, av6_id, error)
   call h5screate_simple_f(2, dims_line, line_id, error)
   
   call h5dcreate_f(spopfile_id, "x", H5T_IEEE_F64LE, x_id, xset_id, error)
   call h5dcreate_f(spopfile_id, "av6", H5T_IEEE_F64LE, av6_id, av6set_id, error)
   call h5dcreate_f(spopfile_id, "spop", H5T_IEEE_F64LE, line_id, lineset_id, error)
   
   num_levels = coolant(k)%cnlev
   allocate(line_data(pdr_ptot, num_levels))

   !experimental [1D models only]
   do p = 1, pdr_ptot
    !  write(16, '(I5, 1X, ES15.7, 1X, ES15.7, 1X)', advance='no') p,pdr(p)%x, pdr(p)%AV(6)
     do ilevel=1,coolant(k)%cnlev
       if (pdr(p)%coolant(k)%pop(ilevel).lt.1d-99) pdr(p)%coolant(k)%pop(ilevel) = 0.0D0
     end do
    !  write(16, '(100ES15.7)', advance='no') pdr(p)%coolant(k)%pop
    !  write(16, *)
     line_data(p, :) = pdr(p)%coolant(k)%pop
   end do
  !  close(16)
   call h5dwrite_f(xset_id, H5T_IEEE_F64LE, x_data, dims_x, error)
   call h5dwrite_f(av6set_id, H5T_IEEE_F64LE, av6_data, dims_x, error) 
   call h5dwrite_f(lineset_id, H5T_IEEE_F64LE, line_data, dims_line, error)

   call h5dclose_f(xset_id, error)
   call h5dclose_f(av6set_id, error)
   call h5dclose_f(lineset_id, error)

   deallocate(line_data)
 enddo
!--------------------------------
!END OUTPUT FOR LEVEL POPULATIONS
!--------------------------------
#endif


!----------------------------
!OUTPUT FOR RTtool
!----------------------------
  !  out_file = trim(adjustl(output))//trim(adjustl(".RTspop"))//".fin"
  !  out_file2 = trim(adjustl(out_file))//"]"
  !  write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
  !  open(unit=16,file=out_file,status='replace')

   out_file_h5 = trim(adjustl(output))//trim(adjustl(".RTspop"))//".h5"
   out_file2 = trim(adjustl(out_file_h5))//"]"
   write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
   access_flags = H5F_ACC_TRUNC_F 
   call h5fcreate_f(out_file_h5, access_flags, rtspopfile_id, error)
   
#ifdef REDUCED
  !  write(16,*) 'REDUCED'
   network_name = 'REDUCED'
#endif
#ifdef MEDIUM
  !  write(16,*) 'MEDIUM'
   network_name = 'MEDIUM'
#endif
#ifdef FULL
  !  write(16,*) 'FULL'
   network_name = 'FULL'
#endif
   do k=1,coo
    !  write(16,*) coolfile(k)
     coo_array(k) = coolfile(k)
   enddo
   write(16,*) 'ENDCOOLFILES'
   string_array(1)=adjustl(network_name)
   call h5screate_simple_f(1, dims_network, network_id, error)
   call h5screate_simple_f(1, dims_cooarray, cooarray_id, error)

   call h5tcopy_f(H5T_C_S1, dtype_id, error)
   call h5tcopy_f(H5T_C_S1, arrtype_id, error)
   str_size = 32
   call h5tset_size_f(dtype_id, str_size, error)
   call h5tset_size_f(arrtype_id, str_size*coo, error)

   call h5dcreate_f(rtspopfile_id, "network", dtype_id, network_id, networkset_id, error)
   call h5dcreate_f(rtspopfile_id, "coolfiles", arrtype_id, cooarray_id, cooarrayset_id, error)
   
   call h5dwrite_f(networkset_id, dtype_id, string_array, dims_network, error)
   call h5dwrite_f(cooarrayset_id, arrtype_id, coo_array, dims_cooarray, error)

   call h5dclose_f(cooarrayset_id, error)
   call h5dclose_f(networkset_id, error)

   total_levels = 0
   do k = 1, coo
      total_levels = total_levels + coolant(k)%cnlev
   end do
    allocate(line_data(pdr_ptot, coo * total_levels))
    line_data = 0.0D0
    
   do p = 1, pdr_ptot
     offset = 0 
     do k = 1, coo
       do ilevel=1,coolant(k)%cnlev
         if (pdr(p)%coolant(k)%pop(ilevel).lt.1d-99) pdr(p)%coolant(k)%pop(ilevel) = 0.0D0
       end do
      line_data(p, offset + 1 : offset + coolant(k)%cnlev) = pdr(p)%coolant(k)%pop(1:coolant(k)%cnlev)
      offset = offset + coolant(k)%cnlev
     end do
    !  do k = 1, coo
    !    write(16, '(100ES15.7)', advance='no') pdr(p)%coolant(k)%pop(1:coolant(k)%cnlev)
    !  end do
    !  write(16, *)
   end do
   dims_line(1) = pdr_ptot
   dims_line(2) = total_levels
   call h5screate_simple_f(2, dims_line, line_id, error)
   
   call h5dcreate_f(rtspopfile_id, "allpop", H5T_IEEE_F64LE, line_id, lineset_id, error)
   call h5dwrite_f(lineset_id, H5T_IEEE_F64LE, line_data, dims_line, error)
   call h5dclose_f(lineset_id, error)
   deallocate(line_data)

  !  close(16)
!--------------------------------
!END OUTPUT FOR LEVEL POPULATIONS
!--------------------------------


return
end subroutine

