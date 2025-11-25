program hdf5_test
    use hdf5
    implicit none

    integer :: error
    integer(hid_t) :: file_id, dset_id, dspace_id
    integer(hsize_t), dimension(2) :: dims = [3, 3]
    real, dimension(3, 3) :: data = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 3])

    ! init HDF5
    call h5open_f(error)

    ! creat a new file
    call h5fcreate_f("test.h5", H5F_ACC_TRUNC_F, file_id, error)

    ! create data space
    call h5screate_simple_f(2, dims, dspace_id, error)

    ! create dataset
    call h5dcreate_f(file_id, "dataset", H5T_NATIVE_REAL, dspace_id, dset_id, error)

    ! write data
    call h5dwrite_f(dset_id, H5T_NATIVE_REAL, data, dims, error)

    ! close resource 
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)
    call h5fclose_f(file_id, error)

    ! close HDF5
    call h5close_f(error)

    print *, "HDF5 create successfully！"
end program hdf5_test