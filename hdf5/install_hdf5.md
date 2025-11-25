# install and use HDF5

## 1 install

1. Download the installation package and decompression it . link :  https://support.hdfgroup.org/downloads/hdf5/hdf5_1_14_6.html

   ```shell
   wget https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_6/downloads/hdf5-1.14.6.tar.gz
   tar -xzvf hdf5-1.14.6.tar.gz 
   ```

2. chose the right gcc

   ```shell
   module load compiler/gcc/7.3.1
   ```

3. go to the folder of decompression and configure: 

   ```shell
   cd hdf5-1.14.5
   ./configure CC=gcc FC=gfortran --prefix=yourpath/hdf5-serial --enable-fortran
   ```

   - The path after prefix is the path you want to install. 
   - Remember to add -- enable fortran to use HDF5 in fortran

4. make the HDF5

   ```shell
   make clean
   make -j4
   make install
   ```

5. add HDF5 to the environment variable

   ```shell
   export LD_LIBRARY_PATH=yourpath/hdf5-serial/lib:$LD_LIBRARY_PATH
   export PATH=yourpath/hdf5-serial/bin:$PATH
   ```

6. test

   create a file : `hdf5_test.f90`

   ```fortran
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
   
   ```

   compile 

   ```shell
   gfortran -o test hdf5_test.f90   -I/yourpath/hdf5-serial/include   -L/yourpath/hdf5-serial/lib   -lhdf5_fortran -lhdf5
   
   ```

   run  

   ```
    ./test
   ```

   you can get a HDF5 file `test.h5` and use `h5dump` see it 

   ```shell
    h5dump test.h5
   ```

   you can see those:

   ```shell
   HDF5 "test.h5" {
   GROUP "/" {
      DATASET "dataset" {
         DATATYPE  H5T_IEEE_F32LE
         DATASPACE  SIMPLE { ( 3, 3 ) / ( 3, 3 ) }
         DATA {
         (0,0): 1, 2, 3,
         (1,0): 4, 5, 6,
         (2,0): 7, 8, 9
         }
      }
   }
   }
   ```