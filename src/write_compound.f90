! This is the F2003 version of the h5_compound.c example source code.
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! Copyright by the Board of Trustees of the University of Illinois.         *
! All rights reserved.                                                      *
!                                                                           *
! This file is part of HDF5.  The full HDF5 copyright notice, including     *
! terms governing use, modification, and redistribution, is contained in    *
! the files COPYING and Copyright.html.  COPYING can be found at the root   *
! of the source code distribution tree; Copyright.html can be found at the  *
! root level of an installed copy of the electronic HDF5 document set and   *
! is linked from the top-level documents page.  It can also be found at     *
! http://hdf.ncsa.uiuc.edu/HDF5/doc/Copyright.html.  If you do not have     *
! access to either file, you may request a copy from hdfhelp@ncsa.uiuc.edu. *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! This example shows how to create an array of a compound datatype which 
! contains an array of type complex and how to write it to hdf5 
! and how to read it back into a compound datatype for hdf5.
!

PROGRAM compound_complex_fortran2003

  USE hdf5
  USE ISO_C_BINDING
  IMPLICIT NONE

  INTEGER, PARAMETER :: r_k8 = KIND(0.0d0)
  INTEGER, PARAMETER :: NMAX = 3

  INTEGER(HID_T)   :: sample_type_id, dset_id, dspace_id, file_id
  INTEGER(HSIZE_T) :: dims(1) = (/NMAX/)
  INTEGER :: error

  COMPLEX(KIND=r_k8), DIMENSION(1:NMAX), TARGET :: samples, read_samples
  INTEGER :: i

  TYPE(C_PTR) :: f_ptr
  INTEGER(8) :: real_size, real_complex_size

  real_size = storage_size(1_r_k8, r_k8) / 8
  real_complex_size = real_size * 2_8  ! a complex is (real,real)

  ! Initialize data
  DO i=1,NMAX
     samples(1:NMAX) = (3.14159_r_k8, 2.71828_r_k8)
  END DO

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)

  ! Create a new file using default properties.
  CALL h5fcreate_f("test.h5", H5F_ACC_TRUNC_F, file_id, error)
  !
  ! Create the memory data type.
  !
  CALL H5Tcreate_f(H5T_COMPOUND_F, real_complex_size, sample_type_id, error)

  ! Create the array type
  ! Then use that array type to insert values into
  CALL H5Tinsert_f( sample_type_id, "r", &
       0_8, h5kind_to_type(r_k8,H5_REAL_KIND), error)
  CALL H5Tinsert_f( sample_type_id, "i", &
       real_size, h5kind_to_type(r_k8,H5_REAL_KIND), error)
  !
  ! Create dataspace
  !
  CALL h5screate_simple_f(1, dims, dspace_id, error)
  !
  ! Create the dataset.
  !
  CALL H5Dcreate_f(file_id, "samples",  sample_type_id, dspace_id, dset_id, error)
  !
  ! Write data to the dataset
  !
  f_ptr = C_LOC(samples(1))
  CALL H5Dwrite_f(dset_id, sample_type_id, f_ptr, error)
  ! Close up
  CALL h5dclose_f(dset_id, error)
  CALL h5sclose_f(dspace_id, error)
  CALL H5Tclose_f(sample_type_id, error)
  CALL h5fclose_f(file_id, error)
  !
  ! Open the file and the dataset.
  !
  CALL H5Fopen_f("test.h5", H5F_ACC_RDONLY_F, file_id, error)
  
  CALL H5Dopen_f(file_id, "samples", dset_id, error)
  !
  ! Create the memory data type.
  !
  CALL H5Tcreate_f(H5T_COMPOUND_F, 16_8, sample_type_id,error)

  CALL H5Tinsert_f( sample_type_id, "r", &
       0_8, h5kind_to_type(r_k8,H5_REAL_KIND), error)
  CALL H5Tinsert_f( sample_type_id, "i", &
       real_size, h5kind_to_type(r_k8,H5_REAL_KIND), error)

  f_ptr = C_LOC(read_samples(1))
  CALL H5Dread_f(dset_id, sample_type_id, f_ptr, error)

  !
  ! Display the fields
  !
  DO i=1,NMAX
     WRITE(*,'(A,3(" (",F8.5,",",F8.5,")"))') "SAMPLES =",read_samples(1:NMAX)
  END DO

  CALL H5Tclose_f(sample_type_id, error)
  CALL H5Dclose_f(dset_id, error)
  CALL H5Fclose_f(file_id, error)

END PROGRAM compound_complex_fortran2003
@Atefeh-FN
Comment

