!******************************************************************::
!
! The module contains functions for reading and writing data in HDF5 
! format
! 
!**********************************************************************
module io

use hdf5
use common

implicit none

contains 

subroutine write2file(A,fname)
complex(dp), intent(in) :: A(:,:)
CHARACTER(LEN=38), intent(in) :: fname
CHARACTER(LEN=38) :: filename
CHARACTER(LEN=3), PARAMETER :: dsetname1 = "A_r" ! Dataset name
CHARACTER(LEN=3), PARAMETER :: dsetname2 = "A_i" ! Dataset name

INTEGER(HID_T) :: file_id       ! File identifier
INTEGER(HID_T) :: dset_id1       ! Dataset identifier
INTEGER(HID_T) :: dset_id2       ! Dataset identifier
INTEGER(HID_T) :: dspace_id     ! Dataspace identifier


INTEGER(HSIZE_T), DIMENSION(2) :: dims  ! Dataset dimensions
INTEGER     ::    rank = 2                       ! Dataset rank

INTEGER     ::   error ! Error flag

filename = fname
dims = (/size(A,1),size(A,2)/)
     
CALL h5open_f(error)

!
! Create a new file using default properties.
!
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

!
! Create the dataspace.
!
CALL h5screate_simple_f(rank, dims, dspace_id, error)

!
! Create and write dataset using default properties.
!
CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, real(A), dims, error)

CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, imag(A), dims, error)

!
! End access to the dataset and release resources used by it.
!
CALL h5dclose_f(dset_id1, error)
CALL h5dclose_f(dset_id2, error)
!
! Terminate access to the data space.
!
CALL h5sclose_f(dspace_id, error)

!
! Close the file.
!
CALL h5fclose_f(file_id, error)

!
! Close FORTRAN interface.
!
CALL h5close_f(error)

end subroutine write2file

!_______________________________________________________________________

subroutine real_write2file(A,fname)
real(dp), intent(in) :: A(:,:)
CHARACTER(LEN=38), intent(in) :: fname
CHARACTER(LEN=38) :: filename 

CHARACTER(LEN=7), PARAMETER :: dsetname1 = "mueller" ! Dataset name


INTEGER(HID_T) :: file_id       ! File identifier
INTEGER(HID_T) :: dset_id1       ! Dataset identifier
INTEGER(HID_T) :: dspace_id     ! Dataspace identifier


INTEGER(HSIZE_T), DIMENSION(2) :: dims  ! Dataset dimensions
INTEGER     ::    rank = 2                       ! Dataset rank

INTEGER     ::   error ! Error flag
!INTEGER     :: i, j

dims = (/size(A,1),size(A,2)/)
filename = fname  
 
CALL h5open_f(error)

!
! Create a new file using default properties.
!
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

!
! Create the dataspace.
!
CALL h5screate_simple_f(rank, dims, dspace_id, error)

!
! Create and write dataset using default properties.
!
CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, A, dims, error)


!
! End access to the dataset and release resources used by it.
!
CALL h5dclose_f(dset_id1, error)

!
! Terminate access to the data space.
!
CALL h5sclose_f(dspace_id, error)

!
! Close the file.
!
CALL h5fclose_f(file_id, error)

!
! Close FORTRAN interface.
!
CALL h5close_f(error)

end subroutine real_write2file
!_---______________________________________________________________

subroutine cmplx_vec_write2file(vec,fname)
complex(dp), intent(in) :: vec(:)
CHARACTER(LEN=38), intent(in) :: fname
CHARACTER(LEN=38) :: filename

!CHARACTER(LEN=7), PARAMETER :: filename = "J.h5" ! File name
CHARACTER(LEN=3), PARAMETER :: dsetname1 = "J_r" ! Dataset name
CHARACTER(LEN=3), PARAMETER :: dsetname2 = "J_i" ! Dataset name
!INTEGER, PARAMETER :: NX = 1128

INTEGER(HID_T) :: file_id       ! File identifier
INTEGER(HID_T) :: dset_id1       ! Dataset identifier
INTEGER(HID_T) :: dset_id2       ! Dataset identifier
INTEGER(HID_T) :: dspace_id     ! Dataspace identifier


INTEGER(HSIZE_T), DIMENSION(1) :: dims           ! Dataset dimensions
INTEGER     ::    rank = 1                       ! Dataset rank

INTEGER     ::   error ! Error flag
!INTEGER     :: i, j

filename = fname
dims = (/size(vec)/)
     
CALL h5open_f(error)

! Create a new file using default properties.
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

! Create the dataspace.
CALL h5screate_simple_f(rank, dims, dspace_id, error)

! Create and write dataset using default properties.

CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, real(vec), dims, error)

CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, imag(vec), dims, error)


! End access to the dataset and release resources used by it.
CALL h5dclose_f(dset_id1, error)
CALL h5dclose_f(dset_id2, error)

! Terminate access to the data space.
CALL h5sclose_f(dspace_id, error)


! Close the file.
CALL h5fclose_f(file_id, error)

! Close FORTRAN interface.
CALL h5close_f(error)

end subroutine cmplx_vec_write2file

!________________________________________________________________________

subroutine real_vec_write2file(vec)
double precision, intent(in) :: vec(:)

CHARACTER(LEN=7), PARAMETER :: filename = "rcs.h5" ! File name
CHARACTER(LEN=3), PARAMETER :: dsetname1 = "rcs" ! Dataset name

INTEGER(HID_T) :: file_id       ! File identifier
INTEGER(HID_T) :: dset_id1       ! Dataset identifier
INTEGER(HID_T) :: dspace_id     ! Dataspace identifier


INTEGER(HSIZE_T), DIMENSION(1) :: dims     ! Dataset dimensions
INTEGER     ::    rank = 1                       ! Dataset rank

INTEGER     ::   error ! Error flag

dims  =  (/size(vec)/)

CALL h5open_f(error)

!
! Create a new file using default properties.
!
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

!
! Create the dataspace.
!
CALL h5screate_simple_f(rank, dims, dspace_id, error)

!
! Create and write dataset using default properties.
!
CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, vec, dims, error)

!
! End access to the dataset and release resources used by it.
!
CALL h5dclose_f(dset_id1, error)

!
! Terminate access to the data space.
!
CALL h5sclose_f(dspace_id, error)

!
! Close the file.
!
CALL h5fclose_f(file_id, error)

!
! Close FORTRAN interface.
!
CALL h5close_f(error)

end subroutine real_vec_write2file


subroutine read_field_points(coord)

CHARACTER(LEN=28), PARAMETER :: file = "field_points.h5" ! File name
CHARACTER(LEN=16), PARAMETER :: dataset1 = "coord"     ! Dataset name


integer(HID_T) :: file_id 
integer(HID_T) :: dataset1_id
integer(HID_T) :: dataspace_id
integer(HSIZE_T), dimension(2) :: coord_dims, dims_out

integer :: error

real(dp), dimension(:,:), allocatable :: coord

call h5open_f(error) ! initialize interface
call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error) ! open file


call h5dopen_f(file_id, dataset1, dataset1_id, error) ! open dataset

call h5dget_space_f(dataset1_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, coord_dims, error)

allocate(coord(coord_dims(1), coord_dims(2)))
call h5dread_f(dataset1_id, H5T_NATIVE_DOUBLE, coord, coord_dims, error)

call h5dclose_f(dataset1_id, error) ! close dataset

call h5fclose_f(file_id, error) ! close file
call h5close_f(error) ! close inteface


end subroutine read_field_points

!***********************************************************



subroutine read_geometry(file, coord, radius, param, angles, tind)
!CHARACTER(LEN=28), intent(in) :: fname

CHARACTER(LEN=38) :: file  ! File name
CHARACTER(LEN=16), PARAMETER :: dataset1 = "coord"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset2 = "radius"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset3 = "param_r"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset4 = "param_i"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset5 = "angles"  
CHARACTER(LEN=16), PARAMETER :: dataset6 = "tind" 

integer(HID_T) :: file_id 
integer(HID_T) :: dataset1_id, dataset2_id, dataset3_id, dataset4_id, &
     dataset5_id, dataset6_id
integer(HID_T) :: dataspace1_id, dataspace2_id, dataspace3_id, &
     dataspace4_id, dataspace5_id, dataspace6_id

integer(HSIZE_T), dimension(2) :: coord_dims, dims_out

integer :: error

real(dp), dimension(:,:), allocatable :: coord, angles
real(dp), dimension(:), allocatable :: radius, param_r, param_i
complex(dp), dimension(:), allocatable :: param
integer, dimension(:), allocatable :: tind

call h5open_f(error) ! initialize interface
call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error) ! open file

if(error == -1) then
   print*, 'Cannot find file:', file
   stop
end if

!____________Read coord _____________________________
call h5dopen_f(file_id, dataset1, dataset1_id, error) ! open dataset
if(error == -1) then
   print*, 'Cannot open dataset coord in file:', file
   stop
end if


call h5dget_space_f(dataset1_id, dataspace1_id, error) 
call H5sget_simple_extent_dims_f(dataspace1_id, dims_out, coord_dims, error)

allocate(coord(coord_dims(1), coord_dims(2)))
call h5dread_f(dataset1_id, H5T_NATIVE_DOUBLE, coord, coord_dims, error)

call h5dclose_f(dataset1_id, error) ! close dataset

!__________Read radius ________________________________
call h5dopen_f(file_id, dataset2, dataset2_id, error) ! open dataset
if(error == -1) then
   print*, 'Cannot open dataset radius in file:', file
   stop
end if


call h5dget_space_f(dataset2_id, dataspace2_id, error) 
call H5sget_simple_extent_dims_f(dataspace2_id, dims_out, coord_dims, error)

allocate(radius(coord_dims(2)))
call h5dread_f(dataset2_id, H5T_NATIVE_DOUBLE, radius, coord_dims, error)

call h5dclose_f(dataset2_id, error) ! close dataset

!__________Read param_r______________________________________

call h5dopen_f(file_id, dataset3, dataset3_id, error) ! open dataset
if(error == -1) then
   print*, 'Cannot open dataset param_r in file:', file
   stop
end if



call h5dget_space_f(dataset3_id, dataspace3_id, error) 
call H5sget_simple_extent_dims_f(dataspace3_id, dims_out, coord_dims, error)

allocate(param_r(coord_dims(2)))
call h5dread_f(dataset3_id, H5T_NATIVE_DOUBLE, param_r, coord_dims, error)

call h5dclose_f(dataset3_id, error) ! close dataset

!__________Read param_i______________________________________

call h5dopen_f(file_id, dataset4, dataset4_id, error) ! open dataset
if(error == -1) then
   print*, 'Cannot open dataset param_i in file:', file
   stop
end if



call h5dget_space_f(dataset4_id, dataspace4_id, error) 
call H5sget_simple_extent_dims_f(dataspace4_id, dims_out, coord_dims, error)

allocate(param_i(coord_dims(2)))
call h5dread_f(dataset4_id, H5T_NATIVE_DOUBLE, param_i, coord_dims, error)

call h5dclose_f(dataset4_id, error) ! close dataset

!______________Read euler angles___________________________________

call h5dopen_f(file_id, dataset5, dataset5_id, error) ! open dataset
if(error == -1) then
   print*, 'Cannot open dataset angles in file:', file
   stop
end if


call h5dget_space_f(dataset5_id, dataspace5_id, error) 
call H5sget_simple_extent_dims_f(dataspace5_id, dims_out, coord_dims, error)

allocate(angles(coord_dims(1), coord_dims(2)))

call h5dread_f(dataset5_id, H5T_NATIVE_DOUBLE, angles, coord_dims, error)
call h5dclose_f(dataset5_id, error) ! close dataset



!__________Read tind ________________________________
call h5dopen_f(file_id, dataset6, dataset6_id, error) ! open dataset
if(error == -1) then
   print*, 'Cannot open dataset radius in file:', file
   stop
end if


call h5dget_space_f(dataset6_id, dataspace6_id, error) 
call H5sget_simple_extent_dims_f(dataspace6_id, dims_out, coord_dims, error)

allocate(tind(coord_dims(2)))
call h5dread_f(dataset6_id, H5T_NATIVE_INTEGER, tind, coord_dims, error)

call h5dclose_f(dataset6_id, error) ! close dataset


!_____________________________________________________________________

call h5fclose_f(file_id, error) ! close file
call h5close_f(error) ! close inteface

allocate(param(size(param_r)))
param = dcmplx(param_r, param_i)




end subroutine read_geometry



subroutine T_write2file(Taa, Tab, Tba, Tbb, fname)
complex(dp), intent(in) :: Taa(:,:)
complex(dp), intent(in) :: Tab(:,:)
complex(dp), intent(in) :: Tba(:,:)
complex(dp), intent(in) :: Tbb(:,:)

CHARACTER(LEN=38), intent(in) :: fname
CHARACTER(LEN=38) :: filename

CHARACTER(LEN=5), PARAMETER :: dsetname1 = "Taa_r" ! Dataset name
CHARACTER(LEN=5), PARAMETER :: dsetname2 = "Taa_i" ! Dataset name

CHARACTER(LEN=5), PARAMETER :: dsetname3 = "Tab_r" ! Dataset name
CHARACTER(LEN=5), PARAMETER :: dsetname4 = "Tab_i" ! Dataset name

CHARACTER(LEN=5), PARAMETER :: dsetname5 = "Tba_r" ! Dataset name
CHARACTER(LEN=5), PARAMETER :: dsetname6 = "Tba_i" ! Dataset name

CHARACTER(LEN=5), PARAMETER :: dsetname7 = "Tbb_r" ! Dataset name
CHARACTER(LEN=5), PARAMETER :: dsetname8 = "Tbb_i" ! Dataset name


INTEGER(HID_T) :: file_id       ! File identifier
INTEGER(HID_T) :: dset_id1       ! Dataset identifier
INTEGER(HID_T) :: dset_id2       ! Dataset identifier
INTEGER(HID_T) :: dset_id3       ! Dataset identifier
INTEGER(HID_T) :: dset_id4       ! Dataset identifier
INTEGER(HID_T) :: dset_id5       ! Dataset identifier
INTEGER(HID_T) :: dset_id6       ! Dataset identifier
INTEGER(HID_T) :: dset_id7       ! Dataset identifier
INTEGER(HID_T) :: dset_id8       ! Dataset identifier

INTEGER(HID_T) :: dspace_id     ! Dataspace identifier


INTEGER(HSIZE_T), DIMENSION(2) :: dims  ! Dataset dimensions
INTEGER     ::    rank = 2                       ! Dataset rank

INTEGER     ::   error ! Error flag
!INTEGER     :: i, j

filename = fname
dims = (/size(Taa,1),size(Taa,2)/)
     
CALL h5open_f(error)

!
! Create a new file using default properties.
!
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

!
! Create the dataspace.
!
CALL h5screate_simple_f(rank, dims, dspace_id, error)

!
! Create and write dataset using default properties.
!
CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, real(Taa), dims, error)

CALL h5dclose_f(dset_id1, error)
!____

CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, imag(Taa), dims, error)

CALL h5dclose_f(dset_id2, error)

!******************************************************************

CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id3, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, real(Tab), dims, error)

CALL h5dclose_f(dset_id3, error)
!____

CALL h5dcreate_f(file_id, dsetname4, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id4, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id4, H5T_NATIVE_DOUBLE, imag(Tab), dims, error)

CALL h5dclose_f(dset_id4, error)

!******************************************************************

CALL h5dcreate_f(file_id, dsetname5, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id5, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id5, H5T_NATIVE_DOUBLE, real(Tba), dims, error)

CALL h5dclose_f(dset_id5, error)
!____

CALL h5dcreate_f(file_id, dsetname6, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id6, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id6, H5T_NATIVE_DOUBLE, imag(Tba), dims, error)

CALL h5dclose_f(dset_id6, error)

!******************************************************************


CALL h5dcreate_f(file_id, dsetname7, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id7, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id7, H5T_NATIVE_DOUBLE, real(Tbb), dims, error)

CALL h5dclose_f(dset_id7, error)
!____

CALL h5dcreate_f(file_id, dsetname8, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id8, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id8, H5T_NATIVE_DOUBLE, imag(Tbb), dims, error)

CALL h5dclose_f(dset_id8, error)

!******************************************************************
! Terminate access to the data space.
!
CALL h5sclose_f(dspace_id, error)

!
! Close the file.
!
CALL h5fclose_f(file_id, error)

!
! Close FORTRAN interface.
!
CALL h5close_f(error)

end subroutine T_write2file



subroutine read_T(file, Taa, Tab, Tba, Tbb)
!CHARACTER(LEN=28), intent(in) :: fname

CHARACTER(LEN=38) :: file  ! File name
CHARACTER(LEN=16), PARAMETER :: dataset1 = "Taa_r"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset2 = "Taa_i"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset3 = "Tab_r"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset4 = "Tab_i"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset5 = "Tba_r"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset6 = "Tba_i"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset7 = "Tbb_r"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset8 = "Tbb_i"     ! Dataset name

integer(HID_T) :: file_id 
integer(HID_T) :: dataset1_id, dataset2_id, dataset3_id, dataset4_id
integer(HID_T) :: dataset5_id, dataset6_id, dataset7_id, dataset8_id
integer(HID_T) :: dataspace_id

integer(HSIZE_T), dimension(2) :: dims_out, dims

integer :: error

real(dp), dimension(:,:), allocatable :: Taa_r, Taa_i
real(dp), dimension(:,:), allocatable :: Tab_r, Tab_i
real(dp), dimension(:,:), allocatable :: Tba_r, Tba_i
real(dp), dimension(:,:), allocatable :: Tbb_r, Tbb_i
complex(dp), dimension(:,:), allocatable :: Taa, Tab, Tba, Tbb 


call h5open_f(error) ! initialize interface
call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error) ! open file

!__________ _____________________________
call h5dopen_f(file_id, dataset1, dataset1_id, error) ! open dataset
call h5dget_space_f(dataset1_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Taa_r(dims(1), dims(2)))
call h5dread_f(dataset1_id, H5T_NATIVE_DOUBLE, Taa_r,dims, error)
call h5dclose_f(dataset1_id, error) ! close dataset
!__________________________________________
call h5dopen_f(file_id, dataset2, dataset2_id, error) ! open dataset
call h5dget_space_f(dataset2_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Taa_i(dims(1), dims(2)))
call h5dread_f(dataset2_id, H5T_NATIVE_DOUBLE, Taa_i,dims, error)
call h5dclose_f(dataset2_id, error) ! close dataset
!__________________________________________
allocate(Taa(size(Taa_r,1),size(Taa_r,1)))
Taa = dcmplx(Taa_r, Taa_i)


!__________ _____________________________
call h5dopen_f(file_id, dataset3, dataset3_id, error) ! open dataset
call h5dget_space_f(dataset3_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tab_r(dims(1), dims(2)))
call h5dread_f(dataset3_id, H5T_NATIVE_DOUBLE, Tab_r,dims, error)
call h5dclose_f(dataset3_id, error) ! close dataset
!__________________________________________
call h5dopen_f(file_id, dataset4, dataset4_id, error) ! open dataset
call h5dget_space_f(dataset4_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tab_i(dims(1), dims(2)))
call h5dread_f(dataset4_id, H5T_NATIVE_DOUBLE, Tab_i,dims, error)
call h5dclose_f(dataset4_id, error) ! close dataset
!__________________________________________
allocate(Tab(size(Tab_r,1),size(Tab_r,1)))
Tab = dcmplx(Tab_r, Tab_i)


!__________ _____________________________
call h5dopen_f(file_id, dataset5, dataset5_id, error) ! open dataset
call h5dget_space_f(dataset5_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tba_r(dims(1), dims(2)))
call h5dread_f(dataset5_id, H5T_NATIVE_DOUBLE, Tba_r,dims, error)
call h5dclose_f(dataset5_id, error) ! close dataset
!__________________________________________
call h5dopen_f(file_id, dataset6, dataset6_id, error) ! open dataset
call h5dget_space_f(dataset6_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tba_i(dims(1), dims(2)))
call h5dread_f(dataset6_id, H5T_NATIVE_DOUBLE, Tba_i,dims, error)
call h5dclose_f(dataset6_id, error) ! close dataset
!__________________________________________
allocate(Tba(size(Tba_r,1),size(Tba_r,1)))
Tba = dcmplx(Tba_r, Tba_i)




!__________ _____________________________
call h5dopen_f(file_id, dataset7, dataset7_id, error) ! open dataset
call h5dget_space_f(dataset7_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tbb_r(dims(1), dims(2)))
call h5dread_f(dataset7_id, H5T_NATIVE_DOUBLE, Tbb_r,dims, error)
call h5dclose_f(dataset7_id, error) ! close dataset
!__________________________________________
call h5dopen_f(file_id, dataset8, dataset8_id, error) ! open dataset
call h5dget_space_f(dataset8_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tbb_i(dims(1), dims(2)))
call h5dread_f(dataset8_id, H5T_NATIVE_DOUBLE, Tbb_i,dims, error)
call h5dclose_f(dataset8_id, error) ! close dataset
!__________________________________________
allocate(Tbb(size(Tbb_r,1),size(Tbb_r,1)))
Tbb = dcmplx(Tbb_r, Tbb_i)



!_______________________________________________________________

call h5fclose_f(file_id, error) ! close file
call h5close_f(error) ! close inteface


end subroutine read_T





subroutine read_T2(file, Taa, Tab, Tba, Tbb)
!CHARACTER(LEN=28), intent(in) :: fname

CHARACTER(LEN=38) :: file  ! File name
CHARACTER(LEN=16), PARAMETER :: dataset1 = "Taa_r"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset2 = "Taa_i"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset3 = "Tab_r"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset4 = "Tab_i"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset5 = "Tba_r"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset6 = "Tba_i"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset7 = "Tbb_r"     ! Dataset name
CHARACTER(LEN=16), PARAMETER :: dataset8 = "Tbb_i"     ! Dataset name
!CHARACTER(LEN=16), PARAMETER :: dataset9 = "Cexts"     ! Dataset name

integer(HID_T) :: file_id 
integer(HID_T) :: dataset1_id, dataset2_id, dataset3_id, dataset4_id
integer(HID_T) :: dataset5_id, dataset6_id, dataset7_id, dataset8_id
!integer(HID_T) :: dataset9_id

integer(HID_T) :: dataspace_id

integer(HSIZE_T), dimension(3) :: dims_out, dims
integer(HSIZE_T), dimension(2) :: dims_out2, dims2

integer :: error

double precision, dimension(:,:,:), allocatable :: Taa_r, Taa_i
double precision, dimension(:,:,:), allocatable :: Tab_r, Tab_i
double precision, dimension(:,:,:), allocatable :: Tba_r, Tba_i
double precision, dimension(:,:,:), allocatable :: Tbb_r, Tbb_i
complex(dp), dimension(:,:,:), allocatable :: Taa, Tab, Tba, Tbb 
!double precision, allocatable, intent(out) :: Cexts(:,:)

call h5open_f(error) ! initialize interface
call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error) ! open file

!__________ _____________________________
call h5dopen_f(file_id, dataset1, dataset1_id, error) ! open dataset
call h5dget_space_f(dataset1_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Taa_r(dims(1), dims(2), dims(3)))
call h5dread_f(dataset1_id, H5T_NATIVE_DOUBLE, Taa_r,dims, error)
call h5dclose_f(dataset1_id, error) ! close dataset
!__________________________________________
call h5dopen_f(file_id, dataset2, dataset2_id, error) ! open dataset
call h5dget_space_f(dataset2_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Taa_i(dims(1), dims(2), dims(3)))
call h5dread_f(dataset2_id, H5T_NATIVE_DOUBLE, Taa_i,dims, error)
call h5dclose_f(dataset2_id, error) ! close dataset
!__________________________________________
allocate(Taa(size(Taa_r,1),size(Taa_r,2),size(Taa_r,3)))
Taa = dcmplx(Taa_r, Taa_i)


!__________ _____________________________
call h5dopen_f(file_id, dataset3, dataset3_id, error) ! open dataset
call h5dget_space_f(dataset3_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tab_r(dims(1), dims(2), dims(3)))
call h5dread_f(dataset3_id, H5T_NATIVE_DOUBLE, Tab_r,dims, error)
call h5dclose_f(dataset3_id, error) ! close dataset
!__________________________________________
call h5dopen_f(file_id, dataset4, dataset4_id, error) ! open dataset
call h5dget_space_f(dataset4_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tab_i(dims(1), dims(2), dims(3)))
call h5dread_f(dataset4_id, H5T_NATIVE_DOUBLE, Tab_i,dims, error)
call h5dclose_f(dataset4_id, error) ! close dataset
!__________________________________________
allocate(Tab(size(Tab_r,1),size(Tab_r,2),size(Tab_r,3)))
Tab = dcmplx(Tab_r, Tab_i)


!__________ _____________________________
call h5dopen_f(file_id, dataset5, dataset5_id, error) ! open dataset
call h5dget_space_f(dataset5_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tba_r(dims(1), dims(2), dims(3)))
call h5dread_f(dataset5_id, H5T_NATIVE_DOUBLE, Tba_r,dims, error)
call h5dclose_f(dataset5_id, error) ! close dataset
!__________________________________________
call h5dopen_f(file_id, dataset6, dataset6_id, error) ! open dataset
call h5dget_space_f(dataset6_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tba_i(dims(1), dims(2), dims(3)))
call h5dread_f(dataset6_id, H5T_NATIVE_DOUBLE, Tba_i,dims, error)
call h5dclose_f(dataset6_id, error) ! close dataset
!__________________________________________
allocate(Tba(size(Tba_r,1),size(Tba_r,2),size(Tba_r,3)))
Tba = dcmplx(Tba_r, Tba_i)




!__________ _____________________________
call h5dopen_f(file_id, dataset7, dataset7_id, error) ! open dataset
call h5dget_space_f(dataset7_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tbb_r(dims(1), dims(2), dims(3)))
call h5dread_f(dataset7_id, H5T_NATIVE_DOUBLE, Tbb_r,dims, error)
call h5dclose_f(dataset7_id, error) ! close dataset
!__________________________________________
call h5dopen_f(file_id, dataset8, dataset8_id, error) ! open dataset
call h5dget_space_f(dataset8_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, dims_out, dims, error)
allocate(Tbb_i(dims(1), dims(2), dims(3)))
call h5dread_f(dataset8_id, H5T_NATIVE_DOUBLE, Tbb_i,dims, error)
call h5dclose_f(dataset8_id, error) ! close dataset
!__________________________________________
allocate(Tbb(size(Tbb_r,1),size(Tbb_r,2),size(Tbb_r,3)))
Tbb = dcmplx(Tbb_r, Tbb_i)



!__________________________________________
!call h5dopen_f(file_id, dataset9, dataset9_id, error) ! open dataset
!call h5dget_space_f(dataset9_id, dataspace_id, error) 
!call H5sget_simple_extent_dims_f(dataspace_id, dims_out2, dims2, error)
!allocate(Cexts(dims2(1),dims2(2)))
!call h5dread_f(dataset9_id, H5T_NATIVE_DOUBLE, Cexts,dims2, error)
!call h5dclose_f(dataset9_id, error) ! close dataset




!_______________________________________________________________

call h5fclose_f(file_id, error) ! close file
call h5close_f(error) ! close inteface




end subroutine read_T2





subroutine read_mat(A, file)

CHARACTER(LEN=38) :: file ! File name
CHARACTER(LEN=16), PARAMETER :: dataset1 = "A"     ! Dataset name


integer(HID_T) :: file_id 
integer(HID_T) :: dataset1_id
integer(HID_T) :: dataspace_id
integer(HSIZE_T), dimension(2) :: A_dims, A_out

integer :: error

real(dp), dimension(:,:), allocatable :: A

call h5open_f(error) ! initialize interface
call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error) ! open file


call h5dopen_f(file_id, dataset1, dataset1_id, error) ! open dataset

call h5dget_space_f(dataset1_id, dataspace_id, error) 
call H5sget_simple_extent_dims_f(dataspace_id, A_out, A_dims, error)

allocate(A(A_dims(1), A_dims(2)))
call h5dread_f(dataset1_id, H5T_NATIVE_DOUBLE, A, A_dims, error)

call h5dclose_f(dataset1_id, error) ! close dataset

call h5fclose_f(file_id, error) ! close file
call h5close_f(error) ! close inteface


end subroutine read_mat

subroutine write2file_mueller(A, crs, fname)
real(dp), intent(in) :: A(:,:)
real(dp), intent(in) :: crs(4)
CHARACTER(LEN=38), intent(in) :: fname
CHARACTER(LEN=38) :: filename 

CHARACTER(LEN=7), PARAMETER :: dsetname1 = "mueller" ! Dataset name
CHARACTER(LEN=14), PARAMETER :: dsetname2 = "cross_sections" ! Dataset name

INTEGER(HID_T) :: file_id       ! File identifier
INTEGER(HID_T) :: dset_id1, dset_id2      ! Dataset identifier
INTEGER(HID_T) :: dspace_id , dspace_id2    ! Dataspace identifier


INTEGER(HSIZE_T), DIMENSION(2) :: dims! Dataset dimensions
INTEGER(HSIZE_T), DIMENSION(1) :: dims2
INTEGER     ::    rank = 2                       ! Dataset rank
INTEGER     ::    rank2 = 1  

INTEGER     ::   error ! Error flag
!INTEGER     :: i, j

dims = (/size(A,1),size(A,2)/)
dims2 = 4
filename = fname  
!print*,crs
!crs=[Cext,Csca,Cabs,Csca_int]

CALL h5open_f(error)

!
! Create a new file using default properties.
!
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

!
! Create the dataspace.
!
CALL h5screate_simple_f(rank, dims, dspace_id, error)

!
! Create and write dataset using default properties.
!
CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id1, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, A, dims, error)


!
! End access to the dataset and release resources used by it.
!
CALL h5dclose_f(dset_id1, error)

!
! Terminate access to the data space.
!
CALL h5sclose_f(dspace_id, error)

!****************************************************

!
! Create the dataspace.
!
CALL h5screate_simple_f(rank2, dims2, dspace_id2, error)

!
! Create and write dataset using default properties.
!
CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id2, &
                       dset_id2, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

CALL h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, crs, dims2, error)

!
! End access to the dataset and release resources used by it.
!
CALL h5dclose_f(dset_id2, error)

!
! Terminate access to the data space.
!
CALL h5sclose_f(dspace_id2, error)


!
! Close the file.
!
CALL h5fclose_f(file_id, error)

!
! Close FORTRAN interface.
!
CALL h5close_f(error)

end subroutine write2file_mueller
!_---______________________________________________________________



end module 
