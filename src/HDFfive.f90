!******************************************************************::
!
! The module contains functions for reading and writing data in HDF5
! format
!
!**********************************************************************
MODULE HDFfive
   use HDF5
   implicit none
   private
   !
   public ::  h5_crtgrp, h5_wrt2file, h5_wrtvec2file, h5_rd_vec, h5_rd_file   

contains

   subroutine h5_crtgrp(filename_, main_grpname, subgrpsname)
      !============================================================
      ! This subroutine creates subgroups in an existing group.
      !============================================================
      !USE HDF5 ! This module contains all necessary modules

      IMPLICIT NONE
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables
      CHARACTER(LEN=64), intent(in) :: filename_  != "groupf.h5" ! File name
      CHARACTER(LEN=32), intent(in) :: main_grpname != "MyGroup"  ! Group name
      CHARACTER(LEN=32), intent(in), optional :: subgrpsname(:)
      ! Local variables
      INTEGER  :: num_sub = 0
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: main_group_id
      INTEGER(HID_T), allocatable :: group_id(:)
      INTEGER     ::   error ! Error flag
      INTEGER     :: idum
      !-----------------------------------
      if (present(subgrpsname)) then
         num_sub = size(subgrpsname)
         allocate (group_id(size(subgrpsname)))
      end if

      ! Initialize FORTRAN interface.

      CALL h5open_f(error)
      ! Open an existing file.
      !
      CALL h5fcreate_f(filename_, H5F_ACC_TRUNC_F, file_id, error)
      ! Open an existing group in the specified file.
      !write(*,*)'error_1_sub=', error
      CALL h5gcreate_f(file_id, main_grpname, main_group_id, error) !(file_id, groupname, group_id, error)
      !  write(*,*)'error_2_sub=', error
      !creating subgroups inside main_group
!

      do idum = 1, num_sub
         ! Create sub_group in main group.
         CALL h5gcreate_f(file_id, subgrpsname(idum), group_id(idum), error)
         !  write(*,*)'error_3_sub=', error
         ! Close the sub_group.
         CALL h5gclose_f(group_id(idum), error)
         !  write(*,*)'error_4_sub=', error
         !write(*,*)'group_id(idum)',group_id(idum)
      end do

! Close the main group.
      !
      CALL h5gclose_f(main_group_id, error)  !(group_id, error)

      ! Terminate access to the file.
      !
      CALL h5fclose_f(file_id, error)

      !
      ! Close FORTRAN interface.
      !
      CALL h5close_f(error)
      if (allocated(group_id)) deallocate(group_id)
   end subroutine h5_crtgrp

!--------------------------------------------------------------------
   SUBROUTINE h5_wrtvec2file(filename_, groupname, dsetname, dset_data, attribute)  !write_real data in a group
      !============================================================
      ! This subroutine writes data in a dataset in an existing group.
      !============================================================

      !USE HDF5 ! This module contains all necessary modules

      IMPLICIT NONE

      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables

      CHARACTER(LEN=64), intent(in) :: filename_ != "groupsf.h5" ! File name
      CHARACTER(LEN=32), intent(in) :: groupname != "MyGroup/Group_A" ! Group name
      CHARACTER(LEN=32), intent(in) :: dsetname  ! = "MyGroup/dset1"  ! Dataset name
      CHARACTER(LEN=256), OPTIONAL, intent(in) :: attribute
      real(8), intent(in)  :: dset_data(:)

      ! Local variables
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: group_id      ! Group identifier
      INTEGER(HID_T) :: dataset_id    ! Dataset identifier
      INTEGER(HID_T) :: dataspace_id ! Data space identifier
      INTEGER ::   error ! Error flag
      INTEGER(HSIZE_T), DIMENSION(1) :: dims  ! (/2,3/)!Datasets dimensions
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      INTEGER     ::   rank = 1 ! Datasets rank

      dims = (/size(dset_data)/)

      !
      ! Initialize FORTRAN interface.
      !
      CALL h5open_f(error)
      !
      ! Open an existing file.
      CALL h5fopen_f(filename_, H5F_ACC_RDWR_F, file_id, error)

      ! Open an existing group in the specified file.

      CALL h5gopen_f(file_id, groupname, group_id, error)  !groupname should be complete!

      !Create the data space for the second dataset.
      !
      CALL h5screate_simple_f(rank, dims, dataspace_id, error)

      ! Create the dataset in group "Group_A" with default properties.
      !
      CALL h5dcreate_f(group_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
                       dataset_id, error)

      ! Write the second dataset.
      !
      data_dims(1) = dims(1)
      CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset_data, dims, error)

      if (present(attribute)) then
         CALL h5_wrt_attr(attribute, dataset_id)
      end if

      ! Close the dataspace for the second dataset.
      !
      CALL h5sclose_f(dataspace_id, error)

      ! Close the second dataset.
      !
      CALL h5dclose_f(dataset_id, error)

      ! Close the group.
      !
      CALL h5gclose_f(group_id, error)

      ! Close the file.
      !
      CALL h5fclose_f(file_id, error)

      !
      ! Close FORTRAN interface.
      !
      CALL h5close_f(error)

   END SUBROUTINE h5_wrtvec2file

!--------------------------------------------------------------------
   SUBROUTINE h5_wrt2file(filename_, groupname, dsetname, dset_data, attribute)  !write_real data in a group
      !============================================================
      ! This subroutine writes data in a dataset in an existing group.
      !============================================================

      !USE HDF5 ! This module contains all necessary modules

      IMPLICIT NONE

      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables

      CHARACTER(LEN=64), intent(in) :: filename_ != "groupsf.h5" ! File name
      CHARACTER(LEN=32), intent(in) :: groupname != "MyGroup/Group_A" ! Group name
      CHARACTER(LEN=32), intent(in) :: dsetname  ! = "MyGroup/dset1"  ! Dataset name
      REAL(8), intent(in)  :: dset_data(:, :)
      CHARACTER(LEN=256), OPTIONAL, intent(in) :: attribute
      ! Local variables
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: group_id      ! Group identifier
      INTEGER(HID_T) :: dataset_id    ! Dataset identifier
      INTEGER(HID_T) :: dataspace_id  ! Data space identifier
      INTEGER ::   error ! Error flag
      INTEGER(HSIZE_T), DIMENSION(2) :: dims  ! (/2,3/)!Datasets dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
      INTEGER     ::   rank = 2 ! Datasets rank
      dims = (/size(dset_data, 1), size(dset_data, 2)/)
      !
      !
      ! Initialize FORTRAN interface.
      !
      CALL h5open_f(error)
      !
      ! Open an existing file.
      CALL h5fopen_f(filename_, H5F_ACC_RDWR_F, file_id, error)
      ! Open an existing group in the specified file.

      CALL h5gopen_f(file_id, groupname, group_id, error)  !groupname should be complete!

      !Create the data space for the second dataset.
      !
      CALL h5screate_simple_f(rank, dims, dataspace_id, error)                     

      ! Create the second dataset in group "Group_A" with default properties.
      !
      CALL h5dcreate_f(group_id, dsetname, H5T_NATIVE_DOUBLE, dataspace_id, &
                       dataset_id, error)

      ! Write the second dataset.
      !
      data_dims(1) = dims(1)
      data_dims(2) = dims(2)
      CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dset_data, data_dims, error)
      if (present(attribute)) then
         CALL h5_wrt_attr(attribute, dataset_id)
      end if
      ! Close the dataspace for the second dataset.
      !
      CALL h5sclose_f(dataspace_id, error)

      ! Close the second dataset.
      !
      CALL h5dclose_f(dataset_id, error)

      ! Close the group.
      !
      CALL h5gclose_f(group_id, error)

      ! Close the file.
      !
      CALL h5fclose_f(file_id, error)

      !
      ! Close FORTRAN interface.
      !
      CALL h5close_f(error)
	!if (allocated(group_id)) deallocate(group_id)
   END SUBROUTINE h5_wrt2file
!-------------------------------------------------------------------

   SUBROUTINE h5_wrt_attr(attribute, dataset_id)

      IMPLICIT NONE

      CHARACTER(LEN=256), intent(in) :: attribute
      INTEGER(HID_T), intent(in) :: dataset_id
      !Local variables
      INTEGER(HID_T) :: dspace_id  ! Data space identifier
      INTEGER(HID_T) :: attr_id       ! attribute identifier
      INTEGER(HID_T) :: type_id     ! Attribute datatype identifier
      INTEGER        :: error ! Error flag
      INTEGER(HSIZE_T), DIMENSION(1) :: dims  !
      CHARACTER(LEN=12) :: attr = 'Description:'
      INTEGER(8)    ::   dp = 256
      dims(1) = 1
      !dp=11
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, error)
      !print *, 'h5tcopy_f returns', type_id
      CALL h5tset_size_f(type_id, dp, error)
      !print *, 'h5tset_size_f returns', error

      ! Create the dataspace.
      !
      CALL h5screate_f(H5S_SCALAR_F, dspace_id, error)
      !print *, 'h5screate_f returns', dspace_id

      !
      ! Create the dataset with default properties.
      !
      CALL h5acreate_f(dataset_id, attr, type_id, dspace_id, &  !dataset_id
                       attr_id, error)

      CALL h5awrite_f(attr_id, type_id, attribute, dims, error)

      CALL h5tclose_f(type_id, error)

      CALL h5aclose_f(attr_id, error)

      CALL h5sclose_f(dspace_id, error)

   END SUBROUTINE h5_wrt_attr

!------------------------------------------------------------------
  SUBROUTINE h5_rd_vec(filename_, groupname, dsetname, dset_data) !, attribute)
      !============================================================
      ! This subroutine reads data in a dataset in an existing group.
      !============================================================

      !USE HDF5 ! This module contains all necessary modules

      IMPLICIT NONE

      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables

    CHARACTER(*), intent(in) :: filename_ ! File name
    !CHARACTER(LEN=4), PARAMETER :: dsetname = "FT"     ! Dataset name
    CHARACTER(*), intent(in) :: groupname != "MyGroup/Group_A" ! Group name
    CHARACTER(*), intent(in) :: dsetname  ! = "MyGroup/dset1"  ! Dataset name
   ! CHARACTER(LEN=256), OPTIONAL, intent(in) :: attribute !Ati: I think this is not needed
    !REAL(8), intent(out)  :: dset_data(:,:)      ! output data
      
      
    INTEGER(HID_T) :: file_id       ! File identifier
    INTEGER(HID_T) :: dataset_id       ! Dataset identifier
    INTEGER(HID_T) :: space_id       ! Dataspace identifier
    INTEGER(HID_T) :: dtype_id       ! Dataspace identifier
    INTEGER(HID_T) :: group_id      ! Group identifier
    INTEGER     ::   error ! Error flag
    INTEGER     ::  i, j, cols, rows
    INTEGER(HSIZE_T) :: npoints
    REAL(KIND = 8), intent(out), DIMENSION(:), ALLOCATABLE :: dset_data
    INTEGER(HSIZE_T), DIMENSION(1) :: dims, maxdims                 
    INTEGER :: rank   
    
  
    !print *, 'Starting HDF5 Fortran Read'

   ! Initialize FORTRAN interface.

   CALL h5open_f(error)
  

   ! Open an existing file.

   CALL h5fopen_f (filename_, H5F_ACC_RDWR_F, file_id, error)   
    
   CALL h5gopen_f(file_id, groupname, group_id, error)  !groupname should be complete!
   ! Open an existing dataset.
    
   CALL h5dopen_f(group_id, dsetname, dataset_id, error)   ! CALL h5dopen_f(file_id, dsetname, dataset_id, error) original
   

   !Get dataspace ID
   CALL h5dget_space_f(dataset_id, space_id,error)
   

   !Get dataspace dims
   CALL h5sget_simple_extent_ndims_f (space_id, rank, error)
  ! print *, 'rank'
 !  print *, rank
   
   CALL h5sget_simple_extent_dims_f(space_id,dims, maxdims, error)
  ! print *, dims
  ! print *, maxdims
   
   if (rank == 0) then
   	dims(1) = 1
   end if
   if (ALLOCATED(dset_data)) DEALLOCATE(dset_data)
   ALLOCATE(dset_data(dims(1)))
   !Get data
    CALL h5dread_f(dataset_id,  H5T_NATIVE_DOUBLE, dset_data, dims,  error)  ! CALL h5dread_f(dataset_id, H5T_NATIVE_INTEGER, dset_data, data_dims, error) 
    CALL h5dclose_f(dataset_id, error)
    CALL h5sclose_f(space_id, error)
    CALL h5gclose_f(group_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
    
	!if (allocated(group_id)) deallocate(group_id)
   END SUBROUTINE h5_rd_vec

!--------------------------------------------------------
SUBROUTINE h5_rd_file(filename_, groupname, dsetname, dset_data) !, attribute)
      !============================================================
      ! This subroutine reads data in a dataset in an existing group.
      !============================================================

      !USE HDF5 ! This module contains all necessary modules
      USE ISO_C_BINDING
      IMPLICIT NONE

      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables

    CHARACTER(*), intent(in) :: filename_ ! File name
    !CHARACTER(LEN=4), PARAMETER :: dsetname = "FT"     ! Dataset name
    CHARACTER(*), intent(in) :: groupname != "MyGroup/Group_A" ! Group name
    CHARACTER(*), intent(in) :: dsetname  ! = "MyGroup/dset1"  ! Dataset name
   ! CHARACTER(LEN=256), OPTIONAL, intent(in) :: attribute !I think this is not needed
   
    INTEGER, PARAMETER :: r_k8 = KIND(0.0d0)
    COMPLEX(KIND = r_k8), intent(out), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: dset_data   ! output data
      
    INTEGER(HID_T) :: file_id       ! File identifier
    INTEGER(HID_T) :: dataset_id       ! Dataset identifier
    INTEGER(HID_T) :: space_id       ! Dataspace identifier
    INTEGER(HID_T) :: dtype_id       ! Dataspace identifier
    INTEGER(HID_T) :: group_id      ! Group identifier
    INTEGER(HID_T) :: memtype      ! Group identifier
    INTEGER(HID_T) :: sample_type_id
    INTEGER     ::   error ! Error flag
    INTEGER     ::  i, j, cols, rows, shell, rank
    
    INTEGER(HSIZE_T), DIMENSION(3) :: dims
    INTEGER(HSIZE_T), DIMENSION(3) :: maxdims                  
    TYPE(C_PTR) :: f_ptr
    INTEGER(8) :: real_size, real_complex_size
    real_size = storage_size(1_r_k8, r_k8) / 8
    real_complex_size = real_size * 2_8  ! a complex is (real,real)
    
    
    
    if (ALLOCATED(dset_data)) DEALLOCATE(dset_data)
    !print *, 'Starting HDF5 Fortran Read'

   ! Initialize FORTRAN interface.
   
   CALL h5open_f(error)


   ! Open an existing file.

   CALL h5fopen_f (filename_, H5F_ACC_RDWR_F, file_id, error)

   CALL h5gopen_f(file_id, groupname, group_id, error)  !groupname should be complete!
   ! Open an existing dataset.

   CALL h5dopen_f(group_id, dsetname, dataset_id, error)   ! CALL h5dopen_f(file_id, dsetname, dataset_id, error) original

   
   !Get dataspace ID
   CALL h5dget_space_f(dataset_id, space_id,error)


   !Get dataspace dims
   CALL h5sget_simple_extent_ndims_f (space_id, rank, error)
  ! print *, 'rank'
  ! print *, rank
      
   CALL h5sget_simple_extent_dims_f(space_id,dims, maxdims, error)

   cols = dims(1)
   rows = dims(2)
   shell= dims(3)
   if (rank == 0)then
       cols = 1
       rows = 1
       shell = 1 
    end if
  
   ALLOCATE(dset_data(cols, rows, shell))
   CALL H5Tcreate_f(H5T_COMPOUND_F, 16_8, sample_type_id,error)   !Creates a new datatype, 16_8: for double complex

   CALL H5Tinsert_f( sample_type_id, "r", &
      	 0_8, h5kind_to_type(r_k8,H5_REAL_KIND), error)
   CALL H5Tinsert_f( sample_type_id, "i", &
       real_size, h5kind_to_type(r_k8,H5_REAL_KIND), error)
   	!Get data
        f_ptr = C_LOC(dset_data(1,1,1))
  	! CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, dset_data, coord_dims, error)  ! CALL h5dread_f(dataset_id, H5T_NATIVE_INTEGER, dset_data, data_dims, error)
   CALL h5dread_f(dataset_id, sample_type_id, f_ptr, error) 

     !print *, dset_data


    !CALL h5sclose_f(space_id, error)

      ! Close the second dataset.
      !
    CALL h5dclose_f(dataset_id, error)
    CALL H5Tclose_f(sample_type_id, error)
      ! Close the group.
      !
    CALL h5gclose_f(group_id, error)

      ! Close the file.
      !
      CALL h5fclose_f(file_id, error)

      !
      ! Close FORTRAN interface.
      !
      CALL h5close_f(error)
	!if (allocated(group_id)) deallocate(group_id)
   END SUBROUTINE h5_rd_file


END MODULE
