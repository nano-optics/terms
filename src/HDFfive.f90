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
   public ::  h5_crtgrp, h5_wrt2file, h5_wrtvec2file

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
      if (allocated(group_id)) deallocate (group_id)
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

END MODULE
