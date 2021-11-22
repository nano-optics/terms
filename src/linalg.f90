module linalg
   !
   !==================================================================
   ! This module contains "wrappers" to drive LAPACK's square-matrix
   ! inversion routines and linear solvers.
   !==================================================================
   !
   implicit none
   !
   private
   !
   public :: invSqrMat, solLinSys
   !
contains
   !
   !******************************************************************
   !
   subroutine invSqrMat(A, trans_, verb_)
      !
      !================================================================
      ! Wrapper function for inverting a complex-valued square matrix
      ! A(n,n), using the ZGETRF and ZGETRI routines in LAPACK. Note
      ! that A is overwritten by inv(A) on output.
      !================================================================
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables:
      complex(8), intent(inout) :: A(:, :)
      integer, intent(in), optional :: verb_
      logical, intent(in), optional :: trans_
      ! Local variables:
      complex(8), allocatable :: work(:)
      real(8) :: rcond, anorm, ZLANGE
      real(8), allocatable :: rwork(:)
      logical                 :: trans
      integer                 :: n, lwork, info, ipiv(size(A, 1)), verb !tic,toc,tps,
      !real                    :: t0,t1
      integer, parameter :: verb0 = 1
      character(*), parameter :: myName = 'invSqrMat'
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      ! call cpu_time(t0)
      !call system_clock(tic)
      !
      verb = 0
      if (present(verb_)) verb = verb_
      !
      trans = .false.
      if (present(trans_)) trans = trans_
      !
      n = size(A, 1)
      if (n /= size(A, 2)) then
         write (*, '(A,A)') myName, '> ERROR: Matrix not square! Stopping..'
         STOP
      end if
      !
      if (trans) A = transpose(A)
      !
      if (verb > verb0) then
         allocate (work(2*n), rwork(2*n))
         anorm = ZLANGE('1', n, n, A, n, rwork(1)) ! the 1-norm of matrix
         write (*, '(A,A,es12.4E3)') &
            myname, '> ZLANGE: Supplied matrix 1-norm= ', anorm
         call ZGECON('1', n, A, n, anorm, rcond, work, rwork, info)
         if (info /= 0) then
            write (*, '(A,A,i3,A)') &
               myname, '> ERROR: ZGECON returned non-zero INFO= ', info
            STOP
         else
            write (*, '(A,A,es12.4E3)') &
               myname, '> ZGECON returned RCOND= ', rcond
         end if
         deallocate (work, rwork)
      end if
      !
      allocate (work(1))
      call ZGETRI(n, A, n, ipiv, work, -1, info) !  Get optimal LWORK
      if (info == 0) then
         lwork = int(work(1))
      else
         write (*, '(A,A,i3,A)') myname, &
            '> ERROR: ZGETRI returned non-zero INFO= ', info, ' when querying'
         STOP
      end if
      deallocate (work)
      allocate (work(lwork))
      !
      call ZGETRF(n, n, A, n, ipiv, info)
      !
      if (info /= 0) then
         write (*, '(A,A,i3)') myname, &
            '> ERROR: ZGETRF returned non-zero INFO= ', info
         STOP
      end if
      !
      call ZGETRI(n, A, n, ipiv, work, lwork, info)
      deallocate (work)
      if (info /= 0) then
         write (*, '(A,A,i3,A)') myname, &
            '> ERROR: ZGETRI returned non-zero INFO= ', info, ' when inverting'
         STOP
      end if
      !
      if (trans) A = transpose(A)
      !
      !call cpu_time(t1)
      !call system_clock(toc,count_rate=tps)
      !if(verb > verb0) write(*,'(A,A,2(1x,es9.3E2))') &
      !     myName,'> Calculation time (CPU & real in s): ', t1-t0, real(toc-tic)/real(tps)
      !
   end subroutine invSqrMat
   !
   !******************************************************************
   !
   subroutine solLinSys(A, X, isol_, verb_)
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables:
      complex(8), intent(inout) :: A(:, :), X(:, :)
      integer, intent(in), optional :: isol_, verb_
      ! Local variables:
      character(*), parameter :: myName = 'solLinSys'
      integer, parameter :: verb0 = 1
      integer :: isol, verb, n, info !tic,tps, toc
      !real :: t0,t1
      real(8) :: rcond, anorm, ZLANGE
      real(8), allocatable :: rwork(:)
      complex(8), allocatable  :: work(:)
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      !call cpu_time(t0)
      !call system_clock(tic)
      !
      isol = 0
      if (present(isol_)) isol = isol_
      !
      verb = 0
      if (present(verb_)) verb = verb_
      !
      n = size(A, 1)
      if (n /= size(A, 2)) then
         write (*, '(A,A)') myName, '> ERROR: size(A,1) != size(A,2)'
         STOP
      end if
      !
      allocate (rwork(2*n), work(2*n))
      anorm = ZLANGE('1', n, n, A, n, rwork(1)) ! the 1-norm of matrix
      if (verb > verb0) write (*, '(A,A,es12.4E3)') &
         myname, '> ZLANGE: Supplied matrix 1-norm= ', anorm
      call ZGECON('1', n, A, n, anorm, rcond, work, rwork, info)
      if (info /= 0) then
         write (*, '(A,A,i3)') &
            myname, '> ERROR: ZGECON returned non-zero INFO= ', info
         STOP
      elseif (verb > verb0) then
         write (*, '(A,A,es12.4E3)') &
            myname, '> ZGECON returned RCOND= ', rcond
      end if
      deallocate (rwork, work)
      !
      if (isol == 0) then
         if (verb > verb0) write (*, '(A,A,i2,A)') myName, &
            '> isol= ', isol, ' triggered wrapper for LAPACK''s ZGESV'
         call solLinSysV(A, X, verb_=verb)
      elseif (isol == 1) then
         if (verb > verb0) write (*, '(A,A,i2,A)') myName, &
            '> isol= ', isol, ' triggered wrapper for LAPACK''s ZGESVX'
         call solLinSysVX(A, X, verb_=verb)
         ! elseif(isol == 2) then
         !    write(*,'(A,A)') myName,'> ERROR: isol= 2 currently usupported'
         !    STOP
         !    if(verb > 0) write(*,'(A,A,i2,A)') myName, &
         !         '> isol= ', isol,' triggered wrapper for LAPACK''s ZGESVXX'
         !    call solLinSysVXX(A,X, verb_ = verb)
         ! elseif(isol == 3) then
         !    if(verb > 0) write(*,'(A,A,i2,A)') myName, &
         !         '> isol= ', isol,' triggered wrapper for LAPACK''s ZGELS'
         !    call solLinSysLS(A,X, verb_ = verb)
         ! elseif(isol == 4) then
         !    if(verb > 0) write(*,'(A,A,i2,A)') myName, &
         !         '> isol= ', isol,' triggered wrapper for LAPACK''s ZGELSS'
         !    call solLinSysLSS(A,X, verb_ = verb)
         ! elseif(isol == 5) then
         !    if(verb > 0) write(*,'(A,A,i2,A)') myName, &
         !         '> isol= ', isol,' triggered wrapper for LAPACK''s ZGELSY'
         !    call solLinSysLSY(A,X, verb_ = verb)
         ! elseif(isol == 6) then
         !    if(verb > 0) write(*,'(A,A,i2,A)') myName, &
         !         '> isol= ', isol,' triggered wrapper for LAPACK''s ZGETSLS'
         !    call solLinSysTSLS(A,X, verb_ = verb)
      else
         write (*, '(A,A,i2)') myName, '> ERROR: Unrecognised isol= ', isol
         STOP
      end if
      !
      !call cpu_time(t1)
      !call system_clock(toc,count_rate=tps)
      ! if(verb > verb0) write(*,'(A,A,2(1x,es9.3E2))') &
      !    myName,'> Calculation time (CPU & real in s): ', t1-t0, real(toc-tic)/real(tps)
      !
   end subroutine solLinSys
   !
   !******************************************************************
   !
   subroutine solLinSysV(A, X, verb_)
      !
      !=============================================================
      ! Wrapper function for solving a complex-valued linear systemt
      ! of equations Ax=b, where A(n,n) is a square matrix, b(n) is
      ! a known vector, and x(n) is a vector to be determined. Use
      ! LAPACK's "simple" driver ZGESV. Note that both A and X are
      ! overwritten on oputput!
      !=============================================================
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables:
      complex(8), intent(inout)  :: A(:, :), X(:, :)
      integer, intent(in), optional :: verb_
      ! Local variables:
      character(*), parameter :: myName = 'solLinSysV'
      integer                 :: n, info, ipiv(size(A, 1))
      integer                 :: verb
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      if (present(verb_)) then
         verb = verb_
      else
         verb = 0
      end if
      !
      n = size(A, 1)
      if (n /= size(A, 2)) then
         write (*, '(A,A)') myName, '> ERROR: size(A,1) != size(A,2)'
         STOP
      elseif (n /= size(X, 1)) then
         write (*, '(A,A)') myName, '> ERROR: size(A,2) !=  size(X,1)'
         STOP
      end if
      !
      call ZGESV(n, size(X, 2), A, n, ipiv, X, n, info)
      if (info /= 0) then
         write (*, '(A,A,i3)') &
            myName, '> ERROR: ZGESV returned non-zero INFO= ', info
         STOP
      end if
      !
   end subroutine solLinSysV
   !
   !******************************************************************
   !
   subroutine solLinSysVX(A, X, verb_)
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables:
      complex(8), intent(inout)  :: A(:, :), X(:, :)
      integer, intent(in), optional :: verb_
      ! Local variables:
      character(*), parameter :: myName = 'solLinSysVX'
      character(len=1)        :: equed
      integer                 :: n, info, ipiv(size(A, 1)), nrhs, verb
      !real                    :: t0,t1
      real(8)                 :: rcond
      real(8), dimension(size(A, 1)) :: r, c
      real(8), dimension(size(X, 2)) :: ferr, berr
      real(8), dimension(2*size(A, 1)) :: rwork
      complex(8), dimension(2*size(A, 1)) :: work
      complex(8), dimension(size(A, 1), size(A, 2)) :: AF
      complex(8), dimension(size(X, 1), size(X, 2)) :: B
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      if (present(verb_)) then
         verb = verb_
      else
         verb = 0
      end if
      !
      n = size(A, 1)
      if (n /= size(A, 2)) then
         write (*, '(A,A)') myName, '> ERROR: size(A,1) != size(A,2)'
         STOP
      elseif (n /= size(X, 1)) then
         write (*, '(A,A)') myName, '> ERROR: size(A,2) !=  size(X,1)'
         STOP
      end if
      nrhs = size(X, 2)
      B = X
      !
      call ZGESVX( &
         'E', 'N', &           ! FACT, TRANS
         n, nrhs, &           ! N, NRHS
         A, n, &              ! matrix A, LDA
         AF, n, ipiv, &       ! matrix AF, LDAF, IPIV
         equed, r, c, &       ! EQUED, R, C
         B, n, X, n, &        ! matrix B, LDB, matrix X, LDX
         rcond, ferr, berr, & ! RCOND, FERR, BERR
         work, rwork, info)  ! WORK, RWORK, and INFO
      !
      if (info /= 0) then
         write (*, '(A,A,es12.4E3)') myName, '> RWORK(1)= ', rwork(1)
         write (*, '(A,A,i6)') myName, '> ZGESVX returned non-zero INFO= ', info
         if (info < 0) then
            write (*, '(A,A,i6,A)') myName, '> ERROR in ', abs(info), 'th argument in ZGESVX'
            STOP
         elseif (abs(info) <= n) then
            write (*, '(A,A,i6,A)') myName, '> ERROR: U(i,i)= 0 for i= ', abs(info), &
               ' => singular matrix'
            write (*, '(A,A,es12.4E3,A,A)') &
               myname, '> ZGESVX returned RCOND= ', rcond, ' EQUED= ', equed
            STOP
         elseif (abs(info) == n + 1) then
            write (*, '(A,A)') myName, &
               '> WARNING: Non-singular matrix with RCOND < machine precision'
         else
            write (*, '(A,A)') myName, '> ERROR: Cannot interpret INFO'
            STOP
         end if
      end if
      !
      if (verb > 0) then
         write (*, '(A,A,es12.4E3,A,A)') &
            myname, '> ZGESVX returned RCOND= ', rcond, ' EQUED= ', equed
         write (*, '(A,A,2(1x,A,es12.4E3))') &
            myname, '> ZGESVX:', 'max(abs(FERR))= ', maxval(abs(ferr)), &
            'max(abs(BERR))= ', maxval(abs(berr))
      end if
      !
   end subroutine solLinSysVX
   !
   !******************************************************************
   !
   ! subroutine solLinSysVXX(A, X, verb_)
   !   !
   !   ! Supposed to use ZGESVXX, which appears to be missing in my system's LAPACK
   !   !
   !   !---------------------------------------------------
   !   ! Start of variable declarations.
   !   !---------------------------------------------------
   !   ! Passed variables:
   !   complex(8), intent(inout)  :: A(:,:), X(:,:)
   !   integer, intent(in), optional :: verb_
   !   ! Local variables:
   !   character(*), parameter :: myName='solLinSysVXX'
   !   character(len=1)        :: equed
   !   integer                 :: n,info,ipiv(size(A,1)), nrhs, verb
   !   real(8)                 :: rcond, rpvgrw, params(3)
   !   real(8), dimension(size(A,1)) :: r, c
   !   real(8), dimension(size(X,2)) :: berr
   !   real(8), dimension(2*size(A,1)) :: rwork
   !   real(8), dimension(size(X,2),3) :: err_bnds_norm, err_bnds_comp
   !   complex(8), dimension(2*size(A,1)) :: work
   !   complex(8), dimension(size(A,1), size(A,2)) :: AF
   !   complex(8), dimension(size(X,1), size(X,2)) :: B
   !   !---------------------------------------------------
   !   ! End of variable declarations. Directives start now
   !   !---------------------------------------------------
   !   !
   !   if(present(verb_)) then
   !      verb = verb_
   !   else
   !      verb = 0
   !   endif
   !   !
   !   n=size(A,1)
   !   if(n /= size(A,2)) then
   !      write(*,'(A,A)') myName,'> ERROR: size(A,1) != size(A,2)'
   !      STOP
   !   elseif(n /= size(X,1)) then
   !      write(*,'(A,A)') myName,'> ERROR: size(A,2) !=  size(X,1)'
   !      STOP
   !   endif
   !   nrhs = size(X,2)
   !   B = X
   !   !
   !   params(1:3) = -1
   !   !
   !   ! Why is ZGESVXX missing?!?
   !   call ZGESVX( &
   !        'E','N', &           ! FACT, TRANS
   !        n, nrhs, &           ! N, NRHS
   !        A, n, &              ! matrix A, LDA
   !        AF, n, ipiv, &       ! matrix AF, LDAF, IPIV
   !        equed, r, c, &       ! EQUED, R, C (none used here)
   !        B, n, X, n, &        ! matrix B, LDB, matrix X, LDX
   !        rcond, rpvgrw, berr, & ! RCOND, RPVGRW, BERR
   !        3, err_bnds_norm, err_bnds_comp, &
   !        3, params, &   ! NPARAMS and PARAMS(3)
   !        work, rwork, info )  ! WORK, RWORK, and INFO
   !   !
   !   if(info /= 0) then
   !      write(*,'(A,A,es12.4E3)') myName,'> RPVGRW= ', rpvgrw
   !      write(*,'(A,A,i6)') myName,'> ZGESVXX returned non-zero INFO= ',info
   !      if(info < 0) then
   !         write(*,'(A,A,i6,A)') &
   !              myName,'> ERROR in ',abs(info),'th argument in ZGESVXX'
   !         STOP
   !      elseif(abs(info) <= n) then
   !         write(*,'(A,A,i6,A)') myName,'> ERROR: U(i,i)= 0 for i= ',abs(info), &
   !              ' => singular matrix'
   !         write(*,'(A,A,es12.4E3,A,A)') &
   !              myname,'> ZGESVXX returned RCOND= ',rcond,' EQUED= ', equed
   !         STOP
   !      elseif(abs(info) <= n + nrhs) then
   !         write(*,'(A,A,i6,A)') myName,&
   !              '> WARNING: Solution to ',abs(info)-n,'th RHS not guaranteed'
   !      else
   !         write(*,'(A,A)') myName,'> ERROR: Cannot interpret INFO, stopping.. '
   !         STOP
   !      endif
   !   endif
   !   !
   !   if(verb > 0) then
   !      write(*,'(A,A,es12.4E3,A,A)') &
   !           myname,'> ZGESVXX returned RCOND= ',rcond,' EQUED= ', equed
   !      write(*,'(A,A,es12.4E3)') &
   !           myname,'> ZGESVXX: max(abs(BERR))= ',maxval(abs(berr))
   !   endif
   !   !
   ! end subroutine solLinSysVXX
   ! !
   ! !******************************************************************
   ! !
   ! subroutine solLinSysLS(A, X, verb_)
   !   !
   !   ! Seek least-squares solutions using ZGELS. Bad results.
   !   ! Note that X should store B on input, and then will be
   !   ! overwritten by the solution on output.
   !   !
   !   !---------------------------------------------------
   !   ! Start of variable declarations.
   !   !---------------------------------------------------
   !   ! Passed variables:
   !   complex(8), intent(inout)  :: A(:,:), X(:,:)
   !   integer, intent(in), optional :: verb_
   !   ! Local variables:
   !   character(*), parameter :: myName='solLinSysLS'
   !   integer                 :: n,info,rank,nrhs, lwork, verb
   !   complex(8), allocatable :: work(:)
   !   !---------------------------------------------------
   !   ! End of variable declarations. Directives start now
   !   !---------------------------------------------------
   !   !
   !   if(present(verb_)) then
   !      verb = verb_
   !   else
   !      verb = 0
   !   endif
   !   !
   !   n=size(A,1)
   !   if(n /= size(A,2)) then
   !      write(*,'(A,A)') myName,'> ERROR: size(A,1) != size(A,2)'
   !      STOP
   !   elseif(n /= size(X,1)) then
   !      write(*,'(A,A)') myName,'> ERROR: size(A,2) !=  size(X,1)'
   !      STOP
   !   endif
   !   nrhs = size(X,2)
   !   !
   !   if (verb > 0) write(*,'(A,A,A)', advance='NO') &
   !        myName,'> Querying for optimal LWORK'
   !   allocate(work(1))
   !   lwork = -1
   !   call ZGELS( 'N', & ! TRANS
   !        n, n, nrhs, & ! M, N, NRHS
   !        A, n, X, n, & ! A, LDA, B, LDB
   !        work, lwork, & ! WORK, LWORK
   !        info ) ! INFO
   !   !
   !   if(info /= 0) then
   !      if(verb > 0) write(*,*)
   !      write(*,'(A,A,i3,A)') myName, &
   !           '> ERROR: ZGELS returned non-zero INFO= ',info,' when querying'
   !      STOP
   !   endif
   !   !
   !   lwork = work(1)
   !   if(verb > 0) write(*,'(A,i9)') '= ', lwork
   !   !
   !   deallocate(work)
   !   allocate(work(lwork))
   !   !
   !   call ZGELS( 'N', & ! TRANS
   !        n, n, nrhs, & ! M, N, NRHS
   !        A, n, X, n, & ! A, LDA, B, LDB
   !        work, lwork, & ! WORK, LWORK
   !        info ) ! RWORK, INFO
   !   !
   !   deallocate(work)
   !   !
   !   if(info /= 0) then
   !      write(*,'(A,A,i3,A)') myName, &
   !           '> ERROR: ZGELS returned non-zero INFO= ',info,' when solving'
   !      STOP
   !   endif
   !   !
   ! end subroutine solLinSysLS
   ! !
   ! !******************************************************************
   ! !
   ! subroutine solLinSysLSS(A,X,verb_)
   !   !
   !   ! Seek solutions by singular value decomposition using ZGELSS.
   !   ! Currently very bad results.
   !   !
   !   !---------------------------------------------------
   !   ! Start of variable declarations.
   !   !---------------------------------------------------
   !   ! Passed variables:
   !   complex(8), intent(inout)  :: A(:,:), X(:,:)
   !   integer, intent(in), optional :: verb_
   !   ! Local variables:
   !   character(*), parameter :: myName='solLinSysLSS'
   !   integer                 :: n,info,rank,nrhs, lwork,verb
   !   real(8)                 :: rcond, s(size(A,1)), rwork(5*size(A,1))
   !   complex(8), allocatable :: work(:)
   !   !---------------------------------------------------
   !   ! End of variable declarations. Directives start now
   !   !---------------------------------------------------
   !   !
   !   if(present(verb_)) then
   !      verb = verb_
   !   else
   !      verb = 0
   !   endif
   !   !
   !   n=size(A,1)
   !   if(n /= size(A,2)) then
   !      write(*,'(A,A)') myName,'> ERROR: size(A,1) != size(A,2)'
   !      STOP
   !   elseif(n /= size(X,1)) then
   !      write(*,'(A,A)') myName,'> ERROR: size(A,2) !=  size(X,1)'
   !      STOP
   !   endif
   !   nrhs = size(X,2)
   !   !
   !   if(verb > 0)write(*,'(A,A,A)', advance='NO') &
   !        myName,'> Querying for optimal LWORK'
   !   allocate(work(1))
   !   lwork = -1
   !   call ZGELSS( &
   !        n, n, nrhs, & ! M, N, NRHS
   !        A, n, X, n, & ! A, LDA, B, LDB
   !        s, -1, rank, & ! S, RCOND [in], RANK
   !        work, lwork, & ! WORK, LWORK
   !        rwork, info ) ! RWORK, INFO
   !   !
   !   if(info /= 0) then
   !      if(verb > 0) write(*,*)
   !      write(*,'(A,A,i3,A)') myName, &
   !           '> ERROR: ZGELSS returned non-zero INFO= ',info,' when querying'
   !      STOP
   !   endif
   !   !
   !   lwork = work(1)
   !   deallocate(work); allocate(work(lwork))
   !   if(verb > 0) write(*,'(A,i9)') '= ', lwork
   !   !
   !   call ZGELSS( &
   !        n, n, nrhs, & ! M, N, NRHS
   !        A, n, X, n, & ! A, LDA, B, LDB
   !        s, -1, rank, & ! S, RCOND [in], RANK
   !        work, lwork, & ! WORK, LWORK
   !        rwork, info ) ! RWORK, INFO
   !   !
   !   deallocate(work)
   !   !
   !   if(info /= 0) then
   !      write(*,'(A,A,i3,A)') myName, &
   !           '> ERROR: ZGELSS returned non-zero INFO= ',info,' when solving'
   !      STOP
   !   elseif(verb > 0) then
   !      write(*,'(A,A,es12.4E3, A, i9)') &
   !           myName,'> ZGELSS returned RCOND=S(N)/S(1)= ',s(n)/s(1), &
   !           ' RANK= ', rank
   !   endif
   !   !
   ! end subroutine solLinSysLSS
   ! !
   ! !******************************************************************
   ! !
   ! subroutine solLinSysLSY(A,X,verb_)
   !   !
   !   ! Seek solutions by singular value decomposition using ZGELSY.
   !   !
   !   !---------------------------------------------------
   !   ! Start of variable declarations.
   !   !---------------------------------------------------
   !   ! Passed variables:
   !   complex(8), intent(inout)  :: A(:,:), X(:,:)
   !   integer, intent(in), optional :: verb_
   !   ! Local variables:
   !   character(*), parameter :: myName='solLinSysLSY'
   !   integer                 :: n,info,rank,nrhs,lwork,verb
   !   integer                 :: jpvt(size(A,1))
   !   real(8), parameter      :: rcond = 1.0d-16
   !   real(8)                 :: rwork(2*size(A,1))
   !   complex(8), allocatable :: work(:)
   !   !---------------------------------------------------
   !   ! End of variable declarations. Directives start now
   !   !---------------------------------------------------
   !   !
   !   if(present(verb_)) then
   !      verb = verb_
   !   else
   !      verb = 0
   !   endif
   !   !
   !   n=size(A,1)
   !   if(n /= size(A,2)) then
   !      write(*,'(A,A)') myName,'> ERROR: size(A,1) != size(A,2)'
   !      STOP
   !   elseif(n /= size(X,1)) then
   !      write(*,'(A,A)') myName,'> ERROR: size(A,2) !=  size(X,1)'
   !      STOP
   !   endif
   !   nrhs = size(X,2)
   !   !
   !   if(verb > 0)write(*,'(A,A,A)', advance='NO') &
   !        myName,'> Querying for optimal LWORK'
   !   allocate(work(1))
   !   lwork = -1
   !   call ZGELSY( &
   !        n, n, nrhs, & ! M, N, NRHS
   !        A, n, X, n, & ! A, LDA, B, LDB
   !        jpvt, rcond, rank, & ! JPVT, RCOND [in], RANK
   !        work, lwork, & ! WORK, LWORK
   !        rwork, info ) ! RWORK, INFO
   !   !
   !   if(info /= 0) then
   !      if(verb > 0) write(*,*)
   !      write(*,'(A,A,i3,A)') myName, &
   !           '> ERROR: ZGELSY returned non-zero INFO= ',info,' when querying'
   !      STOP
   !   endif
   !   !
   !   lwork = work(1)
   !   deallocate(work); allocate(work(lwork))
   !   if(verb > 0) write(*,'(A,i9)') '= ', lwork
   !   !
   !   call ZGELSY( &
   !        n, n, nrhs, & ! M, N, NRHS
   !        A, n, X, n, & ! A, LDA, B, LDB
   !        jpvt, rcond, rank, & ! S, RCOND [in], RANK
   !        work, lwork, & ! WORK, LWORK
   !        rwork, info ) ! RWORK, INFO
   !   !
   !   deallocate(work)
   !   !
   !   if(info /= 0) then
   !      write(*,'(A,A,i3,A)') myName, &
   !           '> ERROR: ZGELSY returned non-zero INFO= ',info,' when solving'
   !      STOP
   !   elseif(verb > 0) then
   !      write(*,'(A, A, i9)') myName,'> ZGELSY returned RANK= ', rank
   !   endif
   !   !
   ! end subroutine solLinSysLSY
   ! !
   ! !******************************************************************
   ! !
   ! subroutine solLinSysTSLS(A,X,verb_)
   !   !
   !   !---------------------------------------------------
   !   ! Start of variable declarations.
   !   !---------------------------------------------------
   !   ! Passed variables:
   !   complex(8), intent(inout)  :: A(:,:), X(:,:)
   !   integer, intent(in), optional :: verb_
   !   ! Local variables:
   !   character(*), parameter :: myName='solLinSysTSLS'
   !   integer                 :: n,info,nrhs,lwork,verb
   !   complex(8), allocatable :: work(:)
   !   !---------------------------------------------------
   !   ! End of variable declarations. Directives start now
   !   !---------------------------------------------------
   !   !
   !   if(present(verb_)) then
   !      verb = verb_
   !   else
   !      verb = 0
   !   endif
   !   !
   !   n=size(A,1)
   !   if(n /= size(A,2)) then
   !      write(*,'(A,A)') myName,'> ERROR: size(A,1) != size(A,2)'
   !      STOP
   !   elseif(n /= size(X,1)) then
   !      write(*,'(A,A)') myName,'> ERROR: size(A,2) !=  size(X,1)'
   !      STOP
   !   endif
   !   nrhs = size(X,2)
   !   !
   !   if(verb > 0)write(*,'(A,A,A)', advance='NO') &
   !        myName,'> Querying for optimal LWORK'
   !   allocate(work(1))
   !   lwork = -1
   !   call ZGETSLS( 'N', & ! TRANS
   !        n, n, nrhs, & ! M, N, NRHS
   !        A, n, X, n, & ! A, LDA, B, LDB
   !        work, lwork, info )  ! WORK, LWORK,  INFO
   !   !
   !   if(info /= 0) then
   !      if(verb > 0) write(*,*)
   !      write(*,'(A,A,i3,A)') myName, &
   !           '> ERROR: ZGETSLS returned non-zero INFO= ',info,' when querying'
   !      STOP
   !   endif
   !   !
   !   lwork = work(1)
   !   deallocate(work); allocate(work(lwork))
   !   if(verb > 0) write(*,'(A,i9)') '= ', lwork
   !   !
   !   call ZGETSLS( 'N', & ! TRANS
   !        n, n, nrhs, & ! M, N, NRHS
   !        A, n, X, n, & ! A, LDA, B, LDB
   !        work, lwork, info ) ! WORK, LWORK , INFO
   !   !
   !   deallocate(work)
   !   !
   !   if(info /= 0) then
   !      write(*,'(A,A,i3,A)') myName, &
   !           '> ERROR: ZGETSLS returned non-zero INFO= ',info,' when solving'
   !      STOP
   !   endif
   !   !
   ! end subroutine solLinSysTSLS
   !
end module linalg
