module miet
   !
   !==================================================================
   ! This module contains routines for calculating one-body T-matrices
   !
   implicit none
   !
   private
   !
   public :: calcMieTMat, calcMieIntCoeffs, calcStoutCoeffs
   !
contains
   !
   !==============================================================
   !
   subroutine calcMieTMat(x, s, tmat, zeropad_)
      !
      ! Calculate the Mie T-matrix (mieTmat) elements for a given
      ! x, s, and, nmax inferred from mieTmat dimensions. The
      ! T-matrix is diagonal, and the off-diagonal elements can be
      ! optionally overwritten (padded) with zeroes. Optionally,
      ! focus on a particular multipolarity band in the range
      ! [nband_(1),nband_(2)].
      !
      !---------------------------------------------------
      ! Start of variable declarations
      !---------------------------------------------------
      ! Passed variables
      complex(8), intent(in) :: x(:) ! size parameters (kOut*R)
      complex(8), intent(in) :: s(size(x)) ! rel. refractive indices (kIn/kOut)
      complex(8), intent(inout) :: tmat(:, :)
      logical, intent(in), optional :: zeropad_
      ! Local variables
      character(*), parameter :: myName = 'calcMieTMat'
      integer :: lmax, pmax, nmax, q, n, m, l, kmax
      complex(8), allocatable :: mie_coeffs(:, :, :)
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      if (size(tmat, 1) /= size(tmat, 2)) then
         write (*, '(A,A)') myname, '> ERROR: Supplied matrix not square'
         STOP
      else
         lmax = size(tmat, 1)
      end if
      !
      pmax = lmax/2
      nmax = int(sqrt(dble(pmax)))
      kmax = size(x) ! maximal shell count
      !
      if (2*nmax*(nmax + 2) /= lmax) then
         write (*, '(A,A)') myname, '> ERROR: 2*nmax*(nmax+2) /= lmax'
         STOP
      end if
      !
      if (present(zeropad_)) then
         if (zeropad_) then ! Pad off-diagonal elements with zeroes
            do l = 1, lmax
               tmat(1:l - 1, l) = 0
               tmat(l + 1:lmax, l) = 0
            end do
         end if
      end if
      !
      allocate (mie_coeffs(nmax, 2, kmax))
      mie_coeffs = 0
      call calcCoatMieCoeffs( &
         x=x(1:kmax), s=s(1:kmax), &
         gammas=mie_coeffs(1:nmax, 1, 1:kmax), &
         deltas=mie_coeffs(1:nmax, 2, 1:kmax))
      !
      l = 0
      do q = 1, 2
         do n = 1, nmax
            do m = -n, n
               l = l + 1
               tmat(l, l) = mie_coeffs(n, q, kmax)
            end do
         end do
      end do
      deallocate (mie_coeffs)
      !
   end subroutine calcMieTMat
   !
   subroutine calcMieCoeffs(x, s, gammas, deltas)
      !
      ! ============================================================
      ! Compute the Mie coefficients as defined by equations
      ! H.46 and H.47 on page 606 of Eric and Pablo's book.
      ! The coefficients are interpreted as magnetic and electric
      ! susceptibilities (Gamma_n and Delta_n, respectively) of
      ! the scattered field. Relation to Mishchenko's Mie
      ! coefficients: a_n = -Delta_n and b_n = -Gamma_n.
      !
      ! INPUT:
      ! -----
      ! x - dimensionless size parameter of the scattering sphere
      ! s - relative refractive index (rri = kIn/kOut)
      ! nmax - maximal value of the multipole index
      !
      ! OUTPUT:
      ! ------
      ! gammas(nmax) - magnetic susceptibilities
      ! deltas(nmax) - electric susceptibilities
      !============================================================
      !
      use swav, only: calcRiccatiBessels
      !
      !---------------------------------------------------
      ! Start of variable declarations
      !---------------------------------------------------
      ! Passed variables
      real(8), intent(in) :: x ! size parameter
      complex(8), intent(in) :: s ! relative refractive index
      complex(8), dimension(:), intent(inout) :: gammas
      complex(8), dimension(size(gammas)), intent(inout) :: deltas
      ! Local variables
      character(*), parameter :: myname = 'calcMieCoeffs'
      integer :: n, nmax
      complex(8) :: z, dum1, dum2, dum3, dum4
      complex(8), dimension(size(gammas)) :: psi1, xi1, dpsi1, dxi1
      complex(8), dimension(size(gammas)) :: psi2, xi2, dpsi2, dxi2
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      nmax = size(gammas)
      !
      z = cmplx(x, 0, 8)
      call calcRiccatiBessels(z, nmax, .true., psi1, dpsi1)
      call calcRiccatiBessels(z, nmax, .false., xi1, dxi1)
      z = x*s
      call calcRiccatiBessels(z, nmax, .true., psi2, dpsi2)
      call calcRiccatiBessels(z, nmax, .false., xi2, dxi2)
      !
      do n = 1, nmax
         dum1 = psi1(n)*dpsi2(n)
         dum2 = psi2(n)*dpsi1(n)
         dum3 = psi2(n)*dxi1(n)
         dum4 = xi1(n)*dpsi2(n)
         gammas(n) = (s*dum1 - dum2)/(dum3 - s*dum4)
         deltas(n) = (dum1 - s*dum2)/(s*dum3 - dum4)
      end do
      !
   end subroutine calcMieCoeffs
   !
   ! ============================================================
   !
   subroutine calcCoatMieCoeffs(x, s, gammas, deltas)
      !
      ! ============================================================
      ! Compute the Mie coefficients as defined by equations
      ! H.110 and H.113 on page 623 of Eric and Pablo's book.
      ! The coefficients are interpreted as magnetic and electric
      ! susceptibilities (Gamma_n and Delta_n, respectively) of
      ! the scattered field.
      !
      ! INPUT:
      ! -----
      ! x(kmax) - size parameter of the concentric spherical surfaces
      !           x_k = kOut_k*r_k, with the index k incereasing
      !           radially outwards (i.e. r_{k} > r_{k-1}).
      ! s(kmax) - relative refractive index for each surface
      !           s_k = kIn_k/kOut_k
      ! nmax - maximal value of the multipole index
      !
      ! OUTPUT:
      ! ------
      ! gammas(nmax,kmax) - magnetic susceptibilities
      ! deltas(nmax,kmax) - electric susceptibilities
      !============================================================

      !
      use swav, only: calcRiccatiBessels
      !
      !---------------------------------------------------
      ! Start of variable declarations
      !---------------------------------------------------
      ! Passed variables
      complex(8), intent(in) :: x(:) ! size parameters
      complex(8), intent(in) :: s(size(x)) ! relative refractive indices
      complex(8), intent(inout) :: gammas(:, :)
      complex(8), intent(inout) :: deltas(size(gammas, 1), size(gammas, 2))
      ! Local variables
      character(*), parameter :: myname = 'calcMieCoatCoeffs'
      integer :: n, nmax, k, kmax
      complex(8) :: z, zdum1, zdum2, N_nk, D_nk
      complex(8), dimension(size(gammas, 1)) :: gammas_km1, deltas_km1
      complex(8), dimension(size(gammas, 1)) :: psi1, xi1, dpsi1, dxi1
      complex(8), dimension(size(gammas, 1)) :: psi2, xi2, dpsi2, dxi2
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      if (size(gammas, 2) /= size(x)) then
         write (*, '(A,A)') myName, '> ERROR: size(gammas,2) /= size(x)'
         STOP
      else
         kmax = size(x)
         nmax = size(gammas, 1)
      end if
      !
      gammas_km1 = 0
      deltas_km1 = 0
      do k = 1, kmax
         !
         !write(*,'(A,i1,3(1x,es15.8E2))') 'TEST k,x(k),s(k)= ', k, x(k), s(k)
         !
         z = x(k)
         call calcRiccatiBessels(z, nmax, .true., psi1, dpsi1)
         call calcRiccatiBessels(z, nmax, .false., xi1, dxi1)
         z = x(k)*s(k)
         call calcRiccatiBessels(z, nmax, .true., psi2, dpsi2)
         call calcRiccatiBessels(z, nmax, .false., xi2, dxi2)
         !
         do n = 1, nmax
            !
            zdum1 = psi2(n) + gammas_km1(n)*xi2(n)
            zdum2 = dpsi2(n) + gammas_km1(n)*dxi2(n)
            N_nk = zdum1*dpsi1(n) - s(k)*psi1(n)*zdum2 ! H.111
            D_nk = s(k)*xi1(n)*zdum2 - zdum1*dxi1(n)   ! H.112
            gammas(n, k) = N_nk/D_nk                  ! H.110
            !
            zdum1 = dpsi2(n) + deltas_km1(n)*dxi2(n)
            zdum2 = psi2(n) + deltas_km1(n)*xi2(n)
            N_nk = psi1(n)*zdum1 - s(k)*zdum2*dpsi1(n) ! H.114
            D_nk = s(k)*zdum2*dxi1(n) - xi1(n)*zdum1   ! H.115
            deltas(n, k) = N_nk/D_nk                  ! H.113
            !
         end do
         !
         gammas_km1 = gammas(:, k)
         deltas_km1 = deltas(:, k)
         !
      end do
      !
   end subroutine calcCoatMieCoeffs
   !
   ! ============================================================
   !
   subroutine calcStoutCoeffs(x, rri, nmax, Cn, Dn)
      !
      ! ============================================================
      ! Compute the coefficients as defined by equation (50) in the
      ! Stout02 paper. Used to calculate absorption...
      !
      ! INPUT:
      ! -----
      ! x - dimensionless size parameter of the scattering sphere
      ! rri - relative refractive index
      ! nmax - maximal value of the multipole index
      !
      ! OUTPUT:
      ! ------
      ! Cn(nmax) - C coefficients
      ! Dn(nmax) - D coefficients
      !============================================================
      !
      use swav, only: calcRiccatiBessels
      !
      !---------------------------------------------------
      ! Start of variable declarations
      !---------------------------------------------------
      ! Passed variables
      complex(8), intent(in) :: x ! size parameter
      complex(8), intent(in) :: rri ! relative refractive index
      integer, intent(in) :: nmax ! maximal order
      real(8), dimension(nmax), intent(out) :: Cn, Dn
      ! Local variables
      character(*), parameter :: myname = 'calcStoutCoeffs'
      integer :: n
      complex(8) :: z, dum1, dum2, dum3, dum4
      complex(8), dimension(nmax) :: psi1, dpsi1, psi2, dpsi2
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      z = x
      call calcRiccatiBessels(z, nmax, .true., psi1, dpsi1)
      z = x*rri
      call calcRiccatiBessels(z, nmax, .true., psi2, dpsi2)
      !
      do n = 1, nmax
         !
         dum1 = conjg(psi2(n))*dpsi2(n)
         dum2 = psi2(n)*dpsi1(n)
         dum3 = dpsi2(n)*psi1(n)
         !
         dum4 = dum2 - rri*dum3
         dum4 = realpart(dum4*conjg(dum4))
         Cn(n) = realpart(cmplx(0, 1, 8)*rri*dum1)/dum4
         !
         dum4 = rri*dum2 - dum3
         dum4 = realpart(dum4*conjg(dum4))
         Dn(n) = realpart(cmplx(0, 1, 8)*conjg(rri)*dum1)/dum4
         !
      end do
      ! write(*,'(A,A,2(1x,es9.2E3))') &
      !      myname,'> min/max(lambdasM)=', &
      !      minval(abs(lambdasM)), maxval(abs(lambdasM))
      ! write(*,'(A,A,2(1x,es9.2E3))') &
      !      myname,'> min/max(lambdasN)=', &
      !      minval(abs(lambdasN)), maxval(abs(lambdasN))
      !
      return
      !
   end subroutine calcStoutCoeffs
   !
   subroutine calcMieIntCoeffs(a, k, scaCoeffs, intCoeffsReg, intCoeffsIrr, csAbs)
      !
      !=============================================================
      ! Calculate the regular and irregular VSW coefficients for the
      ! field in each concentric region of a (coated) Mie scatterer.
      ! The regular coefficients for the outer-most region inside
      ! the scatterer are calculated from the irregular coefficients
      ! (i.e. the scattered field) --- not the regular ones (the
      ! excitation field) --- in the host medium; and then downward
      ! recurrence is used to compute the regular coefficients in all
      ! other sub-surface regions. See eqns. H.117-H.123 of Eric and
      ! Pablo's book for the recurrence details, but beware that the
      ! concentric region and interface indexing is very confusing!
      !
      ! Explainer:
      ! x, s, and Mie coefficients are indexed in reference to
      ! interfaces, whereas the wave coefficients are indexed in
      ! reference to regions. In section H.5.2, index k refers
      ! to an interface as well as the adjacent exterior region.
      ! The indexing is zero-based, with the core region having k=0.
      ! The core has non-zero regular field associated with it, but
      ! the interface for k=0 is non-existend, so the corresponding
      ! Mie coefficients are manually set to zero.
      !
      ! INPUT:
      ! -----
      !
      ! a(1:) - radii of of concentric interfaces (increasing order)
      ! k(0:) - relative refractive index and kMedium
      ! scaCoeffs(:) - scattered field coefficients for the host medium
      !
      ! OUTPUT:
      ! ------
      !
      ! intCoeffsReg(:,:) - Regular field coefficients for each
      !                     concentric region inside the scatterer
      ! intCoeffsIrr(:,:) - Irregular field coefficients for each
      !                     concentric region inside the scatterer
      !
      !=============================================================
      !
      use swav, only: calcRiccatiBessels
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variable
      real(8), intent(in) :: a(:) ! interface radii in increasing order
      complex(8), intent(in) :: k(0:size(a)), scaCoeffs(:)
      complex(8), intent(out) :: intCoeffsReg(size(scaCoeffs), 0:size(a) - 1)
      complex(8), intent(out) :: intCoeffsIrr(size(scaCoeffs), 0:size(a) - 1)
      real(8), intent(out) :: csAbs(0:size(a) - 1)
      ! Note that the concentric region indices are made 0-based!
      ! Local variables
      character(*), parameter :: myName = 'calcMieIntCoeffs'
      integer :: p, n, m, l, pmax, nmax, ns
      complex(8) :: z, zdum1, zdum2, s(size(a)), x(size(a))
      complex(8), allocatable :: psi1(:), dpsi1(:), psi2(:), dpsi2(:)
      complex(8), allocatable :: xi1(:), dxi1(:), xi2(:), dxi2(:)
      complex(8), allocatable :: gammas(:, :), deltas(:, :)
      complex(8), parameter :: imu = cmplx(0, 1, kind(imu))
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      ns = size(a) ! shell count
      pmax = size(scaCoeffs)/2
      nmax = int(sqrt(real(pmax)))
      if (size(scaCoeffs) /= 2*nmax*(nmax + 2)) then
         write (*, '(A,A)') myName, '> ERROR: size(scaCoeffs) /= 2*nmax*(nmax+2)'
         STOP
      end if
      !
      ! Calculate Eric and Pablo's x_k and s_k on pg. 623
      x(1:ns) = a(1:ns)*k(1:ns)
      s(1:ns) = k(0:ns - 1)/k(1:ns)
      !
      allocate (gammas(nmax, 0:ns), deltas(nmax, 0:ns))
      gammas(:, 0) = 0 ! need this if ns=1
      deltas(:, 0) = 0
      allocate (psi1(nmax), dpsi1(nmax), xi1(nmax), dxi1(nmax))
      allocate (psi2(nmax), dpsi2(nmax), xi2(nmax), dxi2(nmax))
      call calcCoatMieCoeffs( &
         x=x, & ! radially increasing order
         s=s, &
         gammas=gammas(:, 1:ns), &
         deltas=deltas(:, 1:ns))
      !
      ! Treat scatterer's outer-most (and possibly the only) concentric region.
      ! Note that the regular internal field is calculated from the (irregular)
      ! scattered field in the host medium.
      z = x(ns)
      call calcRiccatiBessels(z, nmax, .true., psi1, dpsi1)
      call calcRiccatiBessels(z, nmax, .false., xi1, dxi1)
      z = x(ns)*s(ns)
      call calcRiccatiBessels(z, nmax, .true., psi2, dpsi2)
      call calcRiccatiBessels(z, nmax, .false., xi2, dxi2)
      z = imu*s(ns)
      p = 0
      do n = 1, nmax
         !
         zdum1 = (psi2(n) + gammas(n, ns - 1)*xi2(n))*dpsi1(n) &
                 - s(ns)*psi1(n)*(dpsi2(n) + gammas(n, ns - 1)*dxi2(n)) ! H.111
         zdum2 = psi1(n)*(dpsi2(n) + deltas(n, ns - 1)*dxi2(n)) &
                 - s(ns)*(psi2(n) + deltas(n, ns - 1)*xi2(n))*dpsi1(n)   ! H.114
         zdum1 = -z/zdum1
         zdum2 = z/zdum2
         !
         do m = -n, n
            p = p + 1
            intCoeffsReg(p, ns - 1) = zdum1*scaCoeffs(p)
            intCoeffsReg(p + pmax, ns - 1) = zdum2*scaCoeffs(p + pmax)
         end do
         !
      end do
      !
      ! Now treat scatterer's sub-surface layers (if present)
      !
      icoat: do l = ns - 1, 1, -1
         !
         !if(l < ns) then
         z = x(l)
         call calcRiccatiBessels(z, nmax, .true., psi1, dpsi1)
         call calcRiccatiBessels(z, nmax, .false., xi1, dxi1)
         z = x(l)*s(l)
         call calcRiccatiBessels(z, nmax, .true., psi2, dpsi2)
         call calcRiccatiBessels(z, nmax, .false., xi2, dxi2)
         !endif
         !
         ! Use equations H.117-H.123
         p = 0
         z = imu*s(l) ! H.119 and H.122
         do n = 1, nmax
            !
            !
            zdum1 = (psi2(n) + gammas(n, l - 1)*xi2(n))*dxi1(n) &
                    - s(l)*xi1(n)*(dpsi2(n) + gammas(n, l - 1)*dxi2(n)) ! H.120
            zdum2 = s(l)*(psi2(n) + deltas(n, l - 1)*xi2(n))*dxi1(n) &
                    - xi1(n)*(dpsi2(n) + deltas(n, l - 1)*dxi2(n))      ! H.123
            !
            zdum1 = z/zdum1 ! H.118
            zdum2 = z/zdum2 ! H.121
            !
            do m = -n, n
               p = p + 1
               intCoeffsReg(p, l - 1) = zdum1*intCoeffsReg(p, l) ! H.117     1st
               intCoeffsReg(p + pmax, l - 1) = zdum2*intCoeffsReg(p + pmax, l) ! 2nd
               intCoeffsIrr(p, l) = gammas(n, l)*intCoeffsReg(p, l) ! H.108    1st
               intCoeffsIrr(p + pmax, l) = deltas(n, l)*intCoeffsReg(p + pmax, l)! 2nd
            end do
            !
         end do
         !
      end do icoat
      !
      intCoeffsIrr(:, 0) = 0 ! no irregular field inside central core
      !
      ! Calculate partial absorptions using equaiton (29) in MackowskiAM90
      ! The paper uses a different basis, so adapt the equaion by changing:
      ! chi_n -> xi_n
      ! (a_n,b_n) -> intCoeffsReg
      ! (c_n,d_n) -> intCoeffsIrr
      csAbs = 0
      do l = ns, 1, -1
         z = k(l - 1)*a(l)
         call calcRiccatiBessels(z, nmax, .true., psi1, dpsi1)
         call calcRiccatiBessels(z, nmax, .false., xi1, dxi1)
         p = 0; z = 0
         do n = 1, nmax
            do m = -n, n
               p = p + 1
               z = z + &
                   (intCoeffsReg(p + pmax, l - 1)*dpsi1(n) + &
                    intCoeffsIrr(p + pmax, l - 1)*dxi1(n))* &
                   conjg(intCoeffsReg(p + pmax, l - 1)*psi1(n) + &
                         intCoeffsIrr(p + pmax, l - 1)*xi1(n))
               z = z - &
                   (intCoeffsReg(p, l - 1)*psi1(n) + &
                    intCoeffsIrr(p, l - 1)*xi1(n))* &
                   conjg(intCoeffsReg(p, l - 1)*dpsi1(n) + &
                         intCoeffsIrr(p, l - 1)*dxi1(n))
            end do
         end do
         z = imu*k(ns)/k(l - 1)*z
         csAbs(l - 1) = realpart(z)/realpart(k(ns))**2
      end do
      !
      deallocate (gammas, deltas)
      deallocate (xi1, dxi1, psi1, dpsi1, xi2, dxi2, psi2, dpsi2)
      !
   end subroutine calcMieIntCoeffs
   !
end module miet
