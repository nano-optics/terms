module sphmsv

!
   ! ==============================================================
   ! This module contains routines for calculating Stokes incident
   ! vector, Stokes phase matrix and Scattering matrix from a T-matrix
   ! Using Amos' toms644.f to calculate the Bessel and Hankel functions
   ! Last modified: 09/12/2019.
   !
   ! References:
   ! ----------
   !
   ! [1] Mishchenko, Travis and Lacis, "Scattering, absorption, and
   !     emission of light by small particles" (CUP, 2002).
   !
   !
   ! List of routines:
   ! ----------------
   ! calcStokesIncVec> Evaluates the Stokes incident vector
   !
   ! calcStokesPhaseMat> Evaluates the Stokes Phase Matrix for
   !                     the given incident and scattering angles
   !
   ! calcScatMat> Evaluates the Scattering matrix using the T-matrix for
   !              the given incident and scattering angles
   ! ==============================================================
   implicit none
   !
   private
   !
   public :: pi, calcScatMat, calcStokesPhaseMat, calcStokesIncVec
   !
   real(8), parameter :: pi = acos(-1.0d0)  ! pi
   complex(8), parameter :: imu = (0.0d0, 1.0d0) ! the imaginary unit i
   real(8), parameter :: eps0 = 8.8541878128d-12 !vacuum permittivity
   real(8), parameter :: mu0 = 4.0*pi*1.0d-7    !vacuum permibility
contains

   subroutine calcStokesIncVec(ehost_, ipwDirn_, ipwAmpl_, Stokes_Vec, verb_)
      use swav, only: xyz2rtp, calcVTxyz2rtp
      !------------------------------------------------------------------
      !This subroutine calculates the stokes vector for the incident field
      !based on the formula 1.54 on page 17 of Mishchenko's book
      !------------------------------------------------------------------
      ! Passed variables
      complex(8), intent(in) :: ipwAmpl_(3)
      real(8), intent(in) :: ipwDirn_(3), ehost_
      real(8), intent(out) :: Stokes_Vec(4)
      integer, intent(in), optional :: verb_
      ! Local variables
      character(*), parameter :: myName = 'calcStokesIncVec'
      complex(8) ::  E0(3)
      real(8) ::  transform(3, 3), cth, rtp(3)
      integer :: verb

      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      verb = 0
      if (present(verb_)) verb = verb_
      call xyz2rtp(ipwDirn_, rtp, cth)
      ! Transform ipwE0 from cartesian to spherical polars
      call calcVTxyz2rtp(rtp, transform)
      E0(:) = matmul(transform(:, :), ipwAmpl_)

      Stokes_Vec(1) = E0(2)*dconjg(E0(2)) + E0(3)*dconjg(E0(3))
      Stokes_Vec(2) = E0(2)*dconjg(E0(2)) - E0(3)*dconjg(E0(3))
      Stokes_Vec(3) = -2*realpart(E0(2)*dconjg(E0(3)))
      Stokes_Vec(4) = -2*imagpart(E0(2)*dconjg(E0(3)))

      Stokes_Vec = (0.5d0)*sqrt(ehost_*eps0/mu0)*Stokes_Vec

      if (verb > 0) write (*, '(A,A,/)') myname, '> Done!'
   end subroutine calcStokesIncVec
   !------------------------------------------------------
   subroutine calcStokesPhaseMat(SMat, Z, verb_)

      !----------------------------------------------------------------------------------------------------
      !This subroutine calculates the Stokes phase matrix in terms of amplitude scattering matrix elements
      !The equations are based on formulae (2.106)-(2.121) of Mishchenko's book on page 51.
      !-----------------------------------------------------------------------------------------------------
      ! Passed variables

      complex(8), dimension(2, 2), intent(in) :: SMat(:, :)
      real(8), intent(out) :: Z(4, 4)
      integer, intent(in), optional :: verb_
      ! Local variables
      character(*), parameter :: myName = 'calcStokesPhaseMat'
      integer :: verb

      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      verb = 0
      if (present(verb_)) verb = verb_

      Z = 0; 
      Z(1, 1) = dble(0.5d0*(sum((ABS(SMat(:, :)))**2)))
      Z(1, 2) = dble(0.5d0*((ABS(SMat(1, 1)))**2 - (ABS(SMat(1, 2)))**2 &
                            + (ABS(SMat(2, 1)))**2 - (ABS(SMat(2, 2)))**2))
      Z(1, 3) = dble(-realpart(SMat(1, 1)*dconjg(SMat(1, 2)) + SMat(2, 2)*dconjg(SMat(2, 1))))
      Z(1, 4) = dble(-imagpart(SMat(1, 1)*dconjg(SMat(1, 2)) - SMat(2, 2)*dconjg(SMat(2, 1))))
      Z(2, 1) = dble(0.5d0*((ABS(SMat(1, 1)))**2 + (ABS(SMat(1, 2)))**2 &
                            - (ABS(SMat(2, 1)))**2 - (ABS(SMat(2, 2)))**2))
      Z(2, 2) = dble(0.5d0*((ABS(SMat(1, 1)))**2 - (ABS(SMat(1, 2)))**2 &
                            - (ABS(SMat(2, 1)))**2 + (ABS(SMat(2, 2)))**2))
      Z(2, 3) = dble(-realpart(SMat(1, 1)*dconjg(SMat(1, 2)) - SMat(2, 2)*dconjg(SMat(2, 1))))
      Z(2, 4) = dble(-imagpart(SMat(1, 1)*dconjg(SMat(1, 2)) + SMat(2, 2)*dconjg(SMat(2, 1))))
      Z(3, 1) = dble(-realpart(SMat(1, 1)*dconjg(SMat(2, 1)) + SMat(2, 2)*dconjg(SMat(1, 2))))
      Z(3, 2) = dble(-realpart(SMat(1, 1)*dconjg(SMat(2, 1)) - SMat(2, 2)*dconjg(SMat(1, 2))))
      Z(3, 3) = dble(realpart(SMat(1, 1)*dconjg(SMat(2, 2)) + SMat(1, 2)*dconjg(SMat(2, 1))))
      Z(3, 4) = dble(imagpart(SMat(1, 1)*dconjg(SMat(2, 2)) + SMat(2, 1)*dconjg(SMat(1, 2))))
      Z(4, 1) = dble(-imagpart(SMat(2, 1)*dconjg(SMat(1, 1)) + SMat(2, 2)*dconjg(SMat(1, 2))))
      Z(4, 2) = dble(-imagpart(SMat(2, 1)*dconjg(SMat(1, 1)) - SMat(2, 2)*dconjg(SMat(1, 2))))
      Z(4, 3) = dble(imagpart(SMat(2, 2)*dconjg(SMat(1, 1)) - SMat(1, 2)*dconjg(SMat(2, 1))))
      Z(4, 4) = dble(realpart(SMat(2, 2)*dconjg(SMat(1, 1)) - SMat(1, 2)*dconjg(SMat(2, 1))))

      if (verb > 0) write (*, '(A,A,/)') myname, '> Done!'

   end subroutine calcStokesPhaseMat
   !---------------------------------------------------------

   subroutine calcScatMat(tmat, k, spwDirn_, ipwDirn_, SMat, verb_)

      use swav, only: xyz2rtp, calcWignerd0andMore, nm2p

      !
      ! Calculate the Scattering matrix (calcScatMat) elements for a given
      ! Tmatrix using the formulae (5.11)-(5.17) Mishchenko's book.
      !
      !---------------------------------------------------
      ! Start of variable declarations
      !---------------------------------------------------
      ! Passed variables
      complex(8), intent(in) :: tmat(:, :)
      real(8), intent(in) :: k, spwDirn_(3), ipwDirn_(3)
      complex(8), dimension(2, 2), intent(out) :: SMat
      integer, intent(in), optional :: verb_
      ! Local variables
      character(*), parameter :: myName = 'calcScatMat'

      integer :: lmax, pmax, nmax, m, n, nprim, mprim, pprim, p  !kmax,
      complex(8) :: alpha
      real(8) :: x_sca, x_inc, rtp_inc(3), rtp_sca(3)
      integer :: l, ll
      real(8), dimension(0:size(tmat, 1)/2) :: d, pi_sca, pi_inc, tau_sca, tau_inc
      integer :: verb
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      verb = 0
      if (present(verb_)) verb = verb_
      SMat = 0
      if (size(tmat, 1) /= size(tmat, 2)) then
         write (*, '(A,A)') myname, '> ERROR: Supplied matrix not square'
         STOP
      else
         lmax = size(tmat, 1)
      end if
      !
      pmax = lmax/2
      nmax = int(sqrt(dble(pmax)))
      !kmax = size(x) ! maximal shell count
      !
      if (2*nmax*(nmax + 2) /= lmax) then
         write (*, '(A,A)') myname, '> ERROR: 2*nmax*(nmax+2) /= lmax'
         STOP
      end if
      call xyz2rtp(ipwDirn_, rtp_inc, x_inc)
      call xyz2rtp(spwDirn_, rtp_sca, x_sca)

      call calcWignerd0andMore(x_sca, pmax, d, pi_sca, tau_sca)  !x_sca=cos(theta_sca),
      call calcWignerd0andMore(x_inc, pmax, d, pi_inc, tau_inc)  !x_inc=cos(theta_inc),

      do n = 1, nmax
         do m = -n, n
            do nprim = 1, nmax
               do mprim = -nprim, nprim

                  p = n*(n + 1) + m
                  pprim = nprim*(nprim + 1) + mprim

                  call nm2p(n, m, l)
                  call nm2p(nprim, mprim, ll)

                  alpha = (imu**(nprim - n - 1))*((-1)**(m + mprim))*((2*n + 1)*(2*nprim + 1)&
                   &/(n*(n + 1)*nprim*(nprim + 1)))**(0.5d0)

                  ! call calcPi(m,n,tetasca,pi_sca, mprim,nprim,tetainc,pi_inc)
                  ! call calcTau(m,n,tetasca,tau_sca, mprim,nprim,tetainc,tau_inc)

                  SMat(1, 1) = SMat(1, 1) &
                               + alpha*(tmat(p, pprim)*pi_sca(l)*pi_inc(ll) &
                                        + tmat(p, pmax + pprim)*tau_sca(l)*pi_inc(ll) &
                                        + tmat(pmax + p, pprim)*pi_sca(l)*tau_inc(ll) &
                                        + tmat(pmax + p, pmax + pprim)*tau_sca(l)*tau_inc(ll)) &
                               *ZEXP(imu*(m*rtp_sca(3) - mprim*rtp_inc(3)))

                  SMat(1, 2) = SMat(1, 2) &
                               + alpha*(tmat(p, pprim)*pi_sca(l)*tau_inc(ll) &
                                        + tmat(p, pmax + pprim)*tau_sca(l)*tau_inc(ll) &
                                        + tmat(pmax + p, pprim)*pi_sca(l)*pi_inc(ll) &
                                        + tmat(pmax + p, pmax + pprim)*tau_sca(l)*pi_inc(ll)) &
                               *ZEXP(imu*(m*rtp_sca(3) - mprim*rtp_inc(3)))

                  SMat(2, 1) = SMat(2, 1) &
                               + alpha*(tmat(p, pprim)*tau_sca(l)*pi_inc(ll) &
                                        + tmat(p, pmax + pprim)*pi_sca(l)*pi_inc(ll) &
                                        + tmat(pmax + p, pprim)*tau_sca(l)*tau_inc(ll) &
                                        + tmat(pmax + p, pmax + pprim)*pi_sca(l)*tau_inc(ll)) &
                               *ZEXP(imu*(m*rtp_sca(3) - mprim*rtp_inc(3)))

                  SMat(2, 2) = SMat(2, 2) &
                               + alpha*(tmat(p, pprim)*tau_sca(l)*tau_inc(ll) &
                                        + tmat(p, pmax + pprim)*pi_sca(l)*tau_inc(ll) &
                                        + tmat(pmax + p, pprim)*tau_sca(l)*pi_inc(ll) &
                                        + tmat(pmax + p, pmax + pprim)*pi_sca(l)*pi_inc(ll)) &
                               *ZEXP(imu*(m*rtp_sca(3) - mprim*rtp_inc(3)))
               end do
            end do
         end do
      end do
      SMat(1, 1) = (1/k)*SMat(1, 1)
      SMat(1, 2) = (-imu/k)*SMat(1, 2)
      SMat(2, 1) = (imu/k)*SMat(2, 1)
      SMat(2, 2) = (1/k)*SMat(2, 2)
      if (verb > 0) write (*, '(A,A,/)') myname, '> Done!'
   end subroutine calcScatMat

end module sphmsv
