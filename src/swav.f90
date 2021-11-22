module swav
   !
   ! ==============================================================
   ! This module contains routines for calculating and transforming
   ! scalar (SSWs) and vector spherical waves (VSWs). It depends on
   ! Amos' toms644.f to calculate the Bessel and Hankel functions
   ! using recurrence. Last modified: 18/03/2019.
   !
   ! References:
   ! ----------
   ! [1] Stout, Auger and Lfait, J. Mod. Optics 49, 2129-2152 (2002).
   ! [2] Mishchenko, Travis and Lacis, "Scattering, absorption, and
   !     emission of light by small particles" (CUP, 2002).
   ! [3] Chew, J. Electron Waves Applic. 6, 133 (1992).
   !
   ! List of routines:
   ! ----------------
   !
   ! 1-calcVTACs> Evaluates the vector translation-addition
   !            coefficients.
   !
   ! 2-calcSTACs> Calculates the scalar translation-addition
   !            coefficients.
   !
   ! 3-calcVTACsAxial> Calculates the simpler VTACs along the z-axis.
   !
   ! 4-calcSTACsAxial> Calculates the simpler STACs along the z-axis.
   !
   ! 5-calcVSWs> Evaluates the vector spherical waves.
   !
   ! 6-calcSSWs> Evaluates the scalar spherical waves.
   !
   ! 7-calcJCoeffsPW> Evaluates coefficients of scatterer-centred
   !                  VSW expansions for an incident planewave.
   !
   ! 8-offsetCoeffsPW> Translates the VSW coefficients of a plane-
   !                   wave to another origin.
   !
   ! 9-calcWignerBigD> Evaluates the Wigner D-funcctions.
   !
   ! 10-calcWignerLittled> Evaluates the Wigner d-funcctions.
   !
   ! 11-calcWignerd0andMore> Evaluates the Wigner d-functions for n=0,
   !                      and the derivative functions pi and tau.
   !
   ! 12-calcRiccatiBessels> Evaluates the Riccati-Bessel functions psi
   !                     or xi, and their derivatives.
   !
   ! 13-calcSphBessels> Evaluates the spherical Bessel/Hankel funs.
   !
   ! 14-xyz2rtp> (x,y,z) -> (r,theta,phi) for a point.
   !
   ! 15-rtp2xyz> (r,theta,phi) -> (x,y,z) for a point.
   !
   ! 16-calcVTrtp2xyz> Compute the transformation matrix for
   !                (E_r, E_t, E_p) -> (E_x, E_y, E_z)
   !
   ! 17-calcVTxyz2rtp> Compute the transformation matrix for
   !                (E_x, E_y, E_z) -> (E_r, E_t, E_p)
   ! 18-nm2p> Compute a generalised index l=n(n+1)+m, for a unique (n,m)
   !
   ! 19-p2nm> Compute unique (n,m) from a given composite index p
   !
   ! 20-nm2pv2> Some recurrences are defined only for m >= 0, in which case
   !            we shall use a second version of the composite index
   !            p_v2:=n*(n+1)/2+m. This subroutins calculates p for (n,m)
   !
   ! 21-testPmax> Tests pmax for commensurability, i.e. is pmax=nnmax*(nmax+2) and nmax=mmax? If not, then stop program
   ! ==============================================================
   !
   implicit none
   !
   private
   !
   public :: pi, tpi, xyz2rtp, calcVSWs, calcVTACs, calcVTACsAxial, &
             calcRiccatiBessels, calcSphBessels, calcWignerBigD, &
             calcCoeffsPW, offsetCoeffsPW, calcJCoeffsPW, calcVTxyz2rtp, &
             calcWignerd0andMore, nm2p, calcVTrtp2xyz, testPmax, calcAbsMat, calcLamMat
   !, &
   !calcScatMat, calcStokesPhaseMat, calcStokesIncVec
   !

   real(8), parameter :: pi = acos(-1.0d0)  ! pi
   real(8), parameter :: hpi = 0.5d0*pi     ! pi/2
   real(8), parameter :: tpi = 2.0d0*pi     ! 2*pi
   real(8), parameter :: fpi = 4.0d0*pi     ! 4*pi
   complex(8), parameter :: imu = (0.0d0, 1.0d0) ! the imaginary unit i
   !real(8), parameter :: eps0=8.8541878128d-12
   !real(8), parameter :: mu0= 4.0*pi*1.0d-7
contains
   !
   subroutine calcVTACs(r0, k, regt, vtacs)
      !
      ! ============================================================
      ! Compute the irregular (if regt=FALSE) normalised vector
      ! translation-addition coefficients alpha_{nu,mu;n,m} or
      ! the regular (if regt=TRUE) beta_{nu,mu;n,m}, as defined
      ! in Appendix B of Stout02, for a given k*r0 and
      !
      ! 1 <= nu <= nmax, -nu <= mu <= nu;
      ! 1 <= n <= nmax, -n <= m <= n.
      !
      ! Note that beta=Reg(alpha), where Reg() corresponds to
      ! replacing the spherical Hankel functions h_n with the
      ! spherical Bessel functions j_n. For convenience,
      ! we use a composite index l(n,m) = n*(n+1)+m
      !
      ! 1 <= l:=n*(n+1)+m <= pmax:=nmax*(nmax+2),
      !
      ! which is spanned twice to accounting for both the M_nm
      ! and the N_nm wave components.
      !
      ! INPUT:
      ! ------
      ! r0(3) - relative position vector [REAL]
      ! k     - wave number [COMPLEX]
      ! regt  - regular or not [LOGICAL]
      !
      ! IN/OUTPUT:
      ! ------
      ! vtacs(1:2*pmax,1:2*pmax) - the coefficients [COMPLEX]
      ! ============================================================
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables
      real(8), dimension(3), intent(in) :: r0
      complex(8), intent(in) :: k
      logical, intent(in) :: regt
      complex(8), intent(inout) :: vtacs(:, :)
      ! Local variables
      integer :: nmax, pmax, n, m, l, nu, mu, lambda, i
      real(8) :: numuterm1, numuterm2, numuterm3, numuterm4, numuterm5
      real(8) :: nmterm1, nmterm2, nnuterm1, nnuterm2
      complex(8) :: z
      complex(8) :: scoeff(0:size(vtacs, 1)/2, 0:size(vtacs, 2)/2)
      character(*), parameter :: myname = 'calcVTACs'
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      ! Check the passed matrix and infer pmax with nmax
      if (size(vtacs, 1) /= size(vtacs, 2)) then
         write (*, '(A,A)') myname, '> ERROR: Passed matrix not square'
         STOP
      else
         pmax = size(vtacs, 1)/2
         nmax = int(sqrt(dble(pmax)))
         i = 2*nmax*(nmax + 2)
         if (i /= size(vtacs, 1)) then
            if (i < size(vtacs, 1)) then
               write (*, '(A,A)') myname, '> WARNING: 2*nmax*(nmax+2) < size(vtacs,1)'
            else
               write (*, '(A,A)') myname, '> ERROR: 2*nmax*(nmax+2) > size(vtacs,1)'
               STOP
            end if
         end if
      end if
      !
      vtacs = 0 ! Just in case
      !
      call calcSTACs(r0, k, pmax, regt, scoeff)
      !
      lambda = 0
      nuloop: do nu = 1, nmax
         muloop: do mu = -nu, nu
            lambda = lambda + 1
            numuterm1 = sqrt(dble((nu - mu)*(nu + mu + 1)))
            numuterm2 = sqrt(dble((nu + mu)*(nu - mu + 1)))
            numuterm3 = sqrt(dble((nu - mu)*(nu + mu)))
            numuterm4 = sqrt(dble((nu - mu)*(nu - mu - 1)))
            numuterm5 = sqrt(dble((nu + mu)*(nu + mu - 1)))
            l = 0
            nloop: do n = 1, nmax
               nnuterm1 = 0.5d0/sqrt(dble(nu*(nu + 1)*n*(n + 1)))
               nnuterm2 = nnuterm1*sqrt(dble(2*nu + 1)/dble(2*nu - 1))
               mloop: do m = -n, n
                  l = l + 1
                  nmterm1 = sqrt(dble((n - m)*(n + m + 1)))
                  nmterm2 = sqrt(dble((n + m)*(n - m + 1)))
                  !
                  ! A_{nu,mu;n,m} block
                  z = 2*mu*m*scoeff(lambda, l)
                  if (abs(mu + 1) <= nu .and. abs(m + 1) <= n) then
                     z = z + nmterm1*numuterm1*scoeff(lambda + 1, l + 1)
                  end if
                  if (abs(mu - 1) <= nu .and. abs(m - 1) <= n) then
                     z = z + nmterm2*numuterm2*scoeff(lambda - 1, l - 1)
                  end if
                  z = nnuterm1*z
                  vtacs(lambda, l) = z
                  vtacs(lambda + pmax, l + pmax) = z
                  !
                  ! B_{nu,mu;n,m} block
                  z = 0
                  if (abs(mu) < nu) then
                     call nm2p(nu - 1, mu, i)
                     z = z + 2*m*numuterm3*scoeff(i, l)
                  end if
                  if (abs(mu + 1) < nu .and. abs(m + 1) <= n) then
                     call nm2p(nu - 1, mu + 1, i)
                     z = z + nmterm1*numuterm4*scoeff(i, l + 1)
                  end if
                  if (abs(mu - 1) < nu .and. abs(m - 1) <= n) then
                     call nm2p(nu - 1, mu - 1, i)
                     z = z - nmterm2*numuterm5*scoeff(i, l - 1)
                  end if
                  z = -imu*nnuterm2*z
                  vtacs(lambda + pmax, l) = z
                  vtacs(lambda, l + pmax) = z
                  !
               end do mloop
            end do nloop
            !
         end do muloop
      end do nuloop
      !
   end subroutine calcVTACs
   !
   subroutine calcSTACs(r0, k, pmax, regt, scoeff)
      !
      ! ============================================================
      ! Compute the normalised scalar translation-addition
      ! coefficients alphas_{nu,mu;n,m} or betas_{nu,mu;n,m} for
      !
      ! 0 <= nu <= nmax, -nu <= mu <= nu;
      ! 0 <= n <= nmax, -n <= m <= n;
      !
      ! using the recurrence described in Appendix C of Stout02
      ! [J. Mod. Optics 49, 2129 (2002)] for a given k*r0.
      !
      ! Note that betas=Reg(alphas), where Reg() corresponds to
      ! replacing the spherical Hankel functions h_n with the
      ! spherical Bessel functions j_n. Also, the paper Chew92
      ! [J. Electron Waves Applic. 6, 133 (1992)] is very helpful
      ! for understanding the implementation.
      !
      ! For convenience, we use a composite index l(n,m) when m
      ! can be +ve and -ve, and l_{v2}(n,m) for non-negative m:
      !
      ! 0 <= l=n*(n+1)+m <= pmax=nmax*(nmax+2);
      ! 0 <= l_{v2}=n*(n+1)/2+m <= l_{v2}^{max} = nmax*(nmax+3)/2.
      !
      ! INPUT:
      ! ------
      ! r0(3) - relative position vector [REAL]
      ! k     - wave-vector amplitude [COMPLEX]
      ! pmax  - maximal composite index [INTEGER]
      ! regt  - 'take regular part' or not [LOGICAL]
      !
      ! OUTPUT:
      ! ------
      ! scoeff(0:pmax,0:pmax) - coefficients [COMPLEX]
      !
      ! Output correponds to the scalar translation-addition
      ! coefficients alphas (irregular, for regt=FALSE) or betas
      ! (regular, for regt=TRUE) as defined in Appendix C of
      ! Stout02. Recall: pmax=nmax*(nmax+2).
      ! ============================================================
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables
      real(8), dimension(3), intent(in) :: r0
      complex(8), intent(in) :: k
      integer, intent(in) :: pmax
      logical, intent(in) :: regt
      complex(8), dimension(0:pmax, 0:pmax), intent(out) :: scoeff
      ! Local variables
      character(*), parameter :: myName = 'calcSTACs'
      integer :: n, m, l, nmax, pmax2, nu, mu, numax, lambda, lambdamax
      integer :: nu2, ndum1, ndum2, npm, nmm, ipass
      integer :: l1, l1v2, l2v2, lv2, lam
      real(8) :: r, theta, phi, ctheta, rtp(3), ddum
      complex(8) :: z
      complex(8), dimension(:), allocatable :: expo, ang, bes
      complex(8), dimension(:, :), allocatable :: c
      real(8), dimension(:), allocatable :: ap, am, bp, bm, d
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      ! Check the maximal value of the composite index and get nmax
      call testPmax(myname, pmax, nmax)
      !
      scoeff(:, :) = cmplx(0, 0, 8) ! just in case
      !
      ! The recurrence will be applied (twice) to non-negative m, so
      ! the composite index of some arrays is bounded above by pmax2:
      pmax2 = nmax*(nmax + 3)/2
      !
      ! The recurrence involves extraneous elements (see Chew92),
      ! so the maximum value of nu needs to be doubled.
      numax = nmax + nmax + 1
      lambdamax = numax*(numax + 2)
      !
      ! Allocate all internal arrays
      allocate (expo(-numax:numax), d(0:lambdamax), ang(0:lambdamax))
      allocate (bes(0:numax), c(0:lambdamax, 0:pmax2))
      allocate (am(0:lambdamax), ap(0:lambdamax))
      allocate (bm(0:lambdamax), bp(0:lambdamax))
      !
      ! Convert r0 from cartesians to spherical polars
      call xyz2rtp(r0, rtp, ctheta)
      r = rtp(1)
      z = k*r ! complex argument for spherical Bessel functions
      theta = rtp(2) ! between 0 and Pi
      phi = rtp(3) ! between 0 and 2*Pi
      !
      ! Compute angular quantities required for eqn. C3 of Stout02.
      do mu = -numax, numax
         expo(mu) = exp(imu*mu*phi)
      end do
      call calcWignerd0andMore(ctheta, lambdamax, d)
      !
      ! Compute the auxiliary coefficients a_nm^+/- and b_nm^+/-
      ! defined in eqn. C2, as well as the radially-INdependent
      ! part of the quantities in C3 of Stout02.
      do nu = 0, numax
         nu2 = 2*nu
         ndum1 = nu2 + 1
         ddum = sqrt(dble(ndum1))
         ndum2 = ndum1*(nu2 - 1)
         ndum1 = ndum1*(nu2 + 3)
         do mu = -nu, nu
            call nm2p(nu, mu, lambda)
            npm = nu + mu
            nmm = nu - mu
            ap(lambda) = -sqrt(dble((npm + 1)*(nmm + 1))/dble(ndum1))
            am(lambda) = sqrt(dble((npm)*(nmm))/dble(ndum2))
            bp(lambda) = sqrt(dble((npm + 2)*(npm + 1))/dble(ndum1))
            bm(lambda) = sqrt(dble((nmm)*(nmm - 1))/dble(ndum2))
            ! Get alpha0/bhn=beta0/bjn defined in equation C.3
            ang(lambda) = (-1)**(nu)*ddum*d(lambda)*expo(-mu)
         end do
      end do
      !
      ! Apply the recurrence on n to build up c_{nu,mu;n,m} in two
      ! passes: 1st for m >= 0; 2nd for m < 0. See eqn. C1 of Stout02
      passes: do ipass = 1, 2
         !
         call calcSphBessels(z, numax, regt, bes)
         !
         ! Now apply recursion in eqn. C1 of Stout02
         mloop: do m = 0, nmax
            !
            ! Referring to Fig. 2 of Chew92, we first step along the
            ! diagonal and compute c_{nu,mu; m,m} for all nu and mu
            !
            if (m == 0) then ! initialise using eqn. C3 of Stout02.
               do nu = 0, numax
                  do mu = -nu, nu
                     call nm2p(nu, mu, lambda)
                     c(lambda, 0) = ang(lambda)*bes(nu)
                  end do
               end do
            else ! n = m > 0 => use bottom eqn. C1
               call nm2pv2(m, m, lv2)
               call nm2p(m - 1, m - 1, l1)
               call nm2pv2(m - 1, m - 1, l1v2)
               do nu = 0, numax - m ! see Fig.3 in Chew92
                  do mu = -nu, nu
                     call nm2p(nu, mu, lambda)
                     c(lambda, lv2) = 0 ! initialise
                     if (nu > 0 .and. abs(mu - 1) < nu) then
                        call nm2p(nu - 1, mu - 1, lam)
                        c(lambda, lv2) = c(lambda, lv2) + bp(lam)*c(lam, l1v2)
                     end if
                     call nm2p(nu + 1, mu - 1, lam)
                     c(lambda, lv2) = c(lambda, lv2) + bm(lam)*c(lam, l1v2)
                     c(lambda, lv2) = c(lambda, lv2)/bp(l1)
                  end do
               end do
            end if
            !
            ! Now compute c_{nu,mu;n,m} for n > m (and all nu and mu).
            !
            nloop: do n = m + 1, nmax
               call nm2pv2(n, m, lv2)
               call nm2p(n - 1, m, l1)
               call nm2pv2(n - 1, m, l1v2)
               call nm2pv2(n - 2, m, l2v2)
               do nu = 0, numax - n ! see Fig.3 in Chew92
                  do mu = -nu, nu
                     call nm2p(nu, mu, lambda)
                     c(lambda, lv2) = 0 ! initialise
                     if (n - 2 >= m) then
                        c(lambda, lv2) = c(lambda, lv2) - am(l1)*c(lambda, l2v2)
                     end if
                     if (nu - 1 >= abs(mu)) then
                        call nm2p(nu - 1, mu, lam)
                        c(lambda, lv2) = c(lambda, lv2) + ap(lam)*c(lam, l1v2)
                     end if
                     call nm2p(nu + 1, mu, lam)
                     c(lambda, lv2) = c(lambda, lv2) + am(lam)*c(lam, l1v2)
                     c(lambda, lv2) = c(lambda, lv2)/ap(l1)
                  end do
               end do
            end do nloop
            !
         end do mloop
         !
         ! Assign final values to the output scoeff. Could do this
         ! inside mloop and nloop, but post-recurrence asignment
         ! seems cleaner for debugging purposes.
         mloop2: do m = 0, nmax
            nloop2: do n = m, nmax
               call nm2pv2(n, m, lv2)
               if (ipass == 1) then
                  call nm2p(n, m, l)
                  do nu = 0, nmax
                     do mu = -nu, nu
                        call nm2p(nu, mu, lambda)
                        scoeff(lambda, l) = c(lambda, lv2)
                     end do
                  end do
               elseif (ipass == 2 .and. m > 0) then ! use eqn. C4 from Stout02
                  call nm2p(n, -m, l)
                  do nu = 0, nmax
                     do mu = -nu, nu
                        call nm2p(nu, mu, lambda)
                        call nm2p(nu, -mu, lam)
                        if (regt) then
                           ndum1 = (-1)**(mu - m)
                        else
                           !ndum1 = -(-1)**(n+nu+mu-m)
                           ndum1 = (-1)**(n + nu + mu - m)
                        end if
                        scoeff(lambda, l) = ndum1*dconjg(c(lam, lv2))
                     end do
                  end do
               end if
            end do nloop2
         end do mloop2
         !
         ! For the second pass, change the radial argument
         ! as prescribed in eqn. C4 of Stout02.
         if (regt) then
            z = dconjg(z)
         else
            z = -dconjg(z)
         end if
         !
      end do passes
      !
      ! Dellocate internal arrays.
      deallocate (d, expo, ang, bes)
      deallocate (c, am, ap, bm, bp)
      !
   end subroutine calcSTACs
   !
   !-------------------------------------------------------------------
   ! ======= SIMPLER ROUTINES FOR AXIAL TRANSLATION START HERE  =======
   !-------------------------------------------------------------------
   !
   subroutine calcVTACsAxial(r0, k, pmax, regt, flip, vtacs, mqn_)
      !
      ! ============================================================
      ! Compute the irregular (if regt=FALSE) normalised vector
      ! translation-addition coefficients alpha_{nu,mu;n,m}  or
      ! the regular (if regt=TRUE) beta_{nu,mu;n,m}, as defined
      ! in Appendix B of Stout02, for a given k*r0 and
      !
      ! 1 <= nu <= nmax, -nu <= mu <= nu;
      ! 1 <= n <= nmax, m = mu.
      !
      ! Note that beta=Reg(alpha), where Reg() corresponds to
      ! replacing the spherical Hankel functions h_n with the
      ! spherical Bessel functions j_n. For convenience,
      ! we use a composite index l(n,m):
      !
      ! 0 <= l:=n*(n+1)+m <= 2*pmax:=2*nmax*(nmax+2),
      !
      ! where the prefactor of 2 is due to accounting for both
      ! the M_nm and the N_nm wave components.
      !
      ! INPUT:
      ! ------
      ! r0    - the z-axial displacelement distance [+ve REAL]
      ! k     - wave-vector amplitude [COMPLEX]
      ! pmax  - maximal composite index [INTEGER]
      ! regt  - 'take regular part' or not [LOGICAL]
      ! flip  -
      ! mqn_  - change from qnm to mqn indexing [LOGICAL]
      !
      ! OUTPUT:
      ! ------
      ! vtacs(1:2*pmax,1:2*pmax) - the coefficients [COMPLEX]
      !                            with (m,n,q) indexing to produce
      !                            block-diagonal form.
      ! ============================================================
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables
      real(8), intent(in) :: r0
      complex(8), intent(in) :: k
      integer, intent(in) :: pmax
      logical, intent(in) :: regt, flip
      logical, intent(in), optional :: mqn_
      complex(8), dimension(2*pmax, 2*pmax), intent(out) :: vtacs
      ! Local variables
      integer :: n, m, l, nmax, nu, mu, lambda, lam1, lo, ln, lnu
      real(8) :: numuterm1, numuterm2, numuterm3, &
                 numuterm4, numuterm5
      real(8) :: nmterm1, nmterm2
      real(8) :: nnuterm1, nnuterm2
      complex(8) :: z
      complex(8), dimension(0:pmax, 0:pmax) :: scoeff
      character(*), parameter :: myname = 'calcVTACs'
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      ! Check the maximal value of the composite index and nmax
      call testPmax(myname, pmax, nmax)
      !
      vtacs = 0 ! Just in case
      !
      call calcSTACsAxial(r0, k, pmax, regt, flip, scoeff)
      !
      ! Can reorder loops and accumulate (m,q,n) index while
      ! using the nm2p function to compute the (q,n,m) index
      lo = 1
      muloop: do mu = -nmax, nmax
         lnu = lo
         nuloop: do nu = max(1, abs(mu)), nmax
            call nm2p(nu, mu, lambda)
            numuterm1 = sqrt(dble((nu - mu)*(nu + mu + 1)))
            numuterm2 = sqrt(dble((nu + mu)*(nu - mu + 1)))
            numuterm3 = sqrt(dble((nu - mu)*(nu + mu)))
            numuterm4 = sqrt(dble((nu - mu)*(nu - mu - 1)))
            numuterm5 = sqrt(dble((nu + mu)*(nu + mu - 1)))
            !
            ln = lo
            nloop: do n = max(1, abs(mu)), nmax
               nnuterm1 = 0.5d0/sqrt(dble(nu*(nu + 1)*n*(n + 1)))
               nnuterm2 = nnuterm1*sqrt(dble(2*nu + 1)/dble(2*nu - 1))
               !
               m = mu
               !
               call nm2p(n, m, l)
               nmterm1 = sqrt(dble((n - m)*(n + m + 1)))
               nmterm2 = sqrt(dble((n + m)*(n - m + 1)))
               !
               ! A_{nu,mu;n,m} block
               z = 2*mu*m*scoeff(lambda, l)
               if (abs(mu + 1) <= nu .and. abs(m + 1) <= n) then
                  z = z + nmterm1*numuterm1*scoeff(lambda + 1, l + 1)
               end if
               if (abs(mu - 1) <= nu .and. abs(m - 1) <= n) then
                  z = z + nmterm2*numuterm2*scoeff(lambda - 1, l - 1)
               end if
               z = nnuterm1*z
               vtacs(lnu, ln) = z
               vtacs(lnu + 1, ln + 1) = z
               !
               ! B_{nu,mu;n,m} block
               z = 0
               if (abs(mu) < nu) then
                  call nm2p(nu - 1, mu, lam1)
                  z = z + 2*m*numuterm3*scoeff(lam1, l)
               end if
               if (abs(mu + 1) < nu .and. abs(m + 1) <= n) then
                  call nm2p(nu - 1, mu + 1, lam1)
                  z = z + nmterm1*numuterm4*scoeff(lam1, l + 1)
               end if
               if (abs(mu - 1) < nu .and. abs(m - 1) <= n) then
                  call nm2p(nu - 1, mu - 1, lam1)
                  z = z - nmterm2*numuterm5*scoeff(lam1, l - 1)
               end if
               z = -imu*nnuterm2*z
               vtacs(lnu + 1, ln) = z
               vtacs(lnu, ln + 1) = z
               !
               ln = ln + 2
            end do nloop
            !
            lnu = lnu + 2
         end do nuloop
         lo = lo + 2*(nmax - max(1, abs(m)) + 1)
      end do muloop
      !
   end subroutine calcVTACsAxial
   !
   subroutine calcSTACsAxial(r0, k, pmax, regt, flip, stacs)
      !
      ! ============================================================
      ! Compute the normalised scalar translation-addition
      ! coefficients alphas_{nu,mu;n,m} or betas_{nu,mu;n,m} for
      !
      ! 0 <= nu <= nmax, -nu <= mu <= nu;
      ! 0 <= n <= nmax, m = mu
      !
      ! using the recurrence described in Appendix C of Stout02
      ! [J. Mod. Optics 49, 2129 (2002)] for a given k*r0, where
      ! r0 is a non-negative displacement along the z-axis.
      !
      ! Note that betas=Reg(alphas), where Reg() corresponds to
      ! replacing the spherical Hankel functions h_n with the
      ! spherical Bessel functions j_n. Also, the paper Chew92
      ! [J. Electron Waves Applic. 6, 133 (1992)] is very helpful
      ! for understanding the implementation.
      !
      ! For convenience, we use a composite index l(n,m) when m
      ! can be +ve and -ve, and l_{v2}(n,m) for non-negative m:
      !
      ! 0 <= l=n*(n+1)+m <= pmax=nmax*(nmax+2);
      ! 0 <= l_{v2}=n*(n+1)/2+m <= l_{v2}^{max} = nmax*(nmax+3)/2.
      !
      ! INPUT:
      ! ------
      ! r0    - displacement distance [REAL]
      ! k     - wave-vector amplitude [COMPLEX]
      ! pmax  - maximal composite index [INTEGER]
      ! regt  - 'take regular part' or not [LOGICAL]
      !
      ! OUTPUT:
      ! ------
      ! stacs(0:pmax,0:pmax) - coefficients [COMPLEX]
      !
      ! Output correponds to the scalar translation-addition
      ! coefficients alphas (irregular, for regt=FALSE) or betas
      ! (regular, for regt=TRUE) as defined in Appendix C of
      ! Stout02. Recall: pmax=nmax*(nmax+2).
      ! ============================================================
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables
      real(8), intent(in) :: r0
      complex(8), intent(in) :: k
      integer, intent(in) :: pmax
      logical, intent(in) :: regt
      logical, intent(in) :: flip
      complex(8), dimension(0:pmax, 0:pmax), intent(out) :: stacs
      ! Local variables
      character(*), parameter :: myName = 'calcSTACsAxial'
      integer :: n, m, l, nmax, pmax2, nu, mu, numax, lambda, lambdamax
      integer :: nu2, ndum1, ndum2, npm, nmm, ipass
      integer :: l1, l1v2, l2v2, lv2, lam
      real(8) :: r, ddum
      complex(8) :: z
      complex(8), dimension(:), allocatable :: bes
      complex(8), dimension(:, :), allocatable :: c
      real(8), dimension(:), allocatable :: ap, am, bp, bm, ang
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      ! Check the maximal value of the composite index and get nmax
      call testPmax(myname, pmax, nmax)
      !
      stacs(:, :) = cmplx(0, 0, 8) ! just in case
      !
      ! The recurrence will be applied (twice) to non-negative m, so
      ! the composite index of some arrays is bounded above by pmax2:
      pmax2 = nmax*(nmax + 3)/2
      !
      ! The recurrence involves extraneous elements (see Chew92),
      ! so the maximum value of nu needs to be doubled.
      numax = nmax + nmax + 1
      lambdamax = numax*(numax + 2)
      !
      ! Allocate all internal arrays
      allocate (ang(0:numax))
      allocate (bes(0:numax), c(0:lambdamax, 0:pmax2))
      allocate (am(0:lambdamax), ap(0:lambdamax))
      allocate (bm(0:lambdamax), bp(0:lambdamax))
      !
      r = r0
      z = k*r ! complex argument for spherical Bessel functions
      !
      ! Compute angular quantities required for eqn. C3 of Stout02.
      !  >>> This can be simplified for axial translation !!! <<<
      !call calcWignerd0andMore(1.0d0, lambdamax, d)
      !
      ! Compute the auxiliary coefficients a_nm^+/- and b_nm^+/-
      ! defined in eqn. C2, as well as the radially-INdependent
      ! part of the quantities in C3 of Stout02.
      do nu = 0, numax
         nu2 = 2*nu
         ndum1 = nu2 + 1
         ddum = sqrt(dble(ndum1))
         ndum2 = ndum1*(nu2 - 1)
         ndum1 = ndum1*(nu2 + 3)
         do mu = -nu, nu
            call nm2p(nu, mu, lambda)
            npm = nu + mu
            nmm = nu - mu
            ap(lambda) = -sqrt(dble((npm + 1)*(nmm + 1))/dble(ndum1))
            am(lambda) = sqrt(dble((npm)*(nmm))/dble(ndum2))
            bp(lambda) = sqrt(dble((npm + 2)*(npm + 1))/dble(ndum1))
            bm(lambda) = sqrt(dble((nmm)*(nmm - 1))/dble(ndum2))
         end do
         ! Get alpha0/bhn=beta0/bjn defined in equation C.3,
         ! which are non-zero only when mu=0 for axial translation
         ang(nu) = (-1)**(nu)*ddum
         if (flip) ang(nu) = ang(nu)*(-1)**nu
      end do
      !
      !stop
      c(:, :) = 0
      !
      ! Apply the recurrence on n to build up c_{nu,mu;n,m} in two
      ! passes: 1st for m >= 0; 2nd for m < 0. See eqn. C1 of Stout02
      passes: do ipass = 1, 2
         !
         call calcSphBessels(z, numax, regt, bes)
         !
         ! Now apply recursion in eqn. C1 of Stout02
         mloop: do m = 0, nmax
            !
            ! Referring to Fig. 2 of Chew92, we first step along the
            ! diagonal and compute c_{nu,mu; m,m} for all nu and mu
            !
            if (m == 0) then ! initialise using eqn. C3 of Stout02.
               mu = m
               do nu = abs(mu), numax
                  call nm2p(nu, mu, lambda)
                  c(lambda, 0) = ang(nu)*bes(nu)
               end do
            else ! n = m > 0 => use bottom eqn. C1
               call nm2pv2(m, m, lv2)
               call nm2p(m - 1, m - 1, l1)
               call nm2pv2(m - 1, m - 1, l1v2)
               mu = m
               do nu = abs(mu), numax - m ! see Fig.3 in Chew92
                  call nm2p(nu, mu, lambda)
                  c(lambda, lv2) = 0 ! initialise
                  if (nu > 0 .and. abs(mu - 1) < nu) then
                     call nm2p(nu - 1, mu - 1, lam)
                     c(lambda, lv2) = c(lambda, lv2) + bp(lam)*c(lam, l1v2)
                  end if
                  call nm2p(nu + 1, mu - 1, lam)
                  c(lambda, lv2) = c(lambda, lv2) + bm(lam)*c(lam, l1v2)
                  c(lambda, lv2) = c(lambda, lv2)/bp(l1)
               end do
            end if
            !
            ! Now compute c_{nu,mu;n,m} for n > m (and all nu and mu).
            !
            nloop: do n = m + 1, nmax
               call nm2pv2(n, m, lv2)
               call nm2p(n - 1, m, l1)
               call nm2pv2(n - 1, m, l1v2)
               call nm2pv2(n - 2, m, l2v2)
               mu = m
               do nu = max(0, abs(mu)), numax - n ! see Fig.3 in Chew92
                  call nm2p(nu, mu, lambda)
                  c(lambda, lv2) = 0 ! initialise
                  if (n - 2 >= m) then
                     c(lambda, lv2) = c(lambda, lv2) - am(l1)*c(lambda, l2v2)
                  end if
                  if (nu - 1 >= abs(mu)) then
                     call nm2p(nu - 1, mu, lam)
                     c(lambda, lv2) = c(lambda, lv2) + ap(lam)*c(lam, l1v2)
                  end if
                  call nm2p(nu + 1, mu, lam)
                  c(lambda, lv2) = c(lambda, lv2) + am(lam)*c(lam, l1v2)
                  c(lambda, lv2) = c(lambda, lv2)/ap(l1)
               end do
            end do nloop
            !
         end do mloop
         !
         ! Assign final values to the output stacs. Could do this
         ! inside mloop and nloop, but post-recurrence asignment
         ! seems cleaner for debugging purposes.
         mloop2: do m = 0, nmax
            nloop2: do n = m, nmax
               call nm2pv2(n, m, lv2)
               if (ipass == 1) then
                  call nm2p(n, m, l)
                  mu = m
                  do nu = abs(mu), nmax
                     call nm2p(nu, mu, lambda)
                     stacs(lambda, l) = c(lambda, lv2)
                  end do
               elseif (ipass == 2 .and. m > 0) then ! use eqn. C4 from Stout02
                  call nm2p(n, -m, l)
                  mu = -m ! Minus sign should really be here!!!
                  do nu = abs(mu), nmax
                     call nm2p(nu, mu, lambda)
                     call nm2p(nu, -mu, lam)
                     if (regt) then
                        ndum1 = (-1)**(mu - m)
                     else
                        ndum1 = (-1)**(n + nu + mu - m)
                     end if
                     stacs(lambda, l) = ndum1*dconjg(c(lam, lv2))
                  end do
               end if
            end do nloop2
         end do mloop2
         !
         ! For the second pass, change the radial argument
         ! as prescribed in eqn. C4 of Stout02.
         if (regt) then
            z = dconjg(z)
         else
            z = -dconjg(z)
         end if
         !
      end do passes
      !
      ! Dellocate internal arrays.
      deallocate (ang, bes)
      deallocate (c, am, ap, bm, bp)
      !
   end subroutine calcSTACsAxial
   !
   subroutine calcVSWs(r, k, pmax, regt, cart, waves, wavesB)
      !
      ! ============================================================
      ! Evaluate (at r) the normalised vector spherical waves, M_nm
      ! and N_nm for l=n(n+1)+m <= pmax, as defined by equation A4
      ! in Appendix A of Stout02 [J. Mod. Optics 49, 2129 (2002)].
      !
      ! INPUT:
      ! ------
      ! r(3)  - Cartiesian coordinate of a point in 3D [REAL]
      ! k     - wave-vector amplitude [COMPLEX]
      ! pmax  - maximal composite index [INTEGER]
      ! regt  - 'take regular part' or not [LOGICAL]
      ! cart  - convert to cartesian coordinates or not [LOGICAL]
      !
      ! OUTPUT:
      ! ------
      ! waves(2*pmax,3) - [COMPLEX] elements of the abstract column
      !                   vector defined in eqn.B1 of Stout02, here
      !                   taking the form of a (2*pmax by 3) matrix.
      !wavesB(2*pmax,3) - is the same as waves, only by changing the
      ! position of the elements for magnetic field multiplied by -ik
      ! ============================================================
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables
      real(8), dimension(3), intent(in) :: r
      complex(8), intent(in) :: k
      integer, intent(in) :: pmax
      logical, intent(in) :: regt, cart
      complex(8), dimension(3, 2*pmax), intent(out) :: waves
      complex(8), dimension(3, 2*pmax), intent(out), optional :: wavesB
      ! Local variables
      character(*), parameter :: myname = 'calcVSWs'
      integer :: n, m, l, nmax
      real(8) :: rtp(3), cth, npref1, npref2, s_nm, u_nm
      real(8) :: transform(3, 3)
      real(8), dimension(0:pmax) :: d, pi, tau
      complex(8) :: z, y_nm(3), x_nm(3), z_nm(3)
      complex(8), allocatable :: expo(:), bes(:), dbes(:)
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      ! write(*,*)'waves',waves
      !waves=cmplx(0.0d0,0.0d0)
      !wavesB=cmplx(0.0d0,0.0d0)
      ! Check the maximal value of the composite index
      call testPmax(myname, pmax, nmax)
      !
      allocate (expo(-nmax:nmax), bes(nmax), dbes(nmax))
      !
      call xyz2rtp(r, rtp, cth) ! only need the cos(theta)
      call calcWignerd0andMore(cth, pmax, d, pi, tau)
      ! Compute the azimuth-angle exponentials
      do m = -nmax, nmax
         expo(m) = exp(imu*m*rtp(3))
      end do
      ! Compute the derivative of xi_n(kr) w.r.t. kr
      z = k*rtp(1) ! z=kr    !!! PROBLEMS when z=0 !!!
      if (CDABS(z) < 1.0e-30) then
         z = z + 1.0e-30 !! avoid division by 0 below by shifting z slightly away from origin
      end if
      call calcRiccatiBessels(z, nmax, regt, bes, dbes)
      bes(1:nmax) = bes(1:nmax)/z ! Get spherical Bess/Hank
      if (regt .eqv. .false.) then
         !write(*,*)'bes',bes
      end if
      nloop: do n = 1, nmax
         npref1 = sqrt(dble((2*n + 1)*n*(n + 1))/fpi)
         npref2 = npref1/dble(n*(n + 1))
         mloop: do m = -n, n
            call nm2p(n, m, l)
            y_nm(1) = npref1*d(l)*expo(m)  !note: ynm here is sqrt(n(n+1))*ynm in eq. A1
            y_nm(2:3) = 0
            s_nm = (-1)**m*npref2*tau(l)
            u_nm = (-1)**m*npref2*pi(l)
            !
            x_nm(1) = y_nm(2)    !note: x_nm here is = -xnm in eq. A1
            x_nm(2) = u_nm*imu*expo(m)
            x_nm(3) = -s_nm*expo(m)
            !
            z_nm(1) = x_nm(1)
            z_nm(2) = -x_nm(3)
            z_nm(3) = x_nm(2)
            !
            ! First compute M_nm and then N_nm
            waves(1:3, l) = bes(n)*x_nm(1:3) ! M_nm
            waves(1:3, l + pmax) = (bes(n)*y_nm(1:3) + dbes(n)*z_nm(1:3))/z ! N_nm
            if (present(wavesB)) then
               wavesB(1:3, l) = -imu*k*(bes(n)*y_nm(1:3) + dbes(n)*z_nm(1:3))/z ! N_nm, Note:for vector harmonics of B, -ik is multiplied
               wavesB(1:3, l + pmax) = -imu*k*bes(n)*x_nm(1:3) ! M_nm, Note:for vector harmonics of B, -ik is multiplied
            end if

         end do mloop
      end do nloop

      if (cart) then
         transform(1, 1:3) = (/sin(rtp(2))*cos(rtp(3)), &
                               cos(rtp(2))*cos(rtp(3)), &
                               -sin(rtp(3))/)
         transform(2, 1:3) = (/sin(rtp(2))*sin(rtp(3)), &
                               cos(rtp(2))*sin(rtp(3)), &
                               cos(rtp(3))/)
         transform(3, 1:3) = (/cos(rtp(2)), -sin(rtp(2)), 0.0d0/)
         waves = matmul(transform, waves)
         if (present(wavesB)) wavesB = matmul(transform, wavesB)
      end if
      !
      deallocate (expo, bes, dbes)

   end subroutine calcVSWs
   !
   !subroutine calcSSWs(xyz,k,pmax,regt,psi)
   !
   ! ============================================================
   ! Evaluate (at xyz) the spherical waves psi_nm as defined by
   ! equation 13a in Chew92 [J. Electron Waves Applic. 6, 133
   ! (1992)]. Note: l:=n(n+1)+m <= pmax. Also, Chew includes the
   ! Condon-Shortley phase (-1)^m in the definition of the
   ! spherical harmonics in eqn. 3b, implying that the Legendre
   ! polynomials in that equation do not contain the (-1)^m.
   ! In this routine, the (-1)^m is included in the Legendre
   ! polynomials obtained from the Wigner d-function.
   !
   ! INPUT:
   ! ------
   ! xyz(3) - Cartiesian coordinate of a point in 3D [REAL]
   ! k      - wave-vector amplitude [COMPLEX]
   ! pmax   - maximal composite index [INTEGER]
   ! regt   - 'take regular part' or not [LOGICAL]
   !
   ! OUTPUT:
   ! ------
   ! psi(0:pmax) - [COMPLEX] elements of the the spherical waves psi_nm as defined by
   !                         equation 13a in Chew92
   ! ============================================================
   !
   !---------------------------------------------------
   ! Start of variable declarations.
   !---------------------------------------------------
   ! Passed variables
   ! real(8), dimension(3), intent(in) :: xyz
   ! complex(8), intent(in) :: k
   ! integer, intent(in) :: pmax
   ! logical, intent(in) :: regt
   ! complex(8), dimension(0:pmax), intent(out) :: psi
   ! ! Local variables
   ! character(*), parameter :: myname='calc'
   ! integer :: nmax,n,m,l
   ! real(8) :: rtp(3), dum, d(0:pmax)
   ! complex(8), allocatable :: bess(:), expo(:)
   !---------------------------------------------------
   ! End of variable declarations. Directives start now
   !---------------------------------------------------
   !
   ! call testPmax(myname,pmax,nmax)
   !allocate(bess(0:nmax), expo(-nmax:nmax))
   !
   ! call xyz2rtp(xyz,rtp,dum) ! dum = cos(theta)
   !call calcWignerd0andMore(dum,pmax,d)
   !call calcSphBessels(k*rtp(1),nmax,regt,bess)
   ! do m=-nmax,nmax
   !    expo(m) = exp(imu*m*rtp(3))
   ! enddo
   !
   ! do n=0,nmax
   !    dum=sqrt(dble(2*n+1)/fpi)
   !    do m=-n,n
   !       call nm2p(n,m,l)
   !       psi(l) = bess(n) * dum * d(l) * expo(m) !(-1)^m in d(l)
   !       !write(*,'(A,2(1x,i2),2(1x,es17.10e2))') 'psi>',n,m,psi(l)
   !    enddo
   ! enddo
   !
   ! deallocate(bess,expo)
   !
   !end subroutine calcSSWs

!call calcJCoeffsPW( &
   !              ipwE0 = ipwE0(:,i), &
   !             kVec = hostK_G*ipwDirn(:,i), &
   !             ipwCoeffsJ = cJ_(:,1,i), &
   !             xyz = geometry(1:3,:) &
   !             )
   !
   subroutine calcJCoeffsPW(ipwE0, kVec, xyz, ipwCoeffsJ)
      !
      ! ============================================================
      ! Translate the supplied ipwCoeffs to different centres, whose
      ! coordinates are supplied in xyz(3,nc), and the incident
      ! plane wave direction is specified by real kVec.
      ! Followed equation 38 of Stout02 paper.
      ! ============================================================
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables
      complex(8), intent(in) :: ipwE0(3)
      real(8), intent(in) :: kVec(3), xyz(:, :)
      complex(8), intent(inout) :: ipwCoeffsJ(:)
      ! Local variables
      character(*), parameter :: myName = 'calcJCoeffsPW'
      integer :: j, js, lmax, nscat
      complex(8), allocatable :: ipwCoeffs(:)
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      if (size(xyz, 1) /= size(kVec)) then
         write (*, '(A,A)') myName, '> ERROR: size(xyz,1) /= size(kVec)'
         stop
      end if
      !
      nscat = size(xyz, 2)
      lmax = size(ipwCoeffsJ)/nscat
      if (size(ipwCoeffsJ) /= lmax*nscat) then
         write (*, '(A,A)') myName, '> ERROR: size(ipwCoeffsJ) /= lmax*nscat'
         stop
      end if
      !
      allocate (ipwCoeffs(lmax))
      call calcCoeffsPW( &
         ipwE0=ipwE0, &
         ipwDirn=kVec/sqrt(dot_product(kVec, kVec)), &
         ipwCoeffs=ipwCoeffs)
      do j = 1, nscat
         js = (j - 1)*lmax
         ipwCoeffsJ(js + 1:js + lmax) = &
            exp(cmplx(0, dot_product(kVec, xyz(:, j)), 8))*ipwCoeffs
      end do
      deallocate (ipwCoeffs)
      !
   end subroutine calcJCoeffsPW
   !
   subroutine calcCoeffsPW(ipwE0, ipwDirn, ipwCoeffs)
      !
      ! ============================================================
      ! Compute the coefficients for expressing an incident plane-
      ! wave in terms of vector spherical waves M_nm and N_nm, up
      ! to a maximum n_max. Follow equations C.57-59 on pg.377 of
      ! Mishchenko's book.
      !
      ! INPUT:
      ! ------
      ! ipwE0(3) - Plane-wave amplitude vector [COMPLEX]
      ! ipwDirn(3) - normalised direction vector [REAL]
      !
      ! OUTPUT:
      ! ------
      ! ipwCoeffs(2*pmax) - [COMPLEX] coefficients commensurate with
      !                   the abstract column vector defined in
      !                   eqn.B1 of Stout02.
      ! ============================================================
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables
      complex(8), dimension(3), intent(in) :: ipwE0
      real(8), dimension(3), intent(in) :: ipwDirn
      complex(8), dimension(:), intent(out) :: ipwCoeffs
      ! Local variables
      character(*), parameter :: myName = 'calcCoeffsPW'
      integer :: nmax, n, m, l, pmax
      real(8) :: rtp(3), cth, transform(3, 3)
      complex(8) :: prefM, prefN, V(3), E0(3)
      real(8), dimension(0:size(ipwCoeffs)/2) :: d, pi, tau
      complex(8), allocatable :: expo(:)
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      ! write(*,'(A,A,3(1x,"( ",e12.5E2,", ",e12.5E2," )"))') &
      !      myname,'> ipwE0(1:3)=',ipwE0(1:3)
      ! write(*,'(A,A,3(1x,e12.5E2))') &
      !      myname,'> ipwDirn(1:3)=',ipwDirn(1:3)
      !
      pmax = size(ipwCoeffs)
      if ((pmax/2)*2 == pmax) then
         pmax = pmax/2
      else
         write (*, '(A)') myName, '> ERROR: ipwCoeffs has odd length! Stopping..'
         STOP
      end if
      !
      ! Check the maximal value of the composite index
      call testPmax(myName, pmax, nmax)
      !
      allocate (expo(-nmax:nmax))
      !
      call xyz2rtp(ipwDirn, rtp, cth)
      ! Transform ipwE0 from cartesian to spherical polars
      call calcVTxyz2rtp(rtp, transform)
      E0(:) = matmul(transform(:, :), ipwE0(:))
      ! Compute the tau- and pi-functions
      call calcWignerd0andMore(cth, pmax, d, pi, tau)
      ! Compute the azimuth-angle exponentials
      do m = -nmax, nmax
         expo(m) = exp(-imu*m*rtp(3))
      end do
      !
      do n = 1, nmax
         prefM = sqrt(fpi*dble(2*n + 1)/dble(n*(n + 1)))*imu**n
         prefN = prefM/imu
         do m = -n, n
            call nm2p(n, m, l)
            ! Compute C_mn (using eqn C.19) in polars
            V(1:3) = (/cmplx(0, 0, 8), cmplx(0, pi(l), 8), cmplx(-tau(l), 0, 8)/)
            ipwCoeffs(l) = prefM*(-1)**m*expo(m)*dot_product(V(:), E0(:))
            ! Compute B_mn (using eqn C.18) in polars
            V(1:3) = (/cmplx(0, 0, 8), cmplx(tau(l), 0, 8), cmplx(0, pi(l), 8)/)
            ipwCoeffs(l + pmax) = prefN*(-1)**m*expo(m)*dot_product(V(:), E0(:))
         end do
      end do
      !
      deallocate (expo)
      !
      ! TEST>>>>>
      !write(*,*) 'ipwE0=', ipwE0
      !write(*,*) 'ipwCoeffs=', ipwCoeffs
      ! <<<<<TEST
      !
   end subroutine calcCoeffsPW
   !
   !==============================================================
   !
   subroutine offsetCoeffsPW(aJ, a, kVec, xyzr)
      !
      !============================================================
      ! Offset the coefficients (a) for a regular VSW expansion
      ! (centred at the origin) of an incident planewave to
      ! multiple centres specified xyzr, producing the scatterer-
      ! centred coefficients (aJ). kVec is the plane-wave-vector.
      !============================================================
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      complex(8), intent(in) :: a(:)
      real(8), intent(in) :: kVec(3), xyzr(:, :)
      complex(8), intent(out) :: aJ(size(xyzr, 2)*size(a)*(size(a) + 2))
      ! Local variables
      character(*), parameter :: myName = 'calcIncCoeffsJ'
      integer :: i, is, lmax, nmax, nscat
      complex(8) :: zdum !vtacs(size(a),size(a))
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      lmax = size(a)
      nmax = int(sqrt(real(lmax/2)))
      if (lmax /= 2*nmax*(nmax + 2)) then
         write (*, '(A,A)') myname, '> ERROR: lmax /= 2*nmax*(nmax+2)'
         STOP
      end if
      !
      if (size(xyzr, 1) /= 4) then
         write (*, '(A,A)') myname, '> ERROR: size(xyzr,1) /= 4'
         STOP
      end if
      nscat = size(xyzr, 2)
      !
      scat: do i = 1, nscat
         is = (i - 1)*lmax
         if (norm2(xyzr(1:3, i)) > 1.0d-10*xyzr(4, i)) then
            ! call calcOIJ(ij=(/i,0/), OIJ=vtacs, regt = .true.)
            ! aJ(is+1:is+lmax_G) = matmul(vtacs,a_(:))
            zdum = dot_product(kVec(1:3), xyzr(1:3, i)) ! assume planewave
            zdum = exp(cmplx(0, realpart(zdum), 8))
            aJ(is + 1:is + lmax) = zdum*a
         else ! vtacs should be near identity, but blow up...
            aJ(is + 1:is + lmax) = a
         end if
      end do scat
      !
   end subroutine offsetCoeffsPW
   !
   subroutine calcWignerBigD(angles, pmax, bigD)
      !
      ! ============================================================
      ! Compute the Wigner D-functions, D^s_{m,n}(alpha,beta,gamma),
      ! using equation B.38 on pg. 367 of Mishchenko's book (2002).
      !
      ! INPUT:
      ! -----
      ! angles(3) - (alpha,beta,gamma) in Radians, see eqn. (B.39).
      ! pmax - maximum composite multipole index l=s(s+1)+m
      !
      ! OUTPUT:
      ! ------
      ! D(pmax,pmax) - [COMPLEX] values for D^s_{m,n} in block-
      !                              diagonal matrix form.
      ! ============================================================
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables
      real(8), intent(in) :: angles(3)
      integer, intent(in) :: pmax
      complex(8), intent(out), dimension(pmax, pmax) :: bigD
      ! Local variables
      character(*), parameter :: myName = 'calcWignerBigD'
      integer :: nmax, n, m1, m2, l1, l2
      real(8) :: littled(0:pmax, 0:pmax)
      complex(8), allocatable :: expAlpha(:), expGamma(:)
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      call testPmax(myName, pmax, nmax)
      !
      bigD = 0
      !
      call calcWignerLittled(angles(2), pmax, littled)
      allocate (expAlpha(-nmax:nmax), expGamma(-nmax:nmax))
      do n = -nmax, nmax
         expAlpha(n) = exp(cmplx(0, -n*angles(1), 8))
         expGamma(n) = exp(cmplx(0, -n*angles(3), 8))
      end do
      do n = 1, nmax
         do m1 = -n, n
            call nm2p(n, m1, l1)
            do m2 = -n, n
               call nm2p(n, m2, l2)
               bigD(l1, l2) = expAlpha(m1)*littled(l1, l2)*expGamma(m2)
            end do
         end do
      end do
      deallocate (expAlpha, expGamma)
      !
   end subroutine calcWignerBigD
   !
   subroutine calcWignerLittled(theta, pmax, d)
      !
      ! ============================================================
      ! Compute the Wigner d-functions, d^s_{m,n}(theta), using
      ! the recurrence relation in equation B.22 on pg. 365 of
      ! Mishchenko's book (2002). The recurrence expects theta
      ! betweem 0 and Pi, so for theta < 0 the block-diagonal
      ! matrix d^s_{m,n} is transposed in the end, using Eq. (B.6).
      !
      ! INPUT:
      ! -----
      ! theta - angle in Radians
      ! pmax - maximum composite multipole index l=s(s+1)+m
      !
      ! OUTPUT:
      ! ------
      ! d(0:pmax,0:pmax) - [REAL] values for d^s_{m,n} in block-
      !                           diagonal matrix form.
      ! ============================================================
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables
      real(8), intent(in) :: theta
      integer, intent(in) :: pmax
      real(8), intent(out), dimension(0:pmax, 0:pmax) :: d
      ! Local variables
      character(*), parameter :: myName = 'calcWignerSmalld'
      integer :: sp1, s, m, l, ll, lm, llm, msq, smax, n, mmn, mpn, ipass
      real(8) :: x(2), ds, dsm1, xi, chi, dum0, dum1, dum2, y1(2), y2(2), pref
      real(8), allocatable :: sqrtfact(:), sqr(:), halftom
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      if (abs(theta) > pi) then
         write (*, *) myname, '> ERROR: |theta| > pi! Stopping...'
         STOP
      end if
      !
      call testPmax(myName, pmax, smax)
      !
      allocate (sqrtfact(0:2*smax), sqr(0:smax))
      do n = 0, smax
         sqr(n) = dble(n*n)
      end do
      x(1) = cos(theta) ! For 1st pass: m and n of same sign
      y1(1) = sqrt(1.0d0 - x(1))
      y2(1) = sqrt(1.0d0 + x(1))
      x(2) = cos(pi - theta) ! For 2nd pass: m and n of different sign
      y1(2) = sqrt(1.0d0 - x(2))
      y2(2) = sqrt(1.0d0 + x(2))
      !
      d(:, :) = 0.0d0 ! zero all the matrix elements
      !
      ! Outer-most loop over m [0,smax], second loop over n [0,m],
      ! and the third, inner-most loop performs the reccurrence on s
      ! (eqn. B.22 on pg. 365 of Mishchenko's book) with s_min = m
      ! in each case, initialised using equations B.23 and B.24.
      mloop: do m = 0, smax
         msq = m*m
         ! Accumulate all the required factorials on the fly
         if (m == 0) then
            sqrtfact(m) = 1.0d0
            halftom = 1.0d0
         else
            halftom = 0.5d0*halftom
            do n = 2*(m - 1) + 1, 2*m
               sqrtfact(n) = sqrt(dble(n))*sqrtfact(n - 1)
            end do
         end if
         nloop: do n = 0, m
            mmn = m - n
            mpn = m + n
            xi = dble((-1)**mmn)
            !pref=sqrtfact(2*m)/(sqrtfact(mmn)*sqrtfact(mpn)*dble(2**m))
            ! 2**m overflows for large m, so do 0.5d0**m instead...
            pref = halftom*sqrtfact(2*m)/(sqrtfact(mmn)*sqrtfact(mpn))
            ploop: do ipass = 1, 2 ! Two passes: for n*m >= 0 and for n*m < 0
               if (ipass == 2 .and. n == 0) exit ploop
               ! Follow eqn. B.23 to set d^{m-1}_{mn} = 0
               dsm1 = 0.0d0 ! d-function for (s-1)
               ! Use eqn. B.24 to set d^{m}_{mn}
               ds = xi*pref*y1(ipass)**mmn*y2(ipass)**mpn
               call nm2p(m, m, l)
               call nm2p(m, -m, lm)
               call nm2p(m, n, ll)
               call nm2p(m, -n, llm)
               ! Apply symmetry relations from equation B.5
               if (ipass == 1) then
                  d(l, ll) = ds ! d^{m}_{m,n}
                  if (m /= n) d(ll, l) = xi*ds ! d^{m}_{n,m}
                  if (m /= 0 .or. n /= 0) then
                     d(lm, llm) = xi*ds ! d^{m}_{-m,-n}
                     if (m /= n) d(llm, lm) = ds ! d^{m}_{-n,-m}
                  end if
               elseif (ipass == 2) then
                  chi = dble((-1)**(m - n))
                  d(lm, ll) = chi*ds ! d^{m}_{-m,n}
                  if (m /= -n) d(ll, lm) = chi*xi*ds ! d^{m}_{n,-m}
                  if (m /= 0 .or. n /= 0) then
                     d(l, llm) = chi*xi*ds  ! d^{m}_{m,-n}
                     if (m /= -n) d(llm, l) = chi*ds ! d^{m}_{-n,m}
                  end if
               else
                  write (*, '(A,A)') myname, '> ERROR: Bad pass index!'
                  STOP
               end if
               !
               sloop: do sp1 = m + 1, smax ! Reccurence for s > s_min = m
                  s = sp1 - 1
                  if (s == 0) then
                     dum0 = dble(2*s + 1)*x(ipass)
                     dum1 = sqrt(dble(sqr(s) - sqr(m)))
                     dum2 = sqrt(dble(sqr(sp1) - sqr(m)))
                  else
                     dum0 = dble(2*s + 1)*((sqr(s) + s)*x(ipass) - dble(m*n))
                     dum1 = dble(sp1)* &
                            sqrt((sqr(s) - sqr(m))*(sqr(s) - sqr(n)))
                     dum2 = dble(s)* &
                            sqrt((sqr(sp1) - sqr(m))*(sqr(sp1) - sqr(n)))
                  end if
                  call nm2p(sp1, m, l)
                  call nm2p(sp1, -m, lm)
                  call nm2p(sp1, n, ll)
                  call nm2p(sp1, -n, llm)
                  if (ipass == 1) then
                     d(l, ll) = (dum0*ds - dum1*dsm1)/dum2
                     dsm1 = ds
                     ds = d(l, ll) ! d^{s}_{m,n}
                     if (n /= m) d(ll, l) = xi*ds ! d^{s}_{n,m}
                     if (m /= 0 .or. n /= 0) then
                        d(lm, llm) = xi*ds ! d^{s}_{-m,-n}
                        if (n /= m) d(llm, lm) = ds ! d^{s}_{-n,-m}
                     end if
                  elseif (ipass == 2) then
                     chi = dble((-1)**(sp1 - n))
                     d(lm, ll) = (dum0*ds - dum1*dsm1)/dum2
                     dsm1 = ds
                     ds = d(lm, ll) ! d^{s}_{-m,n}
                     d(lm, ll) = chi*ds
                     if (n /= -m) d(ll, lm) = chi*xi*ds ! d^{s}_{n,-m}
                     if (m /= 0 .or. n /= 0) then
                        d(l, llm) = chi*xi*ds  ! d^{s}_{m,-n}
                        if (n /= -m) d(llm, l) = chi*ds ! d^{s}_{-n,m}
                     end if
                  else
                     write (*, '(A,A)') myname, '> ERROR: Bad pass index!'
                     STOP
                  end if
               end do sloop
            end do ploop
         end do nloop
      end do mloop
      !
      deallocate (sqrtfact, sqr)
      !
      if (theta < 0.0d0) d = transpose(d)
      !
   end subroutine calcWignerLittled
   !
   subroutine calcWignerd0andMore(x, pmax, d, pi, tau)
      !
      ! ============================================================
      ! Compute the Wigner d-functions (d^s_{m,n} with n=0) using
      ! the reccurence relation in equation B.22 on pg. 365 of
      ! Mishchenko's book (2002). Optionally, also compute the
      ! pi and tau functions defined on pg. 373. Note:
      !
      ! d^s_{m,0} = (-1)^m*d^s_{0,m} = (-1)^m*d^s_{-m,0}
      !           = sqrt{(s-m)!/(s+m)!}*P_s^m,
      !
      ! where P_s^{m} are associated Legendre functions.
      !
      ! INPUT:
      ! -----
      ! x - cos(theta)
      ! pmax - maximal value of the upper index
      !
      ! OUTPUT:
      ! ------
      ! d(0:pmax) - [REAL] values for d^s_{m,0}
      ! pi(0:pmax) - [OPTIONAL, REAL] values for pi_{m,s}
      ! tau(0:pmax) - [OPTIONAL, REAL] values for tau_{m,s}
      ! ============================================================
      !
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      ! Passed variables
      real(8), intent(in) :: x
      integer, intent(in) :: pmax
      real(8), intent(out), dimension(0:pmax) :: d
      real(8), intent(out), dimension(0:pmax), optional :: pi, tau
      ! Local variables
      character(*), parameter :: myname = 'calcWignerd0andMore'
      integer :: sp1, s, m, l, ll, msq, smax
      real(8) :: ds, dsm1, pis, pism1, xi, ysq, y, dum0, dum1, dum2
      real(8) :: oddfac, evenfac
      logical :: domore
      !---------------------------------------------------
      ! End of variable declarations. Directives start now
      !---------------------------------------------------
      !
      call testPmax(myName, pmax, smax)
      !
      if (present(pi) .and. present(tau)) then
         domore = .true.
      else
         domore = .false.
         if (present(pi) .or. present(tau)) then
            write (*, *) myname, '> ERROR: Deficient optional arguments! Stopping...'
            STOP
         end if
      end if
      !
      ! Outer loop over m, and the inner loop performs the reccurrence
      ! on s (eqn. B.22 on pg. 365 of Mishchenko's book) with s_min = m
      ! in each case, initialised using equations B.23 and B.24. Note
      ! for s_min = m: 2^{-m}/(m!)*sqrt{(2m)!} = sqrt{(2m-1)!!/(2m)!!},
      ! which starts from 1 and slooowly approaches zero as m increases.
      !
      oddfac = 0.0d0  ! initialise (2m-1)!!
      evenfac = 0.0d0 ! initialise (2m)!!
      ysq = 1.0d0 - x*x ! =1-cos^2(theta) = sin^2(theta)
      y = sqrt(ysq)
      mloop: do m = 0, smax
         oddfac = max(1.0d0, oddfac*(2*m - 1))
         evenfac = max(1.0d0, evenfac*(2*m))
         msq = m*m
         xi = dble((-1)**m)
         ! Follow eqn. B.23 to set d^{m-1}_{m0} = 0
         dsm1 = 0.0d0 ! d-function for (s-1)
         ! Use a special case of eqn. B.24 to set d^{m}_{m0}
         ds = xi*y**m*sqrt(oddfac/evenfac)
         call nm2p(m, m, l)
         call nm2p(m, -m, ll) ! NB: ll=l for m=0 !!!!!
         d(l) = ds
         d(ll) = xi*ds
         if (domore) then
            pism1 = 0.0d0 ! pi-function for (m-1)
            if (m == 0) then
               pis = 0.0d0 ! No division by zero for x=+/-1 when m=0
            else
               pis = dble(m)*y**(m - 1)*sqrt(oddfac/evenfac)
            end if
            pi(l) = pis
            pi(ll) = -xi*pis
            tau(ll) = x*xi*pis
            tau(l) = xi*tau(ll)
         end if
         sloop: do sp1 = m + 1, smax ! Apply reccurence for s > s_min = m
            s = sp1 - 1
            dum0 = dble(2*s + 1)*x
            dum1 = sqrt(dble(s**2 - msq))
            dum2 = sqrt(dble(sp1**2 - msq))
            call nm2p(sp1, m, l)
            call nm2p(sp1, -m, ll) ! NB: ll=l for m=0
            d(l) = (dum0*ds - dum1*dsm1)/dum2
            dsm1 = ds
            ds = d(l)
            ! Use d^{s}_{m,0} = (-1)^m*d^{s}_{-m,0} for m < 0
            d(ll) = xi*ds
            if (domore) then
               if (m > 0) then ! Do the same recursion for pi
                  pi(l) = (dum0*pis - dum1*pism1)/dum2
                  ! Use equation B.14 in Walter's thesis to build
                  ! up tau_{s,m} for m > 0
                  tau(l) = (dble(sp1)*x*pi(l) - dum2*pis)/dble(m)
                  tau(ll) = xi*tau(l)
               else
                  pi(l) = 0.0d0
               end if
               pism1 = pis
               pis = pi(l)
               pi(ll) = -xi*pis
            end if
         end do sloop
      end do mloop
      !
      ! Fill tau_{s,m} for m=0 and s > 0 using a different recurrence,
      ! relating Dd^n_{m,0} (D indicates differentiation) to
      ! d^n_{m+1,0} and d^n_{m-1,0} (see Wikipedia), which simplifies
      ! for m=0 and does not involve division by m or sin(theta).
      if (domore) then
         do s = 1, smax
            call nm2p(s, 0, l)
            tau(l) = sqrt(dble(s*(s + 1)))*d(l + 1)
         end do
      end if
      !
   end subroutine calcWignerd0andMore
   !
   subroutine calcRiccatiBessels(z, nmax, regt, f, df)
      !

      ! ============================================================
      ! Compute the Riccati-Bessel functions psi_n (if regt=TRUE) or
      ! xi_n (regt=FRALSE), and their derivatives, for n=1,...,nmax,
      ! using calcSphBessels.
      !
      ! INPUT:
      ! -----
      ! z - scalar argument [COMPLEX]
      ! nmax - maximal order of the functions [INTEGER]
      ! regt - psi (regular) or xi (irregular)? [LOGICAL]
      !
      ! OUTPUT:
      ! ------
      ! f(1:nmax) - Riccati-Bessel functions psi_n(z)=z*j_n(z) or
      !             xi_n(z)=z*h_n(z) for n=1,...,nmax.
      ! df(1:nmax)- The corresponding derivatives w.r.t. z.
      !
      ! NOTE: Will stop the program if ZBESH/ZBESJ gives an error!
      ! ============================================================
      !
      ! Passed variables
      complex(8), intent(in) :: z
      integer, intent(in) :: nmax
      logical, intent(in) :: regt
      complex(8), intent(out) :: f(nmax), df(nmax)
      ! Local variables
      integer :: n
      complex(8) :: funs(0:nmax)
      !
      ! Note: df_{n} = f_{n-1} - n*f_{n}/z = z*g_{n-1} - n*g_{n},
      ! where f_n = z*g_n
      call calcSphBessels(z, nmax, regt, funs)
      do n = nmax, 1, -1
         df(n) = z*funs(n - 1) - n*funs(n)
         f(n) = z*funs(n)
      end do
      !
   end subroutine calcRiccatiBessels
   !
   subroutine calcSphBessels(z, nmax, regt, bes)
      !
      ! ============================================================
      ! A wrapper routine for computing spherical Bessel/Hankel
      ! functions of the first kind for a complex argument z. Use
      ! the standard Amos' TOMS644 routines ZBESJ (if regt=TRUE) or
      ! ZBESH (if regt=FALSE) in the 'toms644.f' file to compute
      ! the usual (non-spherical) Bessel/Hankel functions of
      ! half-integer order n+1/2 for n=0,...,nmax.
      !
      ! INPUT:
      ! -----
      ! z - scalar argument [COMPLEX]
      ! nmax - maximal order of the functions [INTEGER]
      ! regt - Bessel (regular) or Hankel (irregular)? [LOGICAL]
      !
      ! OUTPUT:
      ! ------
      ! bes(0:nmax) - Bessel (J_{n+1/2}) or Hankel (H_{n+1/2})
      !               functions (of 1st kind) for n=0,...,nmax
      !
      ! NOTE: Will stop the program if ZBESH/ZBESJ gives an error!
      ! ============================================================
      !
      ! Passed variables

      complex(8), intent(in) :: z
      integer, intent(in) :: nmax
      logical, intent(in) :: regt
      complex(8), intent(out) :: bes(0:nmax)
      ! Local variables
      character(*), parameter :: myname = 'calcSphBessels'
      integer :: n, ierr
      real(8), dimension(nmax + 1) :: bes_r, bes_i
      character(len=5) :: besRoutine
      complex(8) :: factor
      !
      !write(*,*)'z=', z
      if (regt) then
         call ZBESJ(realpart(z), imagpart(z), 0.5d0, 1, nmax + 1, &
                    bes_r, bes_i, n, ierr)
         besRoutine = 'ZBESJ'
      else
         call ZBESH(realpart(z), imagpart(z), 0.5d0, 1, 1, nmax + 1, &
                    bes_r, bes_i, n, ierr)
         besRoutine = 'ZBESH'
      end if
      if (n /= 0) then
         write (*, '(A,A,A,A,i6)') &
            myname, '> WARNING: ', besRoutine, &
            ' set NZ components to zero. NZ= ', n
      end if
      if (ierr /= 0) then
         write (*, '(A,A,A,A,i2,A,2(1x,es15.8E2))') &
            myname, '> ', besRoutine, &
            ' gave ERROR flag ', ierr, ' for z=', z
         do n = 0, nmax
            write (*, '(f6.1,2(1x,es16.8E3))') &
               real(n) + 0.5, bes_r(n + 1), bes_i(n + 1)
         end do
         STOP
      end if
      !
      factor = sqrt(hpi/z)
      do n = 0, nmax
         bes(n) = factor*cmplx(bes_r(n + 1), bes_i(n + 1), 8)
      end do
      !
      return
      !
   end subroutine calcSphBessels
   !
   subroutine xyz2rtp(xyz, rtp, cth)
      !
      ! ============================================================
      ! Transform the cartesian coordinates (x,y,z) of a point in
      ! 3D space to spherical polar coordinates (r,theta,phi) using
      ! the ISO convention, i.e. for physics: radius r [0,infinity),
      ! inclination theta [0,pi], azimuth phi [0,2*pi]. Formulae:
      !
      !    r     = |xyz|
      !    theta = arccos(z/r)
      !    phi   = arctan(y/x)
      !
      ! INPUT:
      ! -----
      ! xyz(3) - cartesian coordinates [REAL]
      !
      ! OUTPUT:
      ! ------
      ! rtp(3) - spherical polar coordinates [REAL]
      ! cth    - cos(theta) [REAL, OPTIONAL]
      ! ============================================================
      !
      real(8), intent(in) :: xyz(3)  ! (x,y,z)
      real(8), intent(out) :: rtp(3) ! (r,theta,phi)
      real(8), intent(out), optional :: cth ! cos(theta)
      !
      real(8) :: costheta
      !
      rtp(1) = sqrt(dot_product(xyz, xyz)) ! r, non-negative
      if (rtp(1) < 1.0d-15) then
         rtp = 0
         if (present(cth)) cth = 1
         return
      end if
      costheta = xyz(3)/rtp(1)
      if (costheta /= costheta) then ! NaN due to (x,y,z) = (0,0,0)
         rtp(1:3) = 0.0d0
         if (present(cth)) cth = 1.0d0
         return
      end if
      rtp(2) = acos(costheta)       ! theta: [0,Pi]
      rtp(3) = atan2(xyz(2), xyz(1)) ! phi: (-Pi,Pi] or (0,2*Pi]?
      ! Alternative: rtp(3) = acos(xyz(1)/(rtp(1)*sin(rtp(2))))
      if (rtp(3) /= rtp(3)) then ! NaN due to (x,y) = (0,0)
         rtp(3) = 0.0d0 ! Is this sensible?
      end if
      if (present(cth)) cth = costheta
      !
   end subroutine xyz2rtp
   !
   ! subroutine rtp2xyz(rtp,xyz)
   !
   ! ============================================================
   ! The inverse of xyz2rtp: transform the spherical polar
   ! coordinates (r,theta,phi) of a point in 3D space to
   ! cartesian coordinates (x,y,z). Formulae:
   !
   !    x = r*sin(theta)*cos(phi)
   !    y = r*sin(theta)*sin(phi)
   !    z = r*cos(theta)
   !
   ! INPUT:
   ! ------
   ! rtp(3) - spherical polar coordinates [REAL]
   !
   ! OUTPUT:
   ! ------
   ! xyz(3) - cartesian coordinates [REAL]
   ! ============================================================
   !
   ! real(8), intent(in) :: rtp(3)
   !real(8), intent(out) :: xyz(3)
   !
   ! xyz(1) = rtp(1)*sin(rtp(2))
   ! xyz(2) = xyz(1)*sin(rtp(3))
   ! xyz(1) = xyz(1)*cos(rtp(3))
   ! xyz(3) = rtp(1)*cos(rtp(2))
   !
   ! return
   !
   !end subroutine rtp2xyz
   !
   subroutine calcVTrtp2xyz(rtp, transform)
      !
      ! Compute matrix T: (E_r, E_t, E_p) -> (E_x, E_y, E_z),
      ! i.e. E(xyz) = T.E(rtp), at point rtp (in polars)
      !
      real(8), intent(in) :: rtp(3)
      real(8), intent(out) :: transform(3, 3)
      !
      transform(1, 1:3) = &
         (/sin(rtp(2))*cos(rtp(3)), &
           cos(rtp(2))*cos(rtp(3)), &
           -sin(rtp(3))/)
      transform(2, 1:3) = &
         (/sin(rtp(2))*sin(rtp(3)), &
           cos(rtp(2))*sin(rtp(3)), &
           cos(rtp(3))/)
      transform(3, 1:3) = &
         (/cos(rtp(2)), -sin(rtp(2)), 0.0d0/)
      !
      return
      !
   end subroutine calcVTrtp2xyz
   !
   subroutine calcVTxyz2rtp(rtp, transform)
      !
      ! Compute matrix T: (E_X, E_Y, E_Z) -> (E_R, E_T, E_P),
      ! i.e. E(RTP) = T.E(XYZ), at point rtp (in polars).
      !
      real(8), intent(in) :: rtp(3)
      real(8), intent(out) :: transform(3, 3)
      !
      transform(1, 1:3) = &
         (/sin(rtp(2))*cos(rtp(3)), &
           sin(rtp(2))*sin(rtp(3)), &
           cos(rtp(2))/)
      transform(2, 1:3) = &
         (/cos(rtp(2))*cos(rtp(3)), &
           cos(rtp(2))*sin(rtp(3)), &
           -sin(rtp(2))/)
      transform(3, 1:3) = &
         (/-sin(rtp(3)), cos(rtp(3)), 0.0d0/)
      !
      return
      !
   end subroutine calcVTxyz2rtp
   !
   ! ==============================================================
   ! Vector spherical harmonics are spanned by two indices: n and m,r
   ! define a generalised index l:=n(n+1)+m, for which each integral
   ! value of l corresponds to a unique (n,m).
   subroutine nm2p(n, m, l)
      ! Compute composite index l from a given n >= |m| >= 0
      integer, intent(in) :: n, m
      integer, intent(out) :: l
      !
      l = n*(n + 1) + m
      !
      return
      !
   end subroutine nm2p
   !
   subroutine p2nm(p, n, m)
      ! Compute unique (n,m) from a given composite index p
      integer, intent(in) :: p
      integer, intent(out) :: n, m
      !
      n = int(sqrt(real(p))) ! multipole index n
      m = p - n*(n + 1) ! angular index m
      !
   end subroutine p2nm
   !
   ! Some recurrences are defined only for m >= 0, in which
   ! case we shall use a second version of the composite
   ! index p_v2:=n*(n+1)/2+m. Still not sure how to compute
   ! the inverse.. i.e. p_v2 -> (n,m), but it's not needed.
   subroutine nm2pv2(n, m, p)
      ! Compute composite index p_v2 from a given n >= m >= 0
      integer, intent(in) :: n, m
      integer, intent(out) :: p
      !
      p = n*(n + 1)/2 + m
      !
   end subroutine nm2pv2
   !
   subroutine testPmax(name, pmax, nmax)
      ! Test pmax for commensurability, i.e. is pmax=nnmax*(nmax+2)
      ! and nmax=mmax? If not, then stop program. Will be phased out,
      ! eventually.
      character(*), intent(in) :: name
      integer, intent(in) :: pmax
      integer, intent(out) :: nmax
      integer :: m
      !
      call p2nm(pmax, nmax, m)
      if (nmax /= m) then
         write (*, '(A,A)') &
            name, '> ERROR: Incommensurate pmax! Stopping...'
         STOP
      end if
      !
   end subroutine testPmax
   !=====================================================
   ! In this subroutine the absorption matrix, gamma_j is calculated based on
   ! Eq. 49, 50 Stout2002, Note: since we are working with non magnetic medium, mu_j and mu=1

   subroutine calcAbsMat(Xi, ro, mat)

      complex(8), intent(in) :: Xi, ro
      complex(8), intent(inout) :: mat(:, :)
      integer :: lmax, pmax, nmax, n, m, p
      complex(8) :: z
      complex(8), allocatable :: f(:), df(:), f2(:), df2(:)
!---------------------------------
      lmax = size(mat, 1)
      pmax = lmax/2
      nmax = int(sqrt(dble(pmax)))
      allocate (f(nmax), df(nmax), f2(nmax), df2(nmax))
      z = Xi
      call calcRiccatiBessels(z, nmax, .true., f, df)
      z = ro*Xi
      call calcRiccatiBessels(z, nmax, .true., f2, df2)

      do n = 1, nmax
         do m = -n, n
            p = n*(n + 1) + m
            mat(p, p) = realpart(imu*ro*conjg(f2(n))*df2(n))/(abs(f2(n)*df(n) - ro*df2(n)*f(n)))**2  !Cnj
            mat(p + pmax, p + pmax) = realpart(imu*conjg(ro)*conjg(f2(n))*df2(n)) &
                                      /(abs(ro*f2(n)*df(n) - df2(n)*f(n)))**2  !Dnj
         end do
      end do

   end subroutine calcAbsMat
   !========================================================
   ! In this subroutine the Lambda matrix, is calculated based on
   ! Eq. 53 Stout2002, Note: since we are working with non magnetic medium, mu_j and mu=1

   subroutine calcLamMat(Xi, ro, mat)

      complex(8), intent(in) :: Xi, ro
      complex(8), intent(inout) :: mat(:, :)
      integer :: lmax, pmax, nmax, n, m, p
      complex(8) :: z
      complex(8), allocatable :: f(:), df(:), f2(:), df2(:)
!---------------------------------
      lmax = size(mat, 1)
      pmax = lmax/2
      nmax = int(sqrt(dble(pmax)))
      allocate (f(nmax), df(nmax), f2(nmax), df2(nmax))
      z = Xi
      call calcRiccatiBessels(z, nmax, .true., f, df)
      z = ro*Xi
      call calcRiccatiBessels(z, nmax, .true., f2, df2)

      do n = 1, nmax
         do m = -n, n
            p = n*(n + 1) + m
            mat(p, p) = (imu*ro)/(ro*f(n)*df2(n) - f2(n)*df(n))  !Lambda_Mj
            mat(p + pmax, p + pmax) = (imu*ro) &
                                      /(f(n)*df2(n) - ro*f2(n)*df(n))  !Lambda_Nj

         end do
      end do

   end subroutine calcLamMat
end module swav
