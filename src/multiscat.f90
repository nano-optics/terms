module multiscat
    !
    !=================================================================
    ! This module contains routines for solving a multiple scattering
    ! problem using the T-matrix formalism. The routines are split
    ! into four sets: (1) high-level, (2) core, (3) low-level, and
    ! (4) supplementary. Low-level routines should not be accessed
    ! from outside and are therefore made private.
    !
    !------------------------------
    ! List of subroutines/functions:
    !------------------------------
    ! 3: calcStokesScaVec> compute the Stokes phase matrix and
    !                      scattering Stokes vector
    !                |--------------------------------->>>[high-level]
    !
    ! 2: spectrumFF> Compute spectra for (many) fixed-orientation
    !                 and/or orientationally averaged cross-sections
    !                 for a system of multiple scatterers. The 1-body
    !                 T-matrices for individual scatterers are either
    !                 constructed using Mie theory or read from an
    !                 optional argument Tmat.
    !                |--------------------------------->>>[high-level]
    !
    ! 1: mapNF> Compute the scattered field coefficients and map the
    !            local enhancement factor over a set of spacial grid
    !            points for different incidece angles.
    !           |-------------------------------------->>>[high-level]
    !
    ! 0: solve> Solve a multi-scat problem for a given wavelength and
    !           compute various quantities.
    !          |--------------------------------------------->>>[coOre]
    !lac
    ! -1: stageAmat>
    !
    ! -2: calcTIJStout>
    !
    ! -3: calcTIMackowski>
    !
    ! -4: balanceMatJI>
    !
    ! -5: balanceVecJ>
    !
    ! -6: calcCrossSections>
    !
    ! -7: calcOAprops>
    !ical
    ! -8: applyRotTranzRotOnMat>
    !
    !-12: calcField>
    !
    !-13: dumpMatrix>
    !
    !-14: offsetTmat>
    !-15: calcStokesScaVec>  Calculates the Stokes Scattering vector and Stokes phase matrix
    !                        based on the formula 2.102 of Mishchenko's book on page 51.
    !
    !-16: calcLDOC> Calculates the optical chirality and normalized optical
    !                  chirality relative to the oc of RHCP light.
    !
    !-17:calcOaExtField> Calculates the orientation average of the total external electric
    ! field intensity based on formulae 16 of Stout 2008
    !
    !-18: calcOaLDOC> Calculates orientation-averaged local degree of optical chirality
    !----------------------------------------------------------------
    ! =================================================================
    !
    implicit none
    !
    private
    !
    public :: mapNF, spectrumFF, solve, &
              dumpMatrix, calcOAprops, cs_stout_G, &
              dumpTmatCol_G, dumpStagedA_G, dumpPrestagedA_G, &
              dumpScaCoeff_G, dumpIncCoeff_G, &
              isolve_G, balScale_G, balStout_G, &
              transInv_G, tfilename_G, calcStokesScaVec
    !
    !---------------------------------------------------
    ! Start of global variables (_G) declaration
    !---------------------------------------------------
    !
    integer :: verb_G = 1
    !
    real(8), parameter :: tiny1 = 1.0d-8, tiny = 1.0d-12, tiny2 = 1.0d-16
    real(8), parameter :: pi = acos(-1.0d0)
    real(8), parameter :: tpi = 2*pi
    !
    real(8), parameter :: eps0 = 8.8541878128d-12
    real(8), parameter :: mu0 = 4.0*pi*1.0d-7
    integer, parameter :: tunit = 987
    complex(8), parameter :: imu = (0.0d0, 1.0d0) ! the imaginary unit i
    real(8) :: hostK_G    ! Wavenumber for the medium
    real(8), parameter :: sp_light = 2.99792458d8
    !
    complex(8), allocatable :: bes4bal_G(:, :, :) ! Weights for matrix balancing
    !
    ! Some logical triggers set from termsProgram, invisible to PyTERMS
    logical :: dumpTmatCol_G = .false.
    logical :: dumpStagedA_G = .false.
    logical :: dumpPrestagedA_G = .false.
    logical :: dumpScaCoeff_G = .false.
    logical :: dumpIncCoeff_G = .false.
    logical :: balStout_G = .true.
    logical :: transInv_G = .false.
    logical :: cs_stout_G = .false.
    integer :: isolve_G = 0
    real(8) :: balScale_G = 1.0d0
    real(8) :: rtol_G = tiny1
    !
    character(len=64) :: tfilename_G = 'tmat_col.txt', afilename_G = 'alpha_col.txt'
    !
    !---------------------------------------------------
    ! End of global variables declaration
    !---------------------------------------------------
    !
contains

!-------------------------------------------------------------------------------------------------------------------
    subroutine mapNF( &
        ncut, wavelen, inc, ehost, geometry, scheme, &      !Mandatory ins
        field, Bfield, N_OC, orAvextEB_int, oa_ldoc, p_label, & ! in/out
        tfiles_, escat_, nselect_, verb_, noRTR_, dump_oaE2, &
        dump_oaB2, dipoles, HDF5_in) ! optional ins

        ! Map the near field for a multiple scattering problem, after solving it.
        ! Currently restricted to either (i) all scatterers being described using
        ! Mie theory, in which case scatMat(nscat,3); or (ii) all scatterers
        ! described by a single T-matrix, in which case have scatMat(lmax,lmax)
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        integer, intent(in) :: scheme, ncut(3)
        real(8), intent(in) :: wavelen(:), ehost(size(wavelen)), inc(:, :)
        real(8), intent(in) :: geometry(:, :)
        real(8), intent(in), optional :: dipoles(:, :)
        real(8), intent(inout) :: field(:, :, :, :, :)  ! 3 x 5 x npts x nwavelen x nInc
        integer, intent(inout) ::  p_label(:, :)
        real(8), intent(inout), optional :: Bfield(:, :, :, :, :)
        real(8), intent(inout), optional ::  N_OC(:, :, :, :, :)
        real(8), intent(inout), optional :: orAvextEB_int(:, :, :)  !npts x 2x(total,B0,C0) x nwavelen
        real(8), intent(inout), optional :: oa_ldoc(:, :, :)  !npts x (total,A0,B0,D0) x nwavelen
        complex(8), intent(in), optional :: &
            escat_(size(geometry, 2), 4, size(wavelen))
        integer, intent(in), optional :: nselect_(2, 2, 2, size(geometry, 2)), verb_
        logical, intent(in), optional :: noRTR_
        logical, intent(in) :: dump_oaE2, dump_oaB2
        character*256, intent(in), optional :: tfiles_(size(geometry, 2))
	 logical, intent(in), optional :: HDF5_in(:)
        !real(8), intent(inout), optional :: acs_int_(size(geometry,2), 4)
        ! Local variables
        character(*), parameter :: myName = 'mapNF'
        character(256) :: tfiles(size(geometry, 2))
        logical :: scatMiet(size(geometry, 2)), noRTR, multisplit
        logical :: oravE = .false.
        integer :: nmax, lmax, nscat, nwav, ierr, j, js, verb, i, is, iwav, nfi = 0
        integer :: nselect(2, 2, 2, size(geometry, 2)), units(size(geometry, 2)), nd = 0
        complex(8), allocatable ::TIJ(:, :), cJ(:, :, :), cJint(:, :, :)
        !real(8), allocatable :: inc(:,:)
        real(8) :: ipwDirn(3, size(inc, 2))
        complex(8) :: ipwE0(3, size(inc, 2))
        integer :: ngrid
        logical ::  PWinc = .false.
        integer :: tic, toc, tps
        real    :: t0, t1
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        nscat = size(geometry, 2)
        nmax = ncut(1)
        rtol_G = 10.0d0**ncut(3)
        lmax = 2*nmax*(nmax + 2)
        nwav = size(wavelen)
        ngrid = size(field, 3)
        !---------------------------------------------------
        ! Start processing arguments and initialisation
        !---------------------------------------------------
        if (present(dipoles)) then
        else
            PWinc = .true.
        end if
        if (scheme /= 0 .AND. PWinc) then
            if (dump_oaE2 .OR. dump_oaB2) oravE = .true.

        else

            write (*, *) myName, '> WARNING: orientation averaged quantities are not calculated at scheme = 0'
            !endif
        end if

        !else
        !write(*,*) myName,'> WARNING: orAvextE_int is not calculated at scheme=0'
        !endif
        if (size(geometry, 1) /= 8) then
            write (*, *) myName, '> ERROR: size(geometry,1) /= 8'
            STOP
        end if

        if (present(tfiles_)) then

            do j = 1, nscat
                if (len_trim(tfiles_(j)) > 0) then
                    tfiles(j) = tfiles_(j)
                    units(j) = 1000 + j
                    iloop: do i = 1, j - 1
                        ! Check for filename duplicates and, if present,
                        ! store scatterer ID with the same T-matrix file.
                        ! Store -ve of the ID to distinguish from actual
                        ! file handles in "units".
                        if (tfiles(j) == tfiles(i)) then
                            units(j) = -i
                            exit iloop
                        end if
                    end do iloop
                else
                    tfiles(j) = ''
                    units(j) = 0
                end if
            end do
        else
            tfiles = ''
            units = 0
        end if
        if (maxval(geometry(8, :)) >= tiny1) then
            ! non-Mie scatterer(s) present
            if (maxval(units) <= 0) then
                write (*, *) myName, &
                    '> ERROR: Missing tfiles for non-Mie scatterers'
                STOP
            end if
        end if
        !
        if (minval(geometry(8, :)) < tiny1) then
            ! Mie scatterer(s) present
            noRTR = .false. ! dummy
            if (present(escat_)) then
                if (maxval(abs(escat_)) > tiny2) then
                    noRTR = .true. ! non-zero escat_ present
                end if
            end if
            if (.not. noRTR) then
                write (*, *) myName, '> ERROR: Missing non-zero escat_'
                STOP
            end if
        end if
        !
        if (scheme < 0) then
            write (*, *) myName, '> ERROR: scheme < 0'
            STOP
        end if
        ! Make optional arguments compatible with f2py's default
        nselect = 0
        if (present(nselect_)) nselect = nselect_
        verb = 0
        if (present(verb_)) verb = verb_

        if (verb > 2) then
            write (*, '(A,A, 5x,i8)') myName, '> Number of grid-points:', ngrid
            write (*, '(12x,A,24x,A,24x,A)') 'x', 'y', 'z'
            do j = 1, ngrid
                write (*, *) field(1:3, 1, j, 1, 1)
            end do
        end if
        noRTR = .false.
        if (present(noRTR_)) noRTR = noRTR_
        !
        if (PWinc) then  !PWinc
            multisplit = .false.
            nfi = size(inc, 2)
            allocate (cJ(nscat*lmax, 2, nfi))  !, inc(4,nfi)

        elseif (present(dipoles)) then
            nd = size(inc, 2)  !number of dipoles
            allocate (cJ(nscat*lmax, 2, nd))
        end if

        allocate (TIJ(nscat*lmax, nscat*lmax))
        allocate (cJint(nscat*lmax, 4, 2))
        cJint = 0
        cJ = 0
        cJ(1:4, 1, :) = inc(1:4, :)

        do iwav = 1, nwav
            !
            if (verb > 0) then
                write (*, *)
                write (*, '(A,A,f8.2,A)') myname, '> ===== Wavelength: ', &
                    wavelen(iwav), ' (nm) ======================'
                !write(*,*)
            end if
            TIJ = 0
            do j = 1, nscat
                js = (j - 1)*lmax
                if (geometry(8, j) < tiny1) then ! Mie scatterer
                    scatMiet(j) = .true.
                    do i = 1, 1 - nint(geometry(8, j))
                        TIJ(js + i, js + i) = escat_(j, i, iwav)
                    end do
                else
                    scatMiet(j) = .false.
                    if (units(j) > 0) then
                        call readTmatFile( &
                            Tmat=TIJ(js + 1:js + lmax, js + 1:js + lmax), &
                            filename=trim(tfiles(j)), &
                            unit=units(j), &
                            wavelen=wavelen(iwav), &
                            HDF5_in_ = HDF5_in(j))
                    elseif (units(j) < 0) then
                        i = abs(units(j))
                        is = (i - 1)*lmax
                        TIJ(js + 1:js + lmax, js + 1:js + lmax) = &
                            TIJ(is + 1:is + lmax, is + 1:is + lmax)
                    else
                        write (*, '(A,A,i5)') &
                            myName, '> ERROR: units(j)=0 for j=', j
                        STOP
                    end if
                end if
            end do
            !
            if (verb > 0) write (*, *)
            if ((oravE .OR. present(oa_ldoc)) .AND. scheme == 3) &  !! there is currently a mismatch in the results
                write (*, *) myName, '> WARNING: orientation-averaged near-field quantities are inaccurate at Scheme = 3'

            cJ = 0
            cJ(1:4, 1, :) = inc(1:4, :)
            cJint = 0
            call solve( &
                wavelen=wavelen(iwav), &
                ehost=ehost(iwav), &
                geometry=geometry, &
                TIJ=TIJ, &
                scheme_=scheme, &
                cJ_=cJ, &
                cJint_=cJint, &
                !csAbs_ = acs_int, &
                noRTR_=noRTR, &
                nselect_=nselect, &
                verb_=verb, &
                ierr_=ierr, &
                dipoles=dipoles)
            ! endif

            if (oravE) then
                if (verb > 1) then
                    call cpu_time(t0)
                    call system_clock(tic)
                end if
                call calcOaExtField( &
                    r=field(:, 1, :, 1, 1), &
                    geometry=geometry, &
                    TIJ=TIJ, &
                    oaEB2=orAvextEB_int(:, :, iwav), &
                    lambda=wavelen(iwav), &
                    ehost=ehost(iwav), &
                    p_label=p_label, &
                    scatK_=tpi/wavelen(iwav)*sqrt(escat_(:, :, iwav)), &
                    verb_=verb)

                if (verb > 1) then
                    call cpu_time(t1)
                    call system_clock(toc, count_rate=tps)
                    write (*, '(A,A,2(1x,es10.3E2))') &
                        myName, '> Calculation time [calcOaExtField] (CPU & real in s): ', &
                        t1 - t0, dble(toc - tic)/dble(tps)
                    write (*, *)
                end if

            end if
            if (present(oa_ldoc)) then
                if (verb > 1) then
                    call cpu_time(t0)
                    call system_clock(tic)
                end if
                call calcOaLDOC(pol_type=nint(inc(1, 1)), &
                                r=field(:, 1, :, 1, 1), &
                                geometry=geometry, &
                                TIJ=TIJ, &
                                Or_OC=oa_ldoc(:, :, iwav), &
                                lambda=wavelen(iwav), &
                                ehost=ehost(iwav), &
                                p_label=p_label, &
                                scatK_=tpi/wavelen(iwav)*sqrt(escat_(:, :, iwav)), &
                                verb_=verb)
                !
                if (verb > 1) then
                    call cpu_time(t1)
                    call system_clock(toc, count_rate=tps)
                    write (*, '(A,A,2(1x,es10.3E2))') &
                        myName, '> Calculation time [calcOaLDOC] (CPU & real in s): ', &
                        t1 - t0, dble(toc - tic)/dble(tps)
                    write (*, *)
                end if
            end if

            if (verb > 1) then
                call cpu_time(t0)
                call system_clock(tic)
            end if

            inc_loop: do i = 1, nfi

                if (size(field, 1) == 3 .and. size(field, 2) == 5 .and. size(field, 3) > 0) then

                    call parseInc( &
                        inc(:, i), &
                        inc_dirn=ipwDirn(:, i), &
                        inc_ampl=ipwE0(:, i), &
                        verb_=verb)

                    if (present(escat_)) then

                        call calcField( &  !this is the case we have core-shell
                            r=field(:, 1, :, 1, 1), &                       !
                            geometry=geometry, &
                            ipwVec=tpi/wavelen(iwav)*sqrt(ehost(iwav))*ipwDirn(1:3, i), &
                            ipwE0=ipwE0(1:3, i), &
                            scaCJ=cJ(:, 2, i), &
                            intCJreg_=cJint(:, :, 1), intCJirr_=cJint(:, :, 2), &
                            reE=field(:, 2, :, iwav, i), imE=field(:, 3, :, iwav, i), &
                            reB=Bfield(:, 2, :, iwav, i), imB=Bfield(:, 3, :, iwav, i), &
                            scatK_=tpi/wavelen(iwav)*sqrt(escat_(:, :, iwav)), &
                            reE_sca=field(:, 4, :, iwav, i), imE_sca=field(:, 5, :, iwav, i), &
                            reB_sca=Bfield(:, 4, :, iwav, i), imB_sca=Bfield(:, 5, :, iwav, i), &
                            p_label=p_label, &
                            dipoles=dipoles, &
                            verb_=verb)

                    else

                        call calcField( &
                            r=field(:, 1, :, 1, 1), &
                            geometry=geometry, &
                            ipwVec=tpi/wavelen(iwav)*sqrt(ehost(iwav))*ipwDirn(1:3, i), &
                            ipwE0=ipwE0(1:3, i), &
                            scaCJ=cJ(:, 2, i), &
                            reE=field(:, 2, :, iwav, i), imE=field(:, 3, :, iwav, i), &
                            reB=Bfield(:, 2, :, iwav, i), imB=Bfield(:, 3, :, iwav, i), &
                            reE_sca=field(:, 4, :, iwav, i), imE_sca=field(:, 5, :, iwav, i), &
                            reB_sca=Bfield(:, 4, :, iwav, i), imB_sca=Bfield(:, 5, :, iwav, i), &
                            p_label=p_label, &
                            dipoles=dipoles, &
                            verb_=verb)

                    end if
                    if (present(N_OC)) then
                        call calcLDOC(Ef=field(:, 1:3, :, iwav, i), &
                                      Bf=Bfield(:, 1:3, :, iwav, i), &
                                      N_OpC=N_OC(:, :, :, iwav, i), &
                                      verb_=verb)
                    end if
                    Bfield(:, 2:5, :, iwav, i) = (wavelen(iwav)/(tpi*sp_light))*Bfield(:, 2:5, :, iwav, i)
                    !
                elseif (verb > 0) then
                    !
                    write (*, '(A,A)') myName, '> No field evaluation for misshapen array.'
                    !
                end if

            end do inc_loop
            if (verb > 1) then
                call cpu_time(t1)
                call system_clock(toc, count_rate=tps)
                write (*, '(A,A,2(1x,es10.3E2))') &
                    myName, '> Calculation time [NFs and LDOC for all inc.] (CPU & real in s): ', &
                    t1 - t0, dble(toc - tic)/dble(tps)
                write (*, *)
            end if
            if (verb > 2 .AND. present(N_OC)) then
                write (*, *) myName, '> Normalized LDOC for diff inc.'
                do j = 1, ngrid
                    write (*, '(1000(1x, ES24.17E2))') N_OC(1, 1, j, iwav, :)
                end do
                write (*, *)
            end if

            if (verb > 0) write (*, '(A,A,/)') myname, '> Done!'
        end do !iwav

        if (allocated(TIJ)) deallocate (TIJ)
        if (allocated(cJ)) deallocate (cJ)
        if (allocated(cJint)) deallocate (cJint)
!if(allocated(inc)) deallocate(inc)

    end subroutine mapNF

    !----------------------------------------------------------------------------

    !=================================================================
    !
    subroutine spectrumFF( &
        ncut, wavelen, ehost, geometry, scheme, &
        sig_oa_, sig_, sig_abs_, jsig_abs_oa, & ! output (cross-sectine outions)
        escat_, tfiles_, nselect_, noRTR_, verb_, HDF5_in) ! optional ins
        !
        use swav, only: calcJCoeffsPW
        use miet, only: calcMieIntCoeffs
        !
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        integer, intent(in) :: scheme, ncut(3)
        real(8), intent(in) :: wavelen(:)
        real(8), intent(in) :: ehost(size(wavelen))
        real(8), intent(in) :: geometry(:, :)
        real(8), intent(inout), optional :: &
            sig_oa_(:, :, :), & ! itype, n, iwave
            sig_(:, :, :, :), & ! itype, ijones, iwave, iorient
            sig_abs_(:, :, :, :, :) ! ishell, iscat, ijones, iwave, iorient
        complex(8), intent(in), optional :: &
            escat_(size(geometry, 2), 4, size(wavelen))
        ! NOTE: f2py fails for character(*) argument, so chose 256 chars
        character*256, intent(in), optional :: tfiles_(size(geometry, 2))
        integer, intent(in), optional :: nselect_(2, 2, 2, size(geometry, 2)), verb_
        logical, intent(in), optional :: noRTR_
        real(8), intent(inout), optional :: jsig_abs_oa(:, :, :)
        logical, intent(in), optional :: HDF5_in(:)
        ! Local variables
        character(*), parameter :: myName = 'spectrumFF'
        character*256 :: tfiles(size(geometry, 2))
        logical :: scatMiet(size(geometry, 2)), ldum, seek_oa, seek_fi
        logical :: multisplit = .false., noRTR
        real(8) :: hostK, rdum
        complex(8) :: s(0:4)
        integer :: iwav, ierr, j, js, verb, i, is, idum, idums, ns
        integer :: nmax, lmax, nscat, nwav, nfi = 0!, pol_type
        integer :: nselect(2, 2, 2, size(geometry, 2)), units(size(geometry, 2))
        real(8), allocatable :: inc(:, :)
        complex(8), allocatable :: TIJ(:, :), cJ(:, :, :), Tcol(:, :), cDum(:, :, :), &
                                   vtacs(:, :)
        integer :: tic, toc, tic2, toc2, tps
        real    :: t0, t1, t2, t3, t4, t6
        real(8) :: t5, t7
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        nscat = size(geometry, 2)
        nmax = ncut(1)
        lmax = 2*nmax*(nmax + 2)
        nwav = size(wavelen)
        rtol_G = 10.0d0**ncut(3)
        !
        !---------------------------------------------------
        ! Start processing arguments
        !---------------------------------------------------
        !
        if (size(geometry, 1) /= 8) then
            write (*, *) myName, '> ERROR: size(geometry,1) /= 8'
            STOP
        end if
        !
        tfiles = ''; units = 0
        if (present(tfiles_)) then
            do j = 1, nscat
                if (len_trim(tfiles_(j)) > 0) then
                    tfiles(j) = tfiles_(j)
                    units(j) = 1000 + j
                    iloop: do i = 1, j - 1
                        ! Check for filename duplicates and, if present,
                        ! store scatterer ID with the same T-matrix file.
                        ! Store -ve of the ID to distinguish from actual
                        ! file handles in "units".
                        if (tfiles(j) == tfiles(i)) then
                            units(j) = -i
                            exit iloop
                        end if
                    end do iloop
                end if
            end do
        end if
        !
        if (maxval(geometry(8, :)) >= tiny1) then
            ! non-Mie scatterer(s) present
            if (maxval(units) <= 0) then
                write (*, *) myName, &
                    '> ERROR: Missing tfiles for non-Mie scatterers'
                STOP
            end if
        end if
        !
        if (minval(geometry(8, :)) < tiny1) then
            ! Mie scatterer(s) present
            ldum = .false.
            if (present(escat_)) then
                if (maxval(abs(escat_)) > tiny2) then
                    ldum = .true. ! non-zero escat_ present
                end if
            end if
            if (.not. ldum) then
                write (*, *) myName, '> ERROR: Missing non-zero escat_'
                STOP
            end if
        end if
        !
        noRTR = .false.
        if (present(noRTR_)) noRTR = noRTR_
        !
        j = ncut(2)
        j = 2*j*(j + 2)
        allocate (TIJ(nscat*lmax, nscat*lmax), Tcol(j, j), &
                  vtacs(nscat*lmax, nscat*lmax))
        !TIJ=0
        !
        seek_oa = .false.
        if (present(sig_oa_) .and. scheme > 0) then
            if (size(sig_oa_, 1) == 6 .and. &
                size(sig_oa_, 2) == ncut(2) + 1 .and. &
                size(sig_oa_, 3) == size(wavelen)) seek_oa = .true.
        end if
        !
        seek_fi = .false.
        if (present(sig_)) then
            if (size(sig_, 1) == 3 .and. & ! ext, sca, abs
                size(sig_, 2) == 4 .and. & ! 4 Jone's vectors
                size(sig_, 3) == size(wavelen)) seek_fi = .true.
        end if
        !
        if (seek_fi) then
            nfi = size(sig_, 4)
            if (nfi == ncut(2)) then
                if (sum((sig_(1:3, 1, 1, ncut(2)) - sig_(1:3, 1, 1, 1))**2) < 1.0d-10) then
                    ! single incidence duplicated, designed to trigger splitting
                    ! of cross sections into multipolar contributions
                    multisplit = .true.
                    nfi = 1
                end if
            end if
            allocate (cJ(nscat*lmax, 2, 4*nfi), inc(4, 4*nfi))
            !pol_type = nint(sig_(1,2,1,1))
            !if(pol_type /= 1 .and. pol_type /= 2) then
            !   write(*,'(A,A)') myname,'> ERROR: Bad polarisation type'
            !endif
            is = 0
            do i = 1, nfi
                ! if(nint(sig_(1,2,1,i)) /= pol_type) then
                !    write(*,'(A,A)') myname,'> ERROR: variable polarisation'
                !    STOP
                ! endif
                do j = 1, 4 ! loop over all 4 Jone's vectors
                    is = is + 1
                    inc(1, is) = j
                    inc(2:4, is) = sig_(1:3, 1, 1, i)
                end do
            end do
            !
            if (present(sig_abs_)) then
                if (size(sig_abs_, 1) /= 4 .or. &
                    size(sig_abs_, 2) /= nscat .or. &
                    size(sig_abs_, 3) /= 4 .or. &
                    size(sig_abs_, 4) /= nwav .or. &
                    size(sig_abs_, 5) /= nfi) write (*, '(A,A)') myName, &
                    '> ERROR: sig_abs_ is present but misshapen!'
            end if
            !
        elseif (.not. seek_oa) then
            write (*, *) myName, '> ERROR: Nothing to seek.'
            STOP
        end if
        !
        nselect = 0 ! f2py will default to zero
        if (present(nselect_)) nselect = nselect_
        !
        verb = 0 ! default compatible with f2py
        if (present(verb_)) verb = verb_
        !
        !---------------------------------------------------------------
        ! Done processing/initialising arguments. Start calculations.
        !---------------------------------------------------------------
        !
        if (verb > 1) then
            call cpu_time(t0)
            call system_clock(tic)
        end if
        waves: do iwav = 1, nwav
            !
            if (verb > 0) then
                write (*, *)
                write (*, '(A,A,f8.2,A)') myname, '> ===== Wavelength: ', &
                    wavelen(iwav), ' (nm) ======================'
                write (*, *)
            end if
            !
            TIJ = 0
            !
            do j = 1, nscat
                js = (j - 1)*lmax
                if (geometry(8, j) < tiny1) then ! Mie scatterer
                    scatMiet(j) = .true.
                    do i = 1, 1 - nint(geometry(8, j))
                        TIJ(js + i, js + i) = escat_(j, i, iwav)
                    end do
                else
                    scatMiet(j) = .false.
                    if (units(j) > 0) then
                        call readTmatFile( &
                            Tmat=TIJ(js + 1:js + lmax, js + 1:js + lmax), &
                            filename=trim(tfiles(j)), &
                            unit=units(j), &
                            wavelen=wavelen(iwav), &
                            HDF5_in_ = HDF5_in(j))
                    elseif (units(j) < 0) then
                        i = abs(units(j))
                        is = (i - 1)*lmax
                        TIJ(js + 1:js + lmax, js + 1:js + lmax) = &
                            TIJ(is + 1:is + lmax, is + 1:is + lmax)
                    else
                        write (*, '(A,A,i5)') &
                            myName, '> ERROR: units(j)=0 for j=', j
                        STOP
                    end if
                end if
            end do
            !
            if (seek_fi) then
                !
                cJ = 0
                cJ(1:4, 1, :) = inc(1:4, :)
                if (verb > 1) then
                    call cpu_time(t2)
                    call system_clock(tic2)
                end if

                call solve( &
                    wavelen=wavelen(iwav), &
                    ehost=ehost(iwav), &
                    geometry=geometry, &
                    TIJ=TIJ, &
                    scheme_=scheme, &
                    nselect_=nselect, &
                    verb_=verb, &
                    cJ_=cJ, &
                    noRTR_=noRTR, &
                    ierr_=ierr)
                if (verb > 1) then
                    call cpu_time(t3)
                    call system_clock(toc2, count_rate=tps)
                    write (*, '(A,A,2(1x,es10.3E2))') &
                        myName, '> Calculation time [Solve: diff. incs] (CPU & real in s): ', &
                        t3 - t2, dble(toc2 - tic2)/dble(tps)
                end if
                !

            elseif (seek_oa) then
                !
                if (verb > 1) then
                    call cpu_time(t2)
                    call system_clock(tic2)
                end if
                call solve( &
                    wavelen=wavelen(iwav), &
                    ehost=ehost(iwav), &
                    geometry=geometry, &
                    TIJ=TIJ, &
                    scheme_=scheme, &
                    nselect_=nselect, &
                    verb_=verb, &
                    noRTR_=noRTR, &
                    ierr_=ierr)

                if (verb > 1) then
                    call cpu_time(t3)
                    call system_clock(toc2, count_rate=tps)
                    write (*, '(A,A,2(1x,es10.3E2))') &
                        myName, '> Calculation time [Solve: OA] (CPU & real in s): ', &
                        t3 - t2, dble(toc2 - tic2)/dble(tps)
                end if
                !
            else
                write (*, '(A,A)') myname, '> ERROR: Not seeking solution?'
            end if
            !
            if (ierr /= 0) then
                !
                write (*, '(A,A,i2)') &
                    myname, '> ERROR: solve returned ierr= ', ierr
                sig_oa_(:, :, iwav) = -ierr
                STOP
                !
            end if
            !
            ! compute rotationally averaged cross-sections
            !
            ! wavenumber k_host
            hostK = tpi/wavelen(iwav)*sqrt(ehost(iwav))
            if (seek_oa) then
                sig_oa_(1, 1, iwav) = hostK
                Tcol(1, 1) = cmplx(wavelen(iwav), ehost(iwav), kind(Tcol))
                j = lmax*nscat
                if (scheme == 3) then
                    j = lmax
                    ldum = .true. ! need this for single scatterer!
                else
                    ldum = .false.
                end if
                if (ldum .and. cs_stout_G .and. verb > 1) write (*, '(A,A)') myName, &
                    '> NOTE: The orientation averaged absorption cross-section &
                        &for each particle is not be calculated with scheme = 3'
                if (.not. ldum .and. cs_stout_G .and. (minval(geometry(8, :)) < -tiny1) .and. verb > 1) &
                    write (*, '(A,A)') myName, &
                    '> NOTE: The orientation averaged absorption cross-section &
                        &for each particle is not calculated when particles are not homogeneous spheres'
                if (verb > 1) then
                    call cpu_time(t2)
                    call system_clock(tic2)
                end if
                !
                call contractTmat( &
                    Tin=TIJ(:, 1:j), &
                    scatXYZR=geometry(1:4, :), &
                    Tout=Tcol, &
                    rtr=.not. noRTR, &
                    verb_=verb, &
                    mack_=ldum)
                !
                if (verb > 1) then
                    call cpu_time(t3)
                    call system_clock(toc2, count_rate=tps)
                    write (*, '(A,A,2(1x,es10.3E2))') &
                        myName, '> Calculation time [contractTmat] (CPU & real in s): ', &
                        t3 - t2, dble(toc2 - tic2)/dble(tps)
                end if
                if (.not. ldum .and. cs_stout_G .and. (minval(geometry(8, :)) > -tiny1)) then
                    if (verb > 1) then
                        call cpu_time(t2)
                        call system_clock(tic2)
                    end if
                    !
                    call calcOaStout( &
                        TIJ=TIJ, &
                        scatXYZR=geometry(1:4, :), &
                        sigOA=sig_oa_(1:3, 1, iwav), &
                        cdOA_=sig_oa_(4, 1, iwav), &
                        jAbsOA=jsig_abs_oa(:, 1, iwav), &
                        escat=escat_(:, 1, iwav), &
                        ehost=ehost(iwav), &
                        vtacsJK=vtacs, &
                        verb_=verb)
                    !
                    if (verb > 1) then
                        call cpu_time(t3)
                        call system_clock(toc2, count_rate=tps)
                        write (*, '(A,A,2(1x,es10.3E2))') &
                            myName, '> Calculation time [calcOaStout] (CPU & real in s): ', &
                            t3 - t2, dble(toc2 - tic2)/dble(tps)
                    end if
                else
                    !
                    if (verb > 1) then
                        call cpu_time(t2)
                        call system_clock(tic2)
                    end if
                    !
                    call calcOAprops( &
                        Tmat=Tcol, &
                        sigOA=sig_oa_(:, :, iwav), &
                        verb_=verb, &
                        rtol_=rtol_G)

                    !
                    if (verb > 1) then
                        call cpu_time(t3)
                        call system_clock(toc2, count_rate=tps)
                        write (*, '(A,A,2(1x,es10.3E2))') &
                            myName, '> Calculation time [calcOAprops] (CPU & real in s): ', &
                            t3 - t2, dble(toc2 - tic2)/dble(tps)
                    end if

                end if
            end if
            !
            t2 = 0.0d0
            t3 = 0.0d0
            t4 = 0.0d0
            t5 = 0.0d0
            t6 = 0.0d0
            t7 = 0.0d0
            if (seek_fi) then ! fixed-inc cross-sections
                !
                is = 0
                inc_loop: do i = 1, nfi
                    !
                    jones_loop: do j = 1, 4
                        is = is + 1
                        sig_(1, j, iwav, i) = hostK
                        if (multisplit) then
                            js = ncut(2) ! recall that nfi=1
                        else
                            js = i
                        end if
                        if (verb > 1) then
                            call cpu_time(t2)
                            call system_clock(tic2)
                        end if
                        !
                        if (cs_stout_G) then
                            call calcCsStout( &
                                scatXYZR=geometry(1:4, :), &
                                aJ=cJ(:, 1, is), &
                                fJ=cJ(:, 2, is), &
                                sig=sig_(1:3, j, iwav, i:js), &
                                nmax2_=ncut(2), &
                                tol_=rtol_G, &
                                vtacs_t=vtacs, &
                                verb_=verb)
                            if (verb > 1) then
                                call cpu_time(t3)
                                call system_clock(toc2, count_rate=tps)
                                t4 = t4 + t3 - t2
                                t5 = t5 + dble(toc2 - tic2)/dble(tps)
                            end if

                        else
                            call calcCs( &
                                scatXYZR=geometry(1:4, :), &
                                inc=inc(:, is), &
                                fJ=cJ(:, 2, is), &
                                sig=sig_(1:3, j, iwav, i:js), &
                                nmax2_=ncut(2), &
                                tol_=rtol_G, &
                                verb_=verb)

                            !
                            if (verb > 1) then
                                call cpu_time(t3)
                                call system_clock(toc2, count_rate=tps)
                                t4 = t4 + t3 - t2
                                t5 = t5 + dble(toc2 - tic2)/dble(tps)
                            end if
                        end if
                        !
                        if (present(sig_abs_)) then
                            !
                            if (verb > 1) then
                                call cpu_time(t2)
                                call system_clock(tic2)
                            end if
                            ! Seek (partial) absorptions for each Mie scatterer
                            !
                            rdum = tpi/wavelen(iwav)
                            allocate (cDum(lmax, 4, 2))
                            do idum = 1, nscat
                                if (scatMiet(idum)) then
                                    idums = (idum - 1)*lmax
                                    ns = nint(abs(geometry(8, idum))) ! shell count
                                    s(0) = hostK
                                    s(1:1 + ns) = rdum*sqrt(escat_(idum, 1:1 + ns, iwav))
                                    ! reversed interface indexing when callin
                                    call calcMieIntCoeffs( &
                                        a=geometry(4 + ns:4:-1, idum), &
                                        k=s(1 + ns:0:-1), & ! zero-index: hostK
                                        scaCoeffs=cJ(idums + 1:idums + lmax, 2, is), &
                                        intCoeffsReg=cDum(1:lmax, 1 + ns:1:-1, 1), &
                                        intCoeffsIrr=cDum(1:lmax, 1 + ns:1:-1, 2), &
                                        csAbs=sig_abs_(1 + ns:1:-1, idum, j, iwav, i))
                                end if
                            end do
                            deallocate (cDum)
                            if (verb > 1) then
                                call cpu_time(t3)
                                call system_clock(toc2, count_rate=tps)
                                t6 = t6 + t3 - t2
                                t7 = t7 + dble(toc2 - tic2)/dble(tps)
                            end if
                        end if
                        !
                    end do jones_loop
                end do inc_loop
                !
            end if
            !
            if (verb > 1) then
                if (cs_stout_G) then
                    write (*, '(A,A,2(1x,es10.3E2))') &
                        myName, '> Calculation time [calcCsStout] (CPU & real in s): ', &
                        t4, t5
                else
                    write (*, '(A,A,2(1x,es10.3E2))') &
                        myName, '> Calculation time [calcCs] (CPU & real in s): ', &
                        t4, t5
                end if
                if (present(sig_abs_)) then
                    write (*, '(A,A,2(1x,es10.3E2))') &
                        myName, '> Calculation time [partial abs. for Mie scatterer] (CPU & real in s): ', &
                        t6, t7
                end if
            end if
        end do waves
        if (verb > 1) then
            call cpu_time(t1)
            call system_clock(toc, count_rate=tps)
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time (CPU & real in s): ', &
                t1 - t0, dble(toc - tic)/dble(tps)
        end if
        !
        if (allocated(TIJ)) deallocate (TIJ)
        if (allocated(Tcol)) deallocate (Tcol)
        if (allocated(cJ)) deallocate (cJ)
        if (allocated(inc)) deallocate (inc)
        !b
    end subroutine spectrumFF
    !
    !=================================================================
    !
    subroutine solve( &
        wavelen, ehost, geometry, & ! mandatory ins
        TIJ, cJ_, cJint_, csAbs_, & ! all in/outs
        nselect_, scheme_, noRTR_, verb_, & ! optional ins
        ierr_, dipoles) ! optional out, becomes the return value when wrapped by f2py
        !

        !==============================================================
        ! This routine is the crux of TERMS, solving a given multi-scat
        ! problm by operating in a specified scheme. The problem
        ! is defined by a bunch of mandatory inputs, some of which are
        ! modified to produce an output. The only genuine output (i.e.
        ! return vlue) is an integer error code:
        !
        ! ierr_:  0 -> all good
        !         1 -> error in arguments processing
        !         2 -> error in prestaging, staging, or solving/inverting Ax=b
        !==============================================================
        !
        use swav, only: calcSphBessels, calcWignerBigD, calcJCoeffsPW
        use miet, only: calcMieTMat, calcMieIntCoeffs
        use linalg, only: invSqrMat, solLinSys
        ! use swav_dipole, only : calcCoeffsJ4dipole
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables
        ! NOTE: f2py converts subroutines to functions, and variables passed with
        ! intent OUT become (part of) the output of the converted function. The
        ! python function can be forced to behave like a subroutine (i.e so it
        ! modifies the arguments without producing an output) by changing the
        ! declared intent from OUT to INOUT.
        real(8), intent(in) :: wavelen
        real(8), intent(in) :: ehost
        real(8), intent(in) :: geometry(:, :)
        complex(8), intent(inout) :: TIJ(:, :)
        complex(8), intent(inout), optional :: cJ_(:, :, :)
        complex(8), intent(inout), optional :: cJint_(size(TIJ, 1), 4, 2)
        real(8), intent(inout), optional :: csAbs_(size(geometry, 2), 4)
        integer, intent(in), optional :: nselect_(2, 2, 2, size(geometry, 2))
        integer, intent(in), optional :: scheme_, verb_
        logical, intent(in), optional :: noRTR_
        integer, intent(out), optional :: ierr_ ! integer error code
        real(8), intent(in), optional :: dipoles(:, :)
        ! Local variables
        character(*), parameter :: myName = 'solve'
        logical :: cJpresent, rtr, scatMiet(size(geometry, 2))
        logical :: mask_rows(2, 2), mask_cols(2, 2), PWinc = .false.
        integer :: i, j, js, k, ks, n, l, p, q, qq, q2, qq2, hi, lo, lmax, pmax, nmax
        integer :: scheme, nscat, nhi(2, 2), nlo(2, 2), verb
        real(8) :: sigExt, sigSca, rdum, inc(4)
        integer :: tic, toc, tps
        real(8) :: t0, t1, t2(22)
        complex(8) :: zdum, scatK(4, size(geometry, 2)), x(4) = 0, s(0:4) = 0
        ! complex(8) :: mat(size(cJ_,1),size(cJ_,3))
        complex(8), allocatable :: work(:, :, :)
        complex(8), allocatable :: ipwE0(:, :)!, idpE0(:,:)
        real(8), allocatable :: ipwDirn(:, :)!, idpDirn(:,:)
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        if (present(ierr_)) ierr_ = 0 ! presumption of innocence
        !
        !---------------------------------------------------
        ! Start processing arguments and initialise
        !---------------------------------------------------
        if (present(dipoles)) then
        else
            PWinc = .true.
        end if
        if (size(TIJ, 1) /= size(TIJ, 2)) then
            write (*, '(A,A)') myName, '> ERROR: Array TIJ not square!'
            if (present(ierr_)) ierr_ = 1
            return
        end if
        !
        if (size(geometry, 1) /= 8) then
            write (*, '(A,A)') myName, '> ERROR: size(geometry,1) /= 8'
            if (present(ierr_)) ierr_ = 1
            return
        end if
        !
        nscat = size(geometry, 2)
        if (mod(size(TIJ, 1), 2*nscat) /= 0) &
            write (*, '(A,A)') myName, '> WARNING: mod(size(TIJ,1),2*nscat) /= 0'
        lmax = size(TIJ, 1)/nscat
        pmax = lmax/2
        nmax = int(sqrt(real(pmax)))
        if (nmax < 1) then
            write (*, '(A,A)') myName, '> ERROR: nmax < 1'
            if (present(ierr_)) ierr_ = 1
            return
        elseif (pmax /= nmax*(nmax + 2)) then
            write (*, '(A,A)') myName, '> ERROR: pmax /= nmax*(nmax+2)'
            if (present(ierr_)) ierr_ = 1
            return
        end if
        !
        i = minloc(geometry(8, :), 1)
        if (geometry(8, i) < tiny1 .and. nint(geometry(8, i)) > 3) then
            write (*, '(A,A)') myName, &
                '> ERROR: Cannot handle more than 3 shells'
            if (present(ierr_)) ierr_ = 1
            return
        end if
        !
        allocate (bes4bal_G(0:nmax, 0:1, 1:nscat))

        ! Make missing optional arguments compatible with f2py's default
        verb = 0
        if (present(verb_)) verb = verb_
        verb_G = verb ! should be phased out..
        scheme = 0
        if (present(scheme_)) scheme = scheme_
        rtr = .true.
        if (present(noRTR_)) rtr = .not. noRTR_
        !
        t2 = 0.0d0
        ! Set the host wavenumber here
        hostK_G = (tpi/wavelen)*sqrt(ehost)
        !
        cJpresent = .false.
        if (present(cJ_)) then
            if (size(cJ_, 1) == size(TIJ, 1) .and. &
                size(cJ_, 2) == 2) cJpresent = .true.

        end if
        if (cJpresent .AND. PWinc) then
            allocate (ipwDirn(3, size(cJ_, 3)), ipwE0(3, size(cJ_, 3)))
            do i = 1, size(cJ_, 3)
                do j = 1, 4
                    inc(j) = realpart(cJ_(j, 1, i))
                end do
                call parseInc( &
                    inc=inc, &
                    inc_dirn=ipwDirn(:, i), &
                    inc_ampl=ipwE0(:, i), &
                    verb_=verb &
                    )
                call calcJCoeffsPW( &
                    ipwE0=ipwE0(:, i), &
                    kVec=hostK_G*ipwDirn(:, i), &
                    ipwCoeffsJ=cJ_(:, 1, i), &
                    xyz=geometry(1:3, :) &
                    )
            end do
        elseif (cJpresent .AND. (.not. PWinc)) then
            ! allocate(idpDirn(3,size(cJ_,3)),idpE0(3,size(cJ_,3)))
            ! do i = 1,size(cJ_,3)
            !   do j=1,4
            !     inc(j) = realpart(cJ_(j,1,i))
            ! enddo
            ! call parseInc( &
            !    inc = inc, &
            !    inc_dirn = ipwDirn(:,i), &
            !    inc_ampl = ipwE0(:,i), &
            !    verb_ = verb &
            !    )

            !call calcCoeffsJ4dipole( &
            !  xyz_d = dipoles (1:3,i), &
            !  xyz =geometry(1:3,:), &
            !  idCoeffsJ = cJ_(:,1,i), &
            !  ehost=ehost, &
            !  k = cmplx(hostK_G,0,kind(hostK_G)) )
            ! enddo

        end if
        !
        !----------------------------------------------------------------
        ! End of arguments processing/initialisation. Start calculations.
        !----------------------------------------------------------------
        !
        ! >>> Prestaging >>>>>>>>>
        if (verb > 0) write (*, '(A,A)') myname, '> Prestaging...'
        !
        scatK = 0
        rdum = tpi/wavelen
        !
        ! Loop over scatterers to compute the balancing weights and
        ! wave-numbers, and to stage the one-body T-matrices.
        jscat: do j = 1, nscat
            !
            ! Compute balancing weights (wavelength dependent) and scatK
            zdum = cmplx(balScale_G*hostK_G*geometry(4, j), 0, 8)
            !write(*,*) 'TEST: hostK_G= ', hostK_G,' zdum= ', zdum
            !zdum = cmplx(hostK_G*geometry(4,j)*sqrt(1.0d0-1.0d0/geometry(8,j)**2),0,8)
            !zdum = cmplx(hostK_G*geometry(4,j)*(2.0d0-1.0d0/geometry(8,j)),0,8)
            if (verb > 1) then
                call cpu_time(t0)
                call system_clock(tic)
            end if

            call calcSphBessels( &
                z=zdum, &
                nmax=nmax, &
                regt=.false., &
                bes=bes4bal_G(:, 0, j)) ! Irregular
            call calcSphBessels( &
                z=zdum, &
                nmax=nmax, &
                regt=.true., &
                bes=bes4bal_G(:, 1, j)) ! Regular

            if (verb > 1) then
                call cpu_time(t1)
                call system_clock(toc, count_rate=tps)
                t2(1) = t2(1) + t1 - t0
                t2(2) = t2(2) + dble(toc - tic)/dble(tps)
            end if
            ! Convert spherical bessels to riccati-bessels
            bes4bal_G(:, :, j) = zdum*bes4bal_G(:, :, j)
            !
            if (present(nselect_)) then
                nlo = nselect_(1, :, :, j)
                nhi = nselect_(2, :, :, j)
            else ! f2py-safe
                nlo = 0
                nhi = 0
            end if
            do q2 = 1, 2
                do q = 1, 2
                    if (nlo(q, q2) < 0) then
                        mask_rows(q, q2) = .false.
                        nlo(q, q2) = abs(nlo(q, q2))
                    else
                        mask_rows(q, q2) = .true.
                    end if
                    if (nhi(q, q2) < 0) then
                        mask_cols(q, q2) = .false.
                        nhi(q, q2) = abs(nhi(q, q2))
                    else
                        mask_cols(q, q2) = .true.
                    end if
                    if (nlo(q, q2) < 1 .or. nlo(q, q2) > nmax) nlo(q, q2) = 1
                    if (nhi(q, q2) < 1 .or. nhi(q, q2) > nmax) nhi(q, q2) = nmax
                end do
            end do
            !
            ! Fill TIJ's diagonal blocks with either the supplied T-matrix
            ! or a calculated Mie T-matrix for the current wavelength.
            js = (j - 1)*lmax
            if (geometry(8, j) < tiny1) then ! Use Mie theory
                !
                i = nint(abs(geometry(8, j))) ! shell count
                ! 0 -> homogeneous sphere
                !-1 -> singly coated sphere (core with 1 shell)
                !-2 -> doubly coated sphere (core with 2 shells)
                !-3 -> triply coated sphere
                ! WARNING: Tricky indexing of K interfaces and K+1 regions!!
                !
                !write(*,*) 'TEST: shell count= ', i
                scatMiet(j) = .true.
                do l = 1, 1 + i

                    scatK(l, j) = rdum*sqrt(TIJ(js + l, js + l)) !escat(l,j)=TIJ(js+l,js+l)
                    !write(*,*) 'TEST: ',scatK(l,j)
                end do
                ! Calculate dimensionless size parameters (x_k = a_k*kOut) and
                ! relative refractive indices (s_k = kIn/kOut) for all spherical
                ! interfaces (indexed by k), with the indexing reverse to that in
                ! Eric and Pablo's book pp.622-3; i.e. here increasing k corresponds
                ! to decreasing (as opposed to increasing) radius.
                !
                x(1) = cmplx(geometry(4, j)*hostK_G, 0, 8) ! a*kOut for outer interface
                s(1) = scatK(1, j)/hostK_G           ! kIn/kOut for outer interface
                if (i > 0) then ! go inwards
                    x(2:1 + i) = geometry(5:4 + i, j)*scatK(1:i, j)  ! a*kOut
                    s(2:1 + i) = scatK(2:1 + i, j)/scatK(1:i, j) ! kIn/kOut
                end if
                !
                if (verb > 1) then
                    call cpu_time(t0)
                    call system_clock(tic)
                end if
                call calcMieTMat( &
                    x=x(1 + i:1:-1), &  ! reverse order to match Eric and Pablo
                    s=s(1 + i:1:-1), &
                    tmat=TIJ(js + 1:js + lmax, js + 1:js + lmax), &
                    zeropad_=.true.)
                !
                if (verb > 1) then
                    call cpu_time(t1)
                    call system_clock(toc, count_rate=tps)
                    t2(3) = t2(3) + t1 - t0
                    t2(4) = t2(4) + dble(toc - tic)/dble(tps)
                end if
                ! zero multipoles not in range [nlo,nhi] for EE and MM
                !write(*,*) 'TEST: nmax, pmax, lmax=', nmax, pmax, lmax
                do q = 1, 2
                    qq = (q - 1)*pmax
                    if (nlo(q, q) > 1) then ! mask multipoles 1 ... nlo-1
                        n = nlo(q, q) - 1
                        p = n*(n + 2)
                        TIJ(js + qq + 1:js + qq + p, js + qq + 1:js + qq + p) = 0
                    end if
                    if (nhi(q, q) < nmax) then ! mask multipoles nhi+1 ... nmax
                        n = nhi(q, q) !+1 Don't nee dadding 1
                        p = n*(n + 2)
                        TIJ(js + qq + p + 1:js + qq + pmax, js + qq + p + 1:js + qq + pmax) = 0
                    end if
                    if (verb > 1) write (*, '(A,A,i3,A,i1,1x,i2,1x,i2)') &
                        myname, '> Prestaged Mie T-matrix for scatterer ', j, &
                        ' with q, nlo, nhi = ', q, nlo(q, q), nhi(q, q)
                end do
                !
            else ! general case
                !
                scatMiet(j) = .false. ! non-mie
                sigExt = 0
                do k = 1, lmax
                    sigExt = sigExt + realpart(TIJ(js + k, js + k))
                end do ! For sanity check
                if (verb > 1) write (*, '(A,A,i3,A,es10.4e1)') &
                    myname, '> General T-matrix for scatterer ', j, &
                    ' with Re{Tr(T)}= ', sigExt
                !
                ! check if rotation is required
                rotcheck: if (abs(geometry(5, j)) > tiny1 .or. &
                              abs(geometry(6, j)) > tiny1 .or. &
                              abs(geometry(7, j)) > tiny1) then
                    !
                    if (verb > 1) write (*, '(A,A,3(1x,es9.3E1))') &
                        myname, '> Rotating by Euler angles= ', geometry(5:7, j)
                    !
                    allocate (work(lmax, lmax, 1)) ! the temporary rotation matrix R
                    work = 0
                    if (verb > 1) then
                        call cpu_time(t0)
                        call system_clock(tic)
                    end if
                    call calcWignerBigD( &
                        angles=(/-geometry(7, j), -geometry(6, j), -geometry(5, j)/), &
                        ! i.e.: (-gamma, -beta, -alpha), with the Euler angles being
                        ! alpha = phi, beta = theta, gamma = 0 for tilted spheroids
                        pmax=pmax, &
                        bigD=work(1:pmax, 1:pmax, 1))

                    if (verb > 1) then
                        call cpu_time(t1)
                        call system_clock(toc, count_rate=tps)
                        t2(5) = t2(5) + t1 - t0
                        t2(6) = t2(6) + dble(toc - tic)/dble(tps)
                    end if

                    work(pmax + 1:2*pmax, pmax + 1:2*pmax, 1) = &
                        work(1:pmax, 1:pmax, 1) ! this contains inv(R)
                    ! Right-multiply T by block-diagonal inv(R):  T <-- T inv(R)
                    do q = 0, 1
                        qq = q*pmax
                        do n = 1, nmax
                            lo = n*n
                            hi = lo + 2*n
                            TIJ(js + 1:js + lmax, js + lo + qq:js + hi + qq) = matmul( &
                                                                               TIJ(js + 1:js + lmax, js + lo + qq:js + hi + qq), &
                                                                               work(lo:hi, lo:hi, 1))
                        end do
                    end do
                    ! Left-multiply T by block-dagonal R=transpose(conjg(inv(R))),
                    ! i.e. update T <-- RT, which we do by transposing both
                    ! sides and right-multiplying, i.e.
                    !          transpose(T) <-- transpose(T)*conjg(inv(R)),
                    ! so that longer strides are along 1st index when multiplying.
                    TIJ(js + 1:js + lmax, js + 1:js + lmax) = &
                        transpose(TIJ(js + 1:js + lmax, js + 1:js + lmax))
                    do q = 0, 1
                        qq = q*pmax
                        do n = 1, nmax
                            lo = n*n
                            hi = lo + 2*n
                            TIJ(js + 1:js + lmax, js + lo + qq:js + hi + qq) = matmul( &
                                                                               TIJ(js + 1:js + lmax, js + lo + qq:js + hi + qq), &
                                                                               conjg(work(lo:hi, lo:hi, 1)))
                            ! TIJ(js+lo+qq:js+hi+qq, js+1:js+lmax_G) = matmul( &
                            !      conjg(transpose(work(lo:hi,lo:hi,1))), &
                            !      TIJ(js+lo+qq:js+hi+qq,js+1:js+lmax_G))
                        end do
                    end do
                    TIJ(js + 1:js + lmax, js + 1:js + lmax) = &
                        transpose(TIJ(js + 1:js + lmax, js + 1:js + lmax))
                    deallocate (work)
                    !
                    ! Check that trace is unchanged by rotation
                    sigSca = 0
                    do k = 1, lmax
                        sigSca = sigSca + realpart(TIJ(js + k, js + k))
                    end do
                    if (abs(sigSca - sigExt)/abs(sigExt) > tiny) then
                        write (*, '(A,A,es10.4E1)') myName, &
                            '> ERROR: After rotation Re{Trace(T)}= ', sigSca
                        if (present(ierr_)) ierr_ = 2
                        return
                    end if
                    !
                end if rotcheck
                !
                if (verb > 1) write (*, '(A,A,i3)') &
                    myname, '> Prestaged general T-matrix for scatterer ', j

                ! mask multipoles outside the range [nlo,nhi]
                do q = 1, 2
                    qq = (q - 1)*pmax
                    do q2 = 1, 2
                        qq2 = (q2 - 1)*pmax
                        if (nlo(q, q2) > 1) then ! mask multipoles 1 ... nlo-1
                            n = nlo(q, q2) - 1
                            p = n*(n + 2)
                            if (mask_cols(q, q2)) then
                                TIJ(js + qq + 1:js + qq + p, js + qq2 + 1:js + qq2 + pmax) = 0
                            end if
                            if (mask_rows(q, q2)) then
                                TIJ(js + qq + 1:js + qq + pmax, js + qq2 + 1:js + qq2 + p) = 0
                            end if
                        end if
                        if (nhi(q, q2) < nmax) then ! mask multipoles nhi+1 ... nmax
                            n = nhi(q, q2) !+1 Don't need adding 1
                            p = n*(n + 2)
                            if (mask_cols(q, q2)) then
                                TIJ(js + qq + p + 1:js + qq + pmax, js + qq2 + 1:js + qq2 + pmax) = 0
                            end if
                            if (mask_rows(q, q2)) then
                                TIJ(js + qq + 1:js + qq + pmax, js + qq2 + p + 1:js + qq2 + pmax) = 0
                            end if
                        end if
                    end do
                end do
                !
                ! if(nlo > 1) then ! mask multipoles 1 ... nlo-1
                !    k = (nlo-1)*(nlo+1)
                !    TIJ(js+1:js+k,js+1:js+lmax) = 0
                !    TIJ(js+1+pmax:js+k+pmax,js+1:js+lmax) = 0
                !    TIJ(js+1:js+lmax, js+1:js+k) = 0
                !    TIJ(js+1+pmax:js+k+pmax,js+1:js+lmax) = 0
                ! endif
                ! if(nhi < nmax) then ! mask multipoles nhi+1 ... nmax
                !    k = (nhi+1)*(nhi+3)
                !    TIJ(js+k:js+pmax, js+1:js+lmax) = 0
                !    TIJ(js+k+pmax:js+lmax, js+1:js+lmax) = 0
                !    TIJ(js+1:js+lmax, js+k:js+pmax) = 0
                !    TIJ(js+1:js+lmax, js+k+pmax:js+lmax) = 0
                ! endif
                !
                ! if(verb > 1) write(*,'(A,A,i3,A,i3,1x,i3)') &
                !      myname,'> Prestaged general T-matrix for scatterer ', j, &
                !      ' with nlo, nhi = ',nlo,nhi
                !
            end if
            !
        end do jscat
        !
        if (dumpPrestagedA_G) then
            call dumpMatrix(mat=TIJ, &
                            ofile='prestagedA.txt', &
                            verb_=verb)
        end if
        !
        ! <<< End of prestaging <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        !
        !
        ! >>> Start staging and solving/inverting linear system >>>>>>>>>
        if (verb > 0) write (*, '(A,A)') myname, '> Staging and solving/inverting...'
        !
        if (scheme == 0) then
            !
            ! Do not seek the collective T-matrices and just solve the
            ! linear system Ax=b for the given incidence.
            !
            allocate (work(lmax, lmax, nscat))
            !
            if (verb > 1) then
                call cpu_time(t0)
                call system_clock(tic)
            end if
            call stageAmat( &
                A=TIJ, &
                scatXYZ=geometry(1:3, :), &
                scatMiet=scatMiet, &
                rtr=rtr, &
                right_=.false., &
                Tmats_=work, &   !T matrix of each scatterer.
                balance_=balStout_G, &
                verb_=verb)
            if (verb > 1) then
                call cpu_time(t1)
                call system_clock(toc, count_rate=tps)
                t2(7) = t2(7) + t1 - t0
                t2(8) = t2(8) + dble(toc - tic)/dble(tps)
            end if
            !
            if (balStout_G) then
                do i = 1, size(cJ_, 3)
                    do j = 1, nscat
                        js = (j - 1)*lmax
                        if (verb > 1) then
                            call cpu_time(t0)
                            call system_clock(tic)
                        end if
                        call balanceVecJ(j=j, jregt=.true., &
                                         Vec=cJ_(js + 1:js + lmax, 1, i))
                        !
                        if (verb > 1) then
                            call cpu_time(t1)
                            call system_clock(toc, count_rate=tps)
                            t2(9) = t2(9) + t1 - t0
                            t2(10) = t2(10) + dble(toc - tic)/dble(tps)
                        end if
                    end do
                end do
                if (verb > 1) write (*, '(A,A)') myname, '> Staged balanced system(s) Ax=b'
            elseif (verb > 1) then
                write (*, '(A,A)') myname, '> Staged unbalanced system(s) Ax=b'
            end if
            !
            cJ_(:, 2, :) = cJ_(:, 1, :)
            if (verb > 1) then
                call cpu_time(t0)
                call system_clock(tic)
            end if
            call solLinSys( &
                A=TIJ, & ! Note: TIJ is changed!
                X=cJ_(:, 2, :), &   ! excitation coefficients
                isol_=isolve_G, &
                verb_=verb)
            !
            if (verb > 1) then
                call cpu_time(t1)
                call system_clock(toc, count_rate=tps)
                t2(11) = t2(11) + t1 - t0
                t2(12) = t2(12) + dble(toc - tic)/dble(tps)
            end if
            if (verb > 1) then
                if (balStout_G) then
                    write (*, '(A,A)') myname, &
                        '> Solved balanced Ax=b for x without inverting A'
                else
                    write (*, '(A,A)') myname, &
                        '> Solved unbalanced Ax=b for x without inverting A'
                end if
            end if
            !
            ! Have (potentially balanced) excitation coefficients.
            ! Unbalance them first, if required, and then left-multiply
            ! by (unbalanced!) Tmats to get the scattering coefficients.
            ! The incident coefficients are also unbalanced here (if needed).
            do i = 1, size(cJ_, 3) ! incidences
                do j = 1, nscat
                    js = (j - 1)*lmax
                    if (balStout_G) then
                        if (verb > 1) then
                            call cpu_time(t0)
                            call system_clock(tic)
                        end if
                        call balanceVecJ(j=j, jregt=.true., &
                                         Vec=cJ_(js + 1:js + lmax, 1, i), rev_=.true.)
                        call balanceVecJ(j=j, jregt=.true., &
                                         Vec=cJ_(js + 1:js + lmax, 2, i), rev_=.true.)
                        if (verb > 1) then
                            call cpu_time(t1)
                            call system_clock(toc, count_rate=tps)
                            t2(9) = t2(9) + t1 - t0
                            t2(10) = t2(10) + dble(toc - tic)/dble(tps)
                        end if

                    end if
                    if (scatMiet(j)) then
                        do n = 1, lmax
                            cJ_(js + n, 2, i) = work(n, n, j)*cJ_(js + n, 2, i)
                        end do
                    else
                        cJ_(js + 1:js + lmax, 2, i) = &  !scattering coefficients
                            matmul(work(:, :, j), cJ_(js + 1:js + lmax, 2, i))   !Eq. 39
                    end if
                end do
            end do
            !
            deallocate (work)
            !
        elseif (scheme == 1) then
            !
            ! Determine the collective T-matrices by direct inversion
            !
            allocate (work(lmax, lmax, nscat))
            !
            if (verb > 1) then
                call cpu_time(t0)
                call system_clock(tic)
            end if
            call stageAmat( &
                A=TIJ, &
                scatXYZ=geometry(1:3, :), &
                scatMiet=scatMiet, &
                rtr=rtr, &
                right_=.false., &
                Tmats_=work, &
                balance_=balStout_G, &
                verb_=verb)
            if (verb > 1) then
                call cpu_time(t1)
                call system_clock(toc, count_rate=tps)
                t2(7) = t2(7) + t1 - t0
                t2(8) = t2(8) + dble(toc - tic)/dble(tps)
            end if
            if (verb > 1) write (*, '(A,A)') myname, &
                '> Staged matrix A for the balanced system Ax=b'
            !
            if (verb > 1) then
                call cpu_time(t0)
                call system_clock(tic)
            end if
            call invSqrMat( &
                A=TIJ, &
                trans_=transInv_G, &
                verb_=verb)
            !
            if (verb > 1) then
                call cpu_time(t1)
                call system_clock(toc, count_rate=tps)
                t2(13) = t2(13) + t1 - t0
                t2(14) = t2(14) + dble(toc - tic)/dble(tps)
            end if
            if (verb > 1) write (*, '(A,A)') myname, &
                '> Inverted matrix A for the balanced systm Ax=b'
            !
            ! To get the actual TIJ matrices, still need to unbalance and
            ! left-multiply by the block-diagonal matrix with (unbalanced!)
            ! 1-body T-matrices.
            do j = 1, nscat
                js = (j - 1)*lmax
                if (balStout_G) then
                    do k = 1, nscat
                        ks = (k - 1)*lmax
                        if (verb > 1) then
                            call cpu_time(t0)
                            call system_clock(tic)
                        end if
                        call balanceMatJI(j=j, jregt=.true., &
                                          Mat=TIJ(js + 1:js + lmax, ks + 1:ks + lmax), &
                                          i=k, iregt=.true., rev_=.true.)
                        if (verb > 1) then
                            call cpu_time(t1)
                            call system_clock(toc, count_rate=tps)
                            t2(15) = t2(15) + t1 - t0
                            t2(16) = t2(16) + dble(toc - tic)/dble(tps)
                        end if
                    end do
                end if
                if (scatMiet(j)) then
                    do k = 1, nscat*lmax
                        do n = 1, lmax
                            TIJ(js + n:js + n, k) = work(n, n, j)*TIJ(js + n:js + n, k)   !Eq. 48
                        end do
                    end do
                else
                    TIJ(js + 1:js + lmax, :) = &
                        matmul(work(:, :, j), TIJ(js + 1:js + lmax, :))  !Eq. 48
                end if
            end do
            !
            deallocate (work)
            !
            if (cJpresent) then
                do i = 1, size(cJ_, 3)
                    cJ_(:, 2:2, i) = matmul(TIJ, cJ_(:, 1:1, i))   !Eq. 46
                end do
            end if
            !
        elseif (scheme == 2) then
            !write(*,*)'Scheme is=', scheme
            ! Determine collective T-matrices using Stout's scheme
            !
            ! Solve Ax=b by inverting A recursively, using Stout's scheme;
            ! and then compute x=inv(A)b.
            !
            if (verb > 1) write (*, '(A,A)') myname, &
                '> Determine the T-matrix using Stout''s solution scheme'
            !
            if (verb > 1) then
                call cpu_time(t0)
                call system_clock(tic)
            end if
            call calcTIJStout( &
                TIJ=TIJ, &
                scatXYZ=geometry(1:3, :), &
                scatMiet=scatMiet, &
                rtr=rtr)
            if (verb > 1) then
                call cpu_time(t1)
                call system_clock(toc, count_rate=tps)
                t2(17) = t2(17) + t1 - t0
                t2(18) = t2(18) + dble(toc - tic)/dble(tps)
            end if
            if (cJpresent) then
                do i = 1, size(cJ_, 3)
                    cJ_(:, 2:2, i) = matmul(TIJ, cJ_(:, 1:1, i))
                end do
            end if
            !
        elseif (scheme == 3) then
            !
            ! Determine collective T-matrices using Mackowski's scheme
            !
            ! Solve AX=B with multiple right-hand sides, where X contains
            ! Mackowski's TJ matrices and B the offset T1 matrices.
            !
            if (verb > 1) write (*, '(A,A)') myname, &
                '> Determine the T-matrix using Mackowski''s solution scheme'
            !
            if (verb > 1) then
                call cpu_time(t0)
                call system_clock(tic)
            end if
            call calcTIMackowski( &
                TIJ=TIJ, &
                scatXYZ=geometry(1:3, :), &
                scatMiet=scatMiet, &
                rtr=rtr)
            !
            if (verb > 1) then
                call cpu_time(t1)
                call system_clock(toc, count_rate=tps)
                t2(19) = t2(19) + t1 - t0
                t2(20) = t2(20) + dble(toc - tic)/dble(tps)
            end if
            if (cJpresent) then
                allocate (work(lmax, 1, 1))

                do i = 1, size(cJ_, 3)
                    if (PWinc) then
                        call calcJCoeffsPW( &
                            ipwE0=ipwE0(1:3, i), &
                            kVec=hostK_G*ipwDirn(1:3, i), &
                            ipwCoeffsJ=work(:, 1, 1), &
                            xyz=reshape((/0.0d0, 0.0d0, 0.0d0/), (/3, 1/)))  !incidence coefficients in center
                    else
                        !call calcCoeffsJ4dipole( &
                        ! !amp_d, &
                        !  xyz_d = dipoles (1:3, i), &
                        ! idCoeffsJ = work(:,1,1), &
                        ! xyz = reshape((/0.0d0,0.0d0,0.0d0/),(/3,1/)), & !incidence coefficients in center)
                        ! ehost = ehost, &
                        !  k = cmplx(hostK_G,0,kind(hostK_G)))

                    end if
                    do j = 1, nscat
                        js = (j - 1)*lmax
                        cJ_(js + 1:js + lmax, 2, i) = &
                            matmul(TIJ(js + 1:js + lmax, 1:lmax), work(:, 1, 1))
                    end do
                end do
                deallocate (work)
            end if
            !
        else ! Unrecognised scheme
            !
            write (*, '(A,A,i2)') &
                myName, '> ERROR: Unrecognised scheme ', scheme
            if (present(ierr_)) ierr_ = 2
            return
            !
        end if
        !
        if (cJpresent .and. present(cJint_)) then
            !
            cJint_ = 0
            !
            if (any(scatMiet) .and. verb > 0 .and. present(csAbs_)) then
                write (*, '(A,A)') myname, '> Partial absorptions:'
                write (*, '(8x,A,7x,4(A,11x))') 'J', 'k = 0', 'k = 1', 'k = 2', 'k = 3'
            end if
            !
            do j = 1, nscat
                if (scatMiet(j)) then
                    !if(verb > 0) write(*,'(6x,i3)',advance='no') j
                    js = (j - 1)*lmax
                    i = nint(abs(geometry(8, j))) ! shell count
                    ! Again recompute Eric and Pablo's a and k (with reverse
                    ! indexing) from geometry, hostK, and scatK
                    !x(1) = cmplx(geometry(4,j)*hostK_G,0,8)! a*kOut for outer interface
                    !s(1) = scatK(1,j) / hostK_G          ! kIn/kOut for outer interface
                    s(0) = hostK_G
                    s(1:1 + i) = scatK(1:1 + i, j)
                    !if(i > 0) then
                    !x(2:1+i) = geometry(5:4+i,j)*scatK(1:i,j) ! a*kOut
                    !s(2:1+i) = scatK(2:1+i,j)/scatK(1:i,j)  ! kIn/kOut
                    !endif
                    inc(1:4) = 0 ! was used for incidence, now just dummy
                    if (verb > 1) then
                        call cpu_time(t0)
                        call system_clock(tic)
                    end if
                    call calcMieIntCoeffs( &
                        a=geometry(4 + i:4:-1, j), & ! reversed interface indexing
                        k=s(1 + i:0:-1), & ! include hostK_G in the zero-index
                        scaCoeffs=cJ_(js + 1:js + lmax, 2, size(cJ_, 3)), &
                        intCoeffsReg=cJint_(js + 1:js + lmax, 1 + i:1:-1, 1), & ! rev. index
                        intCoeffsIrr=cJint_(js + 1:js + lmax, 1 + i:1:-1, 2), &
                        csAbs=inc(1 + i:1:-1)) ! rev. index
                    !cJint_(js+1:js+lmax,1+i:1:-1,1:2) = cJint_(js+1:js+lmax,1:1+i,1:2)
                    if (verb > 1) then
                        call cpu_time(t1)
                        call system_clock(toc, count_rate=tps)
                        t2(21) = t2(21) + t1 - t0
                        t2(22) = t2(22) + dble(toc - tic)/dble(tps)
                    end if
                    if (present(csAbs_)) csAbs_(j, :) = inc(:)
                    if (verb > 0 .and. present(csAbs_)) then
                        do k = 1 + i, 1, -1
                            write (*, '(1x,es15.8E2)', advance='no') inc(k)
                        end do
                        write (*, *)
                    end if
                end if
            end do
            ! call dumpMatrix(cJint_(:,:,1), 'cJintReg')
            ! call dumpMatrix(cJint_(:,:,2), 'cJintIrr')

            if (dumpScaCoeff_G) then
                call dumpMatrix(cJ_(:, 2, :), 'Sca_Coeff')
            end if
            if (dumpIncCoeff_G) then
                call dumpMatrix(cJ_(:, 1, :), 'Inc_Coeff')
            end if
            !
        end if
        !
        ! <<< End of staging and solving/inverting <<<<<<<<<<<<<<<<<<<<<<<
        !
        if (verb > 0) write (*, '(A,A,/)') myname, '> Done!'
        !
        if (verb > 1) then
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [calcSphBessels (reg. & irreg)] (CPU & real in s): ', &
                t2(1), t2(2)
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [calcMieTMat] (CPU & real in s): ', &
                t2(3), t2(4)
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [calcWignerBigD] (CPU & real in s): ', &
                t2(5), t2(6)
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [stageAmat] (CPU & real in s): ', &
                t2(7), t2(8)
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [balanceVecJ] (CPU & real in s): ', &
                t2(9), t2(10)
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [solLinSys] (CPU & real in s): ', &
                t2(11), t2(12)
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [invSqrMat] (CPU & real in s): ', &
                t2(13), t2(14)
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [balanceMatJI] (CPU & real in s): ', &
                t2(15), t2(16)
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [calcTIJStout] (CPU & real in s): ', &
                t2(17), t2(18)
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [calcTIMackowski] (CPU & real in s): ', &
                t2(19), t2(20)
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [calcMieIntCoeffs] (CPU & real in s): ', &
                t2(21), t2(22)
        end if
        !
        deallocate (bes4bal_G)
        !
    end subroutine solve
    !
    !==================================================================
    !
    subroutine stageAmat(A, scatXYZ, scatMiet, rtr, right_, Tmats_, balance_, verb_)
        !
        ! Stage a prestaged matrix A, which must contain 1-body
        ! T-matrices in the diagonal blocks on input.
        !
        use swav, only: calcVTACS
        !
        !---------------------------------------------------
        ! Start of variable declarations
        !---------------------------------------------------
        ! Passed arguments:
        complex(8), intent(inout) :: A(:, :)
        real(8), intent(in) :: scatXYZ(:, :)
        logical, intent(in) :: scatMiet(size(scatXYZ, 2)), rtr
        logical, intent(in), optional :: right_, balance_
        complex(8), intent(out), optional :: Tmats_( &
                                             size(A, 1)/size(scatXYZ, 2), size(A, 2)/size(scatXYZ, 2), size(scatXYZ, 2))
        integer, intent(in), optional :: verb_
        ! Local variables:
        character(*), parameter :: myName = 'stageAmat'
        logical :: balance, stage, right
        integer :: j, k, js, ks, n, i, lmax, nscat, verb !,ij(2)
        real(8) :: kr(3)
        complex(8) :: off(3, 1)
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        if (size(A, 1) /= size(A, 2)) then
            write (*, '(A,A)') myName, '> ERROR: size(A,1) /= size(A,2)'
            STOP
        else if (size(scatXYZ, 1) /= 3) then
            write (*, '(A,A)') myName, '> ERROR: size(scatXYZ,1) /= 3'
            STOP
        end if
        !
        nscat = size(scatXYZ, 2)
        lmax = size(A, 1)/nscat
        !
        if (size(A, 1) /= nscat*lmax) then
            write (*, '(A,A)') myName, '> ERROR: size(A,1) /= nscat*lmax'
            STOP
        end if
        !
        verb = 0
        if (present(verb_)) verb = verb_
        !
        if (present(right_) .and. present(Tmats_)) then
            stage = .true.
            right = right_
            !write(*,*) 'TEST: stage, right=', stage, right
        elseif (present(right_) .or. present(Tmats_)) then
            write (*, '(A,A)') myName, &
                '> ERROR: Only one of right_ and Tmats_ present'
            STOP
        end if
        !
        if (present(balance_)) then
            balance = balance_
        else
            balance = .false.
        end if
        !
        columns: do k = 1, nscat
            ks = (k - 1)*lmax
            rows: do j = 1, nscat
                js = (j - 1)*lmax
                !
                if (j == k) cycle rows
                !
                if (stage) then
                    if (verb > 1) write (*, *)
                    if (right) then
                        A(js + 1:js + lmax, ks + 1:ks + lmax) = -A(js + 1:js + lmax, js + 1:js + lmax)
                    else
                        A(js + 1:js + lmax, ks + 1:ks + lmax) = -A(ks + 1:ks + lmax, ks + 1:ks + lmax)
                    end if
                    kr(1:3) = hostk_G*(scatXYZ(1:3, j) - scatXYZ(1:3, k))
                    !
                    do i = 1, 3
                        off(i, 1) = cmplx(0, kr(i), 8) ! irregular variant
                    end do
                    ! if(balance) then
                    !    ij = (/j,k/)
                    ! else
                    !    ij = 0
                    ! endif
                    call offsetTmat( &
                        Tmat=A(js + 1:js + lmax, ks + 1:ks + lmax), &
                        off=off(:, 1:1), &
                        miet=scatMiet(k), &
                        rtr=rtr, &
                        right=right)!, & !end here if no balancing inside offsetTmat
                    !balJI_ = ij ) ! balancing inside offsetTmat in dev
                    if (balance) then
                        call balanceMatJI( &
                            j=j, jregt=.not. right, &
                            Mat=A(js + 1:js + lmax, ks + 1:ks + lmax), &
                            i=k, iregt=.not. right)
                        if (verb > 1) write (*, '(A,A,2(1x,i2))') myName, &
                            '> Offset and balanced block (j,k)=', j, k
                    elseif (verb > 1) then
                        write (*, '(A,A,2(1x,i3))') myName, &
                            '> Offset block (j,k)=', j, k
                    end if
                    ! No sure if it's worth diagnosing the half-offset T-matrices...
                    ! call diagnoseTmat( &
                    !      Tmat = A(js+1:js+lmax, ks+1:ks+lmax), &
                    !      mode_ = 0, verb_ = verb)
                else
                    call calcVTACs( &
                        r0=scatXYZ(1:3, j) - scatXYZ(1:3, k), &
                        k=cmplx(hostK_G, 0, kind(hostK_G)), &
                        regt=.false., &
                        vtacs=A(js + 1:js + lmax, ks + 1:ks + lmax))
                end if
                !
            end do rows
        end do columns
        !
        if (stage) then
            if (verb > 1) write (*, *)
            do k = 1, nscat
                ks = (k - 1)*lmax
                ! Move (k,k)th block from A to Tmats_ and replace it with identity
                Tmats_(:, :, k) = A(ks + 1:ks + lmax, ks + 1:ks + lmax)
                A(ks + 1:ks + lmax, ks + 1:ks + lmax) = 0
                do n = 1, lmax
                    A(ks + n, ks + n) = 1 ! Note: remains unity when balanced
                end do
                !
            end do
        end if
        !
        if (dumpStagedA_G) then
            if (balance) then
                call dumpMatrix(mat=A, &
                                ofile='stagedA_bal.txt', &
                                verb_=verb)
            else
                call dumpMatrix(mat=A, &
                                ofile='stagedA.txt', &
                                verb_=verb)
            end if
        end if
        !
    end subroutine stageAmat
    !
    !==============================================================
    !
    subroutine calcTIJStout(TIJ, scatXYZ, scatMiet, rtr)
        !
        ! ============================================================
        ! Compute scatterer-centred T-matrices using the recursive
        ! scheme in Stout02 with matrix balancing from Stout08. The
        ! relevant equations are 33 and 35 in Stout02, and 20, 22
        ! and 24 in Stout08.
        ! ============================================================
        !
        use swav, only: xyz2rtp, calcWignerBigD, calcVTACs, calcVTACsAxial
        use linalg, only: invSqrMat
        !
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables
        complex(8), intent(inout) :: TIJ(:, :)
        real(8), intent(in) :: scatXYZ(:, :)
        logical, intent(in) :: scatMiet(size(scatXYZ, 2)), rtr
        ! Local variables
        character(*), parameter :: myName = 'calcTIJStout'
        integer :: j, k, N, js, ks, Ns, tic, toc, tps, nscat, lmax, pmax, nmax !,i
        logical, dimension(size(scatXYZ, 2) - 1) :: applyRot
        logical :: flip
        real(8) :: dr(3), rtp(3)
        real    :: t1, t2
        complex(8), dimension(size(TIJ, 1)/size(scatXYZ, 2), &
                              size(TIJ, 2)/size(scatXYZ, 2)) :: matSum, mat
        complex(8), dimension(size(TIJ, 1)/size(scatXYZ, 2), &
                              size(TIJ, 2)/size(scatXYZ, 2), size(scatXYZ, 2) - 1) :: &
            vcoeff_mN, vcoeff_Nm, matSumJ
        complex(8), dimension(size(TIJ, 1)/(2*size(scatXYZ, 2)), &
                              size(TIJ, 2)/(2*size(scatXYZ, 2)), size(scatXYZ, 2) - 1) :: bigD
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        if (size(TIJ, 1) /= size(TIJ, 2)) then
            write (*, '(A,A)') myName, '> ERROR: size(TIJ,1) /= size(TIJ,2)'
            STOP
        elseif (size(scatXYZ, 1) /= 3) then
            write (*, '(A,A)') myName, '> ERROR: size(scatXYZ,1) /= 3'
            STOP
        end if
        !
        nscat = size(scatXYZ, 2)
        lmax = size(TIJ, 1)/nscat ! was M
        pmax = lmax/2
        nmax = int(sqrt(dble(pmax)))
        if (size(TIJ, 1) /= 2*pmax*nscat) then
            write (*, '(A,A)') myName, '> ERROR: size(TIJ,1) /= 2*pmax*nscat'
            STOP
        end if
        !
        ! Recursively apply equations 33 and 35 for N > 1.
        Nloop: do N = 2, nscat
            !
            Ns = lmax*(N - 1)
            !
            ! --- Pre-compute all the offsetting matrices in first jLoop ---------
            jLoop1: do j = 1, N - 1
                !
                dr(:) = scatXYZ(1:3, j) - scatXYZ(1:3, N)
                if (rtr) then
                    call xyz2rtp(dr, rtp)
                    if ((abs(rtp(2)) < tiny1 .or. abs(rtp(2) - pi) < tiny1) .and. &
                        (abs(rtp(3)) < tiny1 .or. abs(rtp(3) - tpi) < tiny1)) then
                        applyRot(j) = .false.
                        ! In the absence of rotation, must check the sign of dr(3)
                        ! and assign the axial vtacs_mN and vtacs_Nm correctly, to
                        ! make sure the intended axial translation is consistent
                        ! with the z-axis direction.
                        if (dr(3) < 0) then
                            flip = .true.
                        else
                            flip = .false.
                        end if
                    else
                        applyRot(j) = .true.
                        ! With rotation, the sign of dr(3) is irrelevant since the
                        ! rotation will always align the local axes so that dr points
                        ! in the +ve z-direction, so no need to flip!
                        flip = .false.
                        call calcWignerBigD( &
                            angles=(/rtp(3), rtp(2), 0.0d0/), &
                            pmax=pmax, &
                            bigD=bigD(:, :, j))
                    end if
                    call calcVTACsAxial( &
                        r0=sqrt(dot_product(dr, dr)), &
                        k=cmplx(hostK_G, 0, kind(hostK_G)), &
                        pmax=pmax, &
                        regt=.false., &
                        flip=flip, &
                        vtacs=vcoeff_mN(:, :, j))
                    call calcVTACsAxial( &
                        r0=sqrt(dot_product(dr, dr)), &
                        k=cmplx(hostK_G, 0, kind(hostK_G)), &
                        pmax=pmax, &
                        regt=.false., &
                        flip=.not. flip, &
                        vtacs=vcoeff_Nm(:, :, j))
                else
                    call calcVTACs( &
                        r0=dr, &
                        k=cmplx(hostK_G, 0, kind(hostK_G)), &
                        regt=.false., &
                        vtacs=vcoeff_mN(:, :, j))
                    call calcVTACs( &
                        r0=-dr, &
                        k=cmplx(hostK_G, 0, kind(hostK_G)), &
                        regt=.false., &
                        vtacs=vcoeff_Nm(:, :, j))
                end if
            end do jLoop1
            !
            ! --- Use equation 33a to get T_N^(N,N) -------------------
            !
            matSum(:, :) = 0 ! Initialise sum
            !
            jLoop2: do j = 1, N - 1
                js = lmax*(j - 1)
                matSumJ(:, :, j) = 0 ! Initialise recurring sum
                !
                do k = 1, N - 1
                    ks = lmax*(k - 1)
                    !
                    mat(1:lmax, 1:lmax) = TIJ(js + 1:js + lmax, ks + 1:ks + lmax)
                    !
                    call offsetTmat( &
                        Tmat=mat, &
                        off=vcoeff_mN(:, :, k), &
                        miet=.false., &
                        rtr=rtr, &
                        right=.true., &
                        bigD_=bigD(:, :, k), &
                        useD_=applyRot(k))
                    !
                    matSumJ(:, :, j) = matSumJ(:, :, j) + mat(:, :)
                    !
                end do
                !
                mat(:, :) = matSumJ(:, :, j)
                !
                call offsetTmat( &
                    Tmat=mat, &
                    off=vcoeff_Nm(:, :, j), &
                    miet=.false., &
                    rtr=rtr, &
                    right=.false., &
                    bigD_=bigD(:, :, j), &
                    useD_=applyRot(j))
                !
                matSum = matSum + mat
                !
            end do jLoop2
            !
            if (scatMiet(N)) then ! Right-multiply by diagonal tmat1
                do j = 1, lmax
                    matSum(:, j) = -matSum(:, j)*TIJ(Ns + j, Ns + j)
                end do
            else ! general mult
                matSum = -matmul(matSum, TIJ(Ns + 1:Ns + lmax, Ns + 1:Ns + lmax))
            end if
            !
            do j = 1, lmax ! Compute the bracketed term to be inverted
                matSum(j, j) = 1 + matSum(j, j)
            end do
            !
            ! Balance matSum = inv(T_N^(N,N)) before inverting
            if (balStout_G) then
                call balanceMatJI( &
                    j=N, jregt=.true., Mat=matSum, iregt=.true., i=N, &
                    rev_=.false.)
            end if
            !
            ! Invert the bracketed term
            call invSqrMat( &
                A=matSum, &
                trans_=transInv_G, &
                verb_=verb_G)
            !
            ! Unbalance the inverted matSum = T_N^(N,N)
            if (balStout_G) then
                call balanceMatJI( &
                    j=N, jregt=.true., Mat=matSum, iregt=.true., i=N, &
                    rev_=.true.)
            end if
            !
            ! Store T_N^(N,N)
            if (scatMiet(N)) then ! left-multiply by diagonal tmat1
                matSum = transpose(matsum)
                do j = 1, lmax
                    matSum(:, j) = matSum(:, j)*TIJ(Ns + j, Ns + j) ! righ-mult after transpose
                end do
                TIJ(Ns + 1:Ns + lmax, Ns + 1:Ns + lmax) = transpose(matSum)
            else
                TIJ(Ns + 1:Ns + lmax, Ns + 1:Ns + lmax) = &
                    matmul(TIJ(Ns + 1:Ns + lmax, Ns + 1:Ns + lmax), matSum)
            end if
            !
            ! --- Use equation 33b to get T_N^(N,k) for k<N -----------
            if (verb_G > 1) then
                call cpu_time(t1)
                call system_clock(tic)
            end if
            do k = 1, N - 1
                ks = lmax*(k - 1)
                matSum = 0
                do j = 1, N - 1
                    js = lmax*(j - 1)
                    !
                    mat = TIJ(js + 1:js + lmax, ks + 1:ks + lmax)
                    !
                    call offsetTmat( &
                        Tmat=mat, &
                        off=vcoeff_Nm(:, :, j), &
                        miet=.false., &
                        rtr=rtr, &
                        right=.false., &
                        bigD_=bigD(:, :, j), &
                        useD_=applyRot(j))
                    !
                    matSum = matSum + mat
                    !
                end do
                !
                TIJ(Ns + 1:Ns + lmax, ks + 1:ks + lmax) = matmul( &       ! matmul bottleneck
                                                          TIJ(Ns + 1:Ns + lmax, Ns + 1:Ns + lmax), matSum)
                !
            end do
            !
            ! --- Now update T_N^(j,N) and T_N^(j,k) for j<N and k<N --------
            do j = 1, N - 1
                js = lmax*(j - 1)
                !
                ! --- Use equation 35a to get T_N^(j,N) for j<N --------------
                matSum(:, :) = matSumJ(:, :, j)
                ! Store T_N^(j,N)
                TIJ(js + 1:js + lmax, Ns + 1:Ns + lmax) = matmul( &       ! matmul bottleneck
                                                          matSum, TIJ(Ns + 1:Ns + lmax, Ns + 1:Ns + lmax))
                !
                ! --- Use equation 35b to update T_N^(j,k) for j<N and k<N ---
                do k = 1, N - 1
                    ks = lmax*(k - 1)
                    !
                    ! Store T_N^(j,k)
                    ! Matmul is the killer, and zgemm isn't much better.
                    TIJ(js + 1:js + lmax, ks + 1:ks + lmax) = &
                        TIJ(js + 1:js + lmax, ks + 1:ks + lmax) + matmul( &
                        matSum, TIJ(Ns + 1:Ns + lmax, ks + 1:ks + lmax))
                    ! call zgemm('N','N', lmax, lmax, lmax, cmplx(1,0,8), &
                    !      matSum,M,TIJ(Ns+1:Ns+lmax, ks+1:ks+lmax),lmax,cmplx(1,0,8), &
                    !      TIJ(js+1:js+lmax,ks+1:ks+lmax), lmax)
                end do
                !
            end do
            if (verb_G > 1) then
                call system_clock(toc, tps)
                call cpu_time(t2)
                write (*, '(A,A,2(1x,es10.3E2))') myName, &
                    '> Calculation time [post-inversion] (CPU & real in s): ', t2 - t1, &
                    dble(toc - tic)/dble(tps)
            end if
            !
        end do Nloop
        !
        ! if(verb_G > 0) write(*,'(A,A,2(1x,es9.2e3))') &
        !      myName,'> min/maxval(abs(TIJ))=', &
        !      minval(abs(TIJ)), maxval(abs(TIJ))
        !
        ! Testing>
        ! write(*,*) myname,'> Final TIJ:'
        ! do i=1,3
        !    write(*,*) (TIJ(i,j), j=1,3)
        ! enddo
        ! < end of Testing.
        !
        call flush ()
        !
    end subroutine calcTIJStout
    !
    subroutine calcTIMackowski(TIJ, scatXYZ, scatMiet, rtr)
        !
        use linalg, only: solLinSys
        !
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables
        complex(8), intent(inout) :: TIJ(:, :)
        real(8), intent(in) :: scatXYZ(:, :)
        logical, intent(in) :: scatMiet(size(scatXYZ, 2)), rtr
        ! Local variables
        character(*), parameter :: myName = 'calcTIMackowski'
        integer :: k, ks, j, nscat, lmax
        complex(8) :: off(3, 1)
        complex(8), dimension(size(TIJ, 1), size(TIJ, 2)/size(scatXYZ, 2)) :: B
        complex(8), allocatable :: X(:, :, :)
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        if (size(TIJ, 1) /= size(TIJ, 2)) then
            write (*, '(A,A)') myName, '> ERROR: size(TIJ,1) /= size(TIJ,2)'
            STOP
        elseif (size(scatXYZ, 1) /= 3) then
            write (*, '(A,A)') myName, '> ERROR: size(scatXYZ,1) /= 3'
            STOP
        end if
        !
        nscat = size(scatXYZ, 2)
        lmax = size(TIJ, 1)/nscat ! was M
        if (size(TIJ, 1) /= lmax*nscat) then
            write (*, '(A,A)') myName, '> ERROR: size(TIJ,1) /= lmax*nscat'
            STOP
        end if
        !
        allocate (X(lmax, lmax, nscat))
        call stageAmat( &
            A=TIJ, &
            scatXYZ=scatXYZ, &
            rtr=rtr, &
            scatMiet=scatMiet, &
            right_=.true., &
            Tmats_=X, &
            balance_=balStout_G, &
            verb_=verb_G)
        columns: do k = 1, nscat
            ks = (k - 1)*lmax
            off(1:3, 1) = cmplx(hostK_G, 0, 8)*scatXYZ(1:3, k)  ! regular
            call offsetTmat( &
                Tmat=X(:, :, k), &
                off=off, &
                miet=scatMiet(k), &
                rtr=rtr, &
                right=.true.)
            if (balStout_G) then
                ! Balance each column vector of X in accordance with the staged A
                do j = 1, lmax
                    call balanceVecJ( &
                        j=k, jregt=.false., &
                        Vec=X(:, j, k), &
                        rev_=.false.)
                end do
                ! Balance columns and row of X
                ! call balanceMatJI( &
                !      j = k, jregt = .false., &
                !      Mat = X(:,:,k), &
                !      i = k, iregt = .true., &
                !      rev_ = .false. )
            end if
            B(ks + 1:ks + lmax, :) = X(:, :, k)
        end do columns
        deallocate (X)
        allocate (X(size(B, 1), size(B, 2), 1))
        !
        ! Solve AX=B for X using LAPACK
        X(:, :, 1) = B
        call solLinSys( &
            A=TIJ, & ! Note: TIJ is changed!
            X=X(:, :, 1), &
            isol_=isolve_G, &
            verb_=verb_G)
        !
        TIJ = 0
        if (balStout_G) then
            do k = 1, nscat
                ks = (k - 1)*lmax
                ! Unbalance each column vector of X
                do j = 1, lmax
                    call balanceVecJ( &
                        j=k, jregt=.false., &
                        Vec=X(ks + 1:ks + lmax, j, 1), &
                        rev_=.true.)
                end do
                ! Unbalance columns and rows of X
                ! call balanceMatJI( &
                !         j = k, jregt = .false., &
                !         Mat = X(ks+1:ks+lmax,1:lmax,1), &
                !         i = k, iregt = .true., &
                !         rev_ = .true. )
            end do
        end if
        TIJ(:, 1:lmax) = X(:, 1:lmax, 1)
        TIJ(:, lmax + 1:) = 0
        !
        deallocate (X)
        !
    end subroutine calcTIMackowski
    !
    !==============================================================
    !
    subroutine balanceMatJI(j, jregt, Mat, iregt, i, rev_, mnq_)
        !
        !============================================================
        ! "Balance" a single matrix (M) using two weights (pointed
        ! at by j and i). M is supposed to relate two vectors of VSW
        ! coefficients, cj (centred at j) and ci (centred at i), such
        ! that cj = M ci. Logicals jregt and iregt specify whether
        ! cj and ci are regular or not. The optional rev_ is .false.
        ! by default, but making it .true. will trigger the reverse
        ! of balancing - "unbalancing". The optinal mnq_ is .false.
        ! by default, but making it true will change the indexing
        ! convention from (q,n,m) to (m,n,q), which is used to make
        ! the z-axial VTACs block-diagonal.
        !============================================================
        !
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variable
        complex(8), intent(inout) :: Mat(:, :)
        integer, intent(in) :: i, j ! weight pointers
        logical, intent(in) :: iregt, jregt ! regularity indicators
        logical, intent(in), optional :: rev_, mnq_
        ! Local variables
        character(*), parameter :: myName = 'balanceMatJI'
        logical :: rev, mnq
        integer :: l, lp, n, m, np, mp, q, qp, ireg, jreg, lmax, nmax
        complex(8) :: zbal
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        lmax = size(Mat, 1)
        if (size(Mat, 2) /= lmax) then
            write (*, '(A,A)') myName, '> ERROR: size(Mat,1) /= size(Mat,2)'
            STOP
        else
            nmax = int(sqrt(dble(lmax/2)))
            if (lmax /= 2*nmax*(nmax + 2)) then
                write (*, '(A,A)') myName, '> ERROR: lmax /= 2*nmax*(nmax+2)'
            end if
        end if
        !
        if (present(rev_)) then
            rev = rev_
        else
            rev = .false.
        end if
        !
        if (present(mnq_)) then
            mnq = mnq_
        else
            mnq = .false.
        end if
        !
        if (iregt) then
            ireg = 1
        else
            ireg = 0
        end if
        !
        if (jregt) then
            jreg = 1
        else
            jreg = 0
        end if
        !
        if (mnq) then
            !
            ! Non-standard (m,n,q) indexing
            lp = 0
            do mp = -nmax, nmax
                do np = max(1, abs(mp)), nmax
                    do qp = 1, 2
                        lp = lp + 1
                        !
                        l = 0
                        do m = -nmax, nmax
                            do n = max(1, abs(m)), nmax
                                !
                                if (rev) then
                                    zbal = bes4bal_G(np, ireg, i)/bes4bal_G(n, jreg, j)
                                else
                                    zbal = bes4bal_G(n, jreg, j)/bes4bal_G(np, ireg, i)
                                end if
                                !
                                do q = 1, 2
                                    l = l + 1
                                    !
                                    Mat(l, lp) = zbal*Mat(l, lp)
                                    !
                                end do
                            end do
                        end do
                        !
                    end do
                end do
            end do

            !
        else
            ! Default (q,n,m) indexing
            lp = 0
            do qp = 1, 2
                do np = 1, nmax
                    do mp = -np, np
                        lp = lp + 1
                        !
                        l = 0
                        do q = 1, 2
                            do n = 1, nmax
                                !
                                if (rev) then
                                    zbal = bes4bal_G(np, ireg, i)/bes4bal_G(n, jreg, j)
                                else
                                    zbal = bes4bal_G(n, jreg, j)/bes4bal_G(np, ireg, i)
                                end if
                                !
                                do m = -n, n
                                    l = l + 1
                                    !
                                    Mat(l, lp) = zbal*Mat(l, lp)
                                    !
                                end do
                            end do
                        end do
                        !
                    end do
                end do
            end do
            !
        end if
        !
    end subroutine balanceMatJI
    !
    !===============================================================
    !
    subroutine balanceVecJ(j, jregt, Vec, rev_)
        !
        !============================================================
        ! "Balance" a single vector (V) using two weights (pointed
        ! at by j). V is supposed to contain VSW coefficients centred
        ! at j. Logical jregt specifies whether the VSWs are  regular
        ! or not. The optional rev_ is .false. by default, but making
        ! it .true. triggers the reverse of balancing - "unbalancing"
        !============================================================
        !
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variable
        complex(8), intent(inout) :: Vec(:)
        integer, intent(in) :: j ! weight pointers
        logical, intent(in) :: jregt ! regularity indicators
        logical, intent(in), optional :: rev_
        ! Local variables
        character(*), parameter :: myName = 'balanceVecJ'
        logical :: rev
        integer :: l, n, m, q, jreg, lmax, nmax
        complex(8) :: zbal
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        lmax = size(Vec)
        nmax = int(sqrt(dble(lmax/2)))
        if (lmax /= 2*nmax*(nmax + 2)) then
            write (*, '(A,A)') myName, '> ERROR: lmax /= 2*nmax*(nmax+2)'
        end if
        !
        if (present(rev_)) then
            rev = rev_
        else
            rev = .false.
        end if
        !
        if (jregt) then
            jreg = 1
        else
            jreg = 0
        end if
        !
        l = 0
        do q = 1, 2
            do n = 1, nmax
                !
                if (rev) then
                    zbal = 1/bes4bal_G(n, jreg, j)
                else
                    zbal = bes4bal_G(n, jreg, j)
                end if
                !
                do m = -n, n
                    l = l + 1
                    !
                    Vec(l) = zbal*Vec(l)
                    !
                end do
            end do
        end do
        !
        return
        !
    end subroutine balanceVecJ
    !
    !===============================================================
    !
    subroutine calcCsStout(scatXYZR, aJ, fJ, sig, nmax2_, tol_, verb_, vtacs_t)
        !
        !============================================================
        ! Calculate the extinction, scattering, and absorption cross-
        ! sections from the incident and scattered coefficients.
        ! Depending on the shape of the array sig, each cross section
        ! is either just a sum total, or resolved into contributions
        ! from the multipole orders.
        !============================================================
        !
        use swav, only: calcVTACs
        !
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables:
        real(8), intent(in) :: scatXYZR(:, :)
        complex(8), intent(in) :: aJ(:), fJ(size(aJ))
        complex(8), intent(in) :: vtacs_t(:, :)
        real(8), intent(inout) :: sig(:, :) ! 3 by 1 or 3 by nmax
        real(8), intent(in), optional :: tol_
        integer, intent(in), optional :: verb_, nmax2_
        ! Local variables
        character(*), parameter :: myName = 'calcCsStout'
        integer :: i, j, k, n, ks, js, is, ps
        integer :: verb, lmax, pmax, nmax, nscat, nmax2, pmax2, lmax2
        real(8) :: hostK, hostK2, rdum, tol
        real(8), dimension(size(scatXYZR, 2)) :: sigExtJ, sigScaJ, sigAbsJ
        complex(8) :: zdummy!, zdum(size(aJ)/size(scatXYZR,2),1)
        complex(8), allocatable :: vtacs(:, :), my_fJ(:), zdum(:)
        integer :: nconvExt(size(scatXYZR, 2)), nconvSca(size(scatXYZR, 2))
        logical :: multisplit

        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------

        nscat = size(scatXYZR, 2)
        lmax = size(aJ)/nscat
        pmax = lmax/2
        nmax = int(sqrt(dble(pmax)))
        if (size(aJ) /= nscat*lmax) then
            write (*, '(A,A)') myname, '> ERROR: size(aJ) /= nscat*lmax'
            STOP
        elseif (lmax /= 2*nmax*(nmax + 2)) then
            write (*, '(A,A)') myname, '> ERROR: lmax /= 2*nmax*(nmax+2)'
            STOP
        elseif (size(scatXYZR, 1) /= 4) then
            write (*, '(A,A)') myname, '> ERROR: size(scatXYZR,1) /= 3'
            STOP
        end if
        !
        nmax2 = nmax
        if (present(nmax2_)) nmax2 = nmax2_
        if (nmax2 < nmax) then
            write (*, '(A,A)') myname, '> ERROR: nmax2 < nmax'
            write (*, *) nmax2, nmax
            STOP
        end if
        pmax2 = nmax2*(nmax2 + 2)
        lmax2 = 2*pmax2
        !
        if (size(sig, 1) /= 3) then
            write (*, '(A,A)') myname, '> ERROR: size(sig,1) /= 3'
            STOP
        elseif (size(sig, 2) == 1) then
            multisplit = .false.
        elseif (size(sig, 2) == nmax2) then
            multisplit = .true.
        else
            write (*, '(A,A)') myname, '> ERROR: Bad size(sig,2)'
            STOP
        end if
        !
        allocate (vtacs(lmax2, lmax2), my_fJ(nscat*lmax2), zdum(nscat*lmax2))
        !
        tol = 0
        if (present(tol_)) then
            if (tol_ >= 1) then
                write (*, '(A,A)') myname, '> ERROR: tol >= 1'
                STOP
            elseif (tol_ >= 0) then
                tol = tol_
            else
                write (*, '(A,A)') myname, '> ERROR: tol < 0'
                STOP
            end if
        end if
        !
        verb = 0
        if (present(verb_)) verb = verb_
        !
        hostK = sig(1, 1)
        hostK2 = hostK**2 ! host wavenumber squared
        !
        ! Pad fJ with zeroes...
        do i = 1, nscat
            is = (i - 1)*lmax
            js = (i - 1)*lmax2
            my_fJ(js + 1:js + lmax) = fJ(is + 1:is + lmax)
            my_fJ(js + lmax + 1:js + lmax2) = 0
        end do
        !
        sig = 0
        !
        ! Compute the extinction cross-section:
        ! >>> old way without convergence check
        ! jscat1: do j=1,nscat
        !    js=(j-1)*lmax
        !    sigExtJ(j) = -realpart(dot_product(aJ(js+1:js+lmax),fJ(js+1:js+lmax)))/hostK2
        !    ! Note that DOT_PRODUCT=SUM(CONJG(VECTOR_A)*VECTOR_B) for complex vectors
        !    ! It can trigger IEEE_UNDERFLOW_FLAG and IEEE_DENORMAL compiler warnings
        !    ! here, if nmax is very high and high-n contributions become nearly zero.
        !    sig(1) = sig(1) + sigExtJ(j)
        ! enddo jscat1
        ! <<< old way
        sigExtJ = 0
        nconvExt = 0
        jscat1: do j = 1, nscat
            js = (j - 1)*lmax
            nloop1: do n = 1, nmax
                if (nconvExt(j) > 0) exit nloop1 ! is this exit important?
                ps = n*(n + 1)
                zdummy = dot_product(aJ(js + ps - n:js + ps + n), fJ(js + ps - n:js + ps + n))
                ps = ps + pmax
                zdummy = zdummy + &
                         dot_product(aJ(js + ps - n:js + ps + n), fJ(js + ps - n:js + ps + n))
                rdum = realpart(zdummy)
                sigExtJ(j) = sigExtJ(j) + rdum
                if (abs(rdum) > 0 .and. abs(sigExtJ(j)) > 0) then
                    if (abs(rdum/sigExtJ(j)) < tol) nconvExt(j) = n
                end if
                if (multisplit) then
                    sig(1, n) = sig(1, n) + rdum
                else
                    sig(1, 1) = sig(1, 1) + rdum
                end if
            end do nloop1
        end do jscat1
        sigExtJ = -sigExtJ/hostK2
        sig(1, :) = -sig(1, :)/hostK2
        !sig(1) = sum(sigExtJ)
        !
        ! Compute the scattering cross-sections
        sigScaJ = 0
        nconvSca = 0
        jscat2: do j = 1, nscat
            js = (j - 1)*lmax2
            zdum = 0
            kscat: do k = 1, nscat
                ks = (k - 1)*lmax2
                if (k /= j) then
                    vtacs = vtacs_t(js + 1:js + lmax, ks + 1:ks + lmax)

                    !  call calcVTACs( &
                    !       r0 = scatXYZR(1:3,j)-scatXYZR(1:3,k), &
                    !       k = cmplx(hostK,0,kind(hostK)), &
                    !      regt = .true., &
                    !      vtacs = vtacs )

                    zdum(1:lmax2) = zdum(1:lmax2) + &
                                    matmul(vtacs(1:lmax2, 1:lmax2), my_fJ(ks + 1:ks + lmax2))

                else
                    zdum(1:lmax2) = zdum(1:lmax2) + my_fJ(ks + 1:ks + lmax2)
                end if
            end do kscat
            ! old way, before convergence testing for scattering
            !sigScaJ(j) = dot_product(fJ(js+1:js+lmax),zdum(1:lmax,1))/hostK2
            !
            ! if( abs(imagpart(sigScaJ(j)))/abs(sigScaJ(j)) > tiny ) then
            !    write(*,'(A,A,es10.3E2,A,i3)') myname, &
            !         '> WARNING: |Im{sigScaJ}|/|sigScaJ|= ', &
            !         abs(imagpart(sigScaJ(j)))/abs(sigScaJ(j)),' for J= ',j
            ! endif
            !
            nloop2: do n = 1, nmax2
                if (nconvSca(j) > 0) exit nloop2 ! does this help?!?
                ps = n*(n + 1)
                zdummy = dot_product(my_fJ(js + ps - n:js + ps + n), zdum(ps - n:ps + n))
                ps = ps + pmax2
                zdummy = zdummy + &
                         dot_product(my_fJ(js + ps - n:js + ps + n), zdum(ps - n:ps + n))
                rdum = realpart(zdummy)
                sigScaJ(j) = sigScaJ(j) + rdum
                if (abs(rdum) > 0 .and. abs(sigScaJ(j)) > 0) then
                    if (abs(rdum/sigScaJ(j)) < tol) nconvSca(j) = n
                end if
                if (multisplit) then
                    sig(2, n) = sig(2, n) + rdum
                else
                    sig(2, 1) = sig(2, 1) + rdum
                end if

            end do nloop2
        end do jscat2
        sigScaJ = sigScaJ/hostK2
        sig(2, :) = sig(2, :)/hostK2
        !sig(2) = sum(sigScaJ)
        sig(3, :) = sig(1, :) - sig(2, :)
        !
        ! Use energy conservation to infer the absorption cross-section
        if (verb > 0) then
            write (*, '(A,A,3(5x,A,3x),1x,A)') myname, '> J', 'csExt(J)', 'csSca(J)', 'csAbs(J)', 'nConv'
        end if
        !sig(3) = 0
        do j = 1, nscat
            sigAbsJ(j) = sigExtJ(j) - sigScaJ(j)
            !sig(3) = sig(3) + sigAbsJ(j)
            if (verb > 0) then
                write (*, '(13x,i4,3(1x,es15.8e2), 2(1x, i2))') &
                    J, sigExtJ(j), sigScaJ(j), sigAbsJ(j), &
                    nConvExt(j), nConvSca(j)
            end if

        end do
        !
        if (verb > 0) then
            write (*, '(13x,51("-"))')
            write (*, '(13x,A,3(1x,es15.8e2))') 'Sum', sum(sigExtJ(:)), &
                sum(sigScaJ(:)), sum(sigAbsJ(:))
            write (*, *)
        end if
        !
        ! Sanity checks
        if (abs(sum(sigExtJ(:)) - sum(sig(1, :))) > tol) then
            write (*, '(A,A)') myname, 'WARNING: Total extinction inconsistent'
        elseif (abs(sum(sigScaJ(:)) - sum(sig(2, :))) > tol) then
            write (*, '(A,A)') myname, 'WARNING: Total scattering inconsistent'
        elseif (abs(sum(sigAbsJ(:)) - sum(sig(3, :))) > tol) then
            write (*, '(A,A)') myname, 'WARNING: Total absorption inconsistent'
        end if
        !
        deallocate (vtacs, my_fJ, zdum)

    end subroutine calcCsStout
    !
    !===============================================================
    subroutine calcCs(scatXYZR, inc, fJ, sig, nmax2_, tol_, verb_)
        !
        !============================================================
        ! Calculate the extinction, scattering, and absorption cross-
        ! sections from the incident and scattered coefficients.
        ! Depending on the shape of the array sig, each cross section
        ! is either just a sum total, or resolved into contributions
        ! from the multipole orders.
        !============================================================
        !
        use swav, only: calcVTACs, calcCoeffsPW
        !
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables:
        real(8), intent(in) :: scatXYZR(:, :), inc(4)
        complex(8), intent(in) :: fJ(:)
        real(8), intent(inout) :: sig(:, :) ! 3 by 1 or 3 by nmax
        real(8), intent(in), optional :: tol_
        integer, intent(in), optional :: verb_, nmax2_
        ! Local variables
        character(*), parameter :: myName = 'calcCs'
        integer :: i, n, is, is2, tic, toc, tps
        integer :: verb, lmax, pmax, nmax, nscat, lmax2, pmax2, nmax2
        real(8) ::  hostK, hostK2
        real    :: t0, t1, t2
        real(8) :: csExt, csSca, csExt_n, csSca_n, ipwDirn(3), tol, t3
        complex(8) :: ipwAmpl(3)
        complex(8), allocatable :: my_fJ(:), a(:), f(:), vtacs(:, :)
        integer :: nconvExt, nconvSca
        logical :: multisplit
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        t2 = 0.0d0
        t3 = 0.0d0
        nscat = size(scatXYZR, 2)
        lmax = size(fJ)/nscat
        pmax = lmax/2
        nmax = int(sqrt(dble(pmax)))
        if (size(fJ) /= nscat*lmax) then
            write (*, '(A,A)') myname, '> ERROR: size(fJ) /= nscat*lmax'
            STOP
        elseif (lmax /= 2*nmax*(nmax + 2)) then
            write (*, '(A,A)') myname, '> ERROR: lmax /= 2*nmax*(nmax+2)'
            STOP
        elseif (size(scatXYZR, 1) /= 4) then
            write (*, '(A,A)') myname, '> ERROR: size(scatXYZR,1) /= 3'
            STOP
        end if
        !
        nmax2 = nmax
        if (present(nmax2_)) nmax2 = nmax2_
        if (nmax2 < nmax) then
            write (*, '(A,A)') myname, '> ERROR: nmax2 < nmax'
            write (*, *) nmax2, nmax
            STOP
        end if
        pmax2 = nmax2*(nmax2 + 2)
        lmax2 = 2*pmax2
        !
        if (size(sig, 1) /= 3) then
            write (*, '(A,A)') myname, '> ERROR: size(sig,1) /= 3'
            STOP
        elseif (size(sig, 2) == 1) then
            multisplit = .false.
        elseif (size(sig, 2) == nmax2) then
            multisplit = .true.
        else
            write (*, '(A,A)') myname, '> ERROR: Bad size(sig,2)'
            STOP
        end if
        !
        allocate (vtacs(lmax2, lmax2), a(lmax2), f(lmax2))
        allocate (my_fJ(nscat*lmax2))
        !
        tol = 0
        if (present(tol_)) then
            if (tol_ >= 1) then
                write (*, '(A,A)') myname, '> ERROR: tol >= 1'
                STOP
            elseif (tol_ >= 0) then
                tol = tol_
            else
                write (*, '(A,A)') myname, '> ERROR: tol < 0'
                STOP
            end if
        end if
        !
        verb = 0
        if (present(verb_)) verb = verb_
        !
        hostK = sig(1, 1)
        hostK2 = hostK**2 ! host wavenumber squared
        !
        call parseInc( &
            inc=inc, &
            inc_dirn=ipwDirn, &
            inc_ampl=ipwAmpl, &
            verb_=verb)
        call calcCoeffsPW( &
            ipwE0=ipwAmpl, &
            ipwDirn=ipwDirn, &
            ipwCoeffs=a)
        !
        ! Collapse the scattered-field coefficients onto common origin
        f = 0; my_fJ = 0
        do i = 1, nscat
            is = (i - 1)*lmax
            is2 = (i - 1)*lmax2
            my_fJ(is2 + 1:is2 + lmax) = fJ(is + 1:is + lmax)
            if (hostK2*sum(scatXYZR(1:3, i)**2) > tiny) then
                if (verb > 1) then
                    call cpu_time(t0)
                    call system_clock(tic)
                end if
                call calcVTACs( &
                    r0=-scatXYZR(1:3, i), &
                    k=cmplx(hostK, 0, kind(hostK)), &
                    regt=.true., &
                    vtacs=vtacs)
                if (verb > 1) then
                    call cpu_time(t1)
                    call system_clock(toc, count_rate=tps)
                    t2 = t2 + t1 - t0
                    t3 = t3 + dble(toc - tic)/dble(tps)
                end if
            else
                vtacs = 0
                do n = 1, lmax
                    vtacs(n, n) = 1
                end do
            end if
            f = f + matmul(vtacs(1:lmax2, 1:lmax2), my_fJ(is2 + 1:is2 + lmax2))
        end do
        !
        sig = 0
        !
        ! Accumulate the extinction cross-sections
        csExt = 0; nconvExt = 0
        ext_loop: do n = 1, nmax2
            if (nconvExt > 0) exit ext_loop
            is = n*(n + 1)
            csExt_n = realpart(dot_product(a(is - n:is + n), f(is - n:is + n)))
            is = is + pmax2
            csExt_n = csExt_n + realpart(dot_product(a(is - n:is + n), f(is - n:is + n)))
            csExt = csExt + csExt_n
            if (abs(csExt) > 0 .and. abs(csExt_n) > 0) then
                if (abs(csExt_n/csExt) < tol) nconvExt = n
            end if
            if (multisplit) then
                sig(1, n) = csExt_n
            else
                sig(1, 1) = sig(1, 1) + csExt_n
            end if
        end do ext_loop
        sig(1, :) = -sig(1, :)/hostK2
        ! why is it twice what i expect?!?!?!?!?
        !
        ! Accumulate the scattering cross-sections
        csSca = 0; nconvSca = 0
        sca_loop: do n = 1, nmax2
            if (nconvSca > 0) exit sca_loop
            is = n*(n + 1)
            csSca_n = realpart(dot_product(f(is - n:is + n), f(is - n:is + n)))
            is = is + pmax2
            csSca_n = csSca_n + realpart(dot_product(f(is - n:is + n), f(is - n:is + n)))
            csSca = csSca + csSca_n
            if (abs(csSca) > 0 .and. abs(csSca_n) > 0) then
                if (abs(csSca_n/csSca) < tol) nconvSca = n
            end if
            if (multisplit) then
                sig(2, n) = csSca_n
            else
                sig(2, 1) = sig(2, 1) + csSca_n
            end if
        end do sca_loop
        sig(2, :) = sig(2, :)/hostK2
        !
        ! Use energy conservation for absorption cross-section
        sig(3, :) = sig(1, :) - sig(2, :)
        !
        if (verb > 0) then
            write (*, '(A,A,es15.8e2,A,i2)') &
                myname, '> csExt= ', sum(sig(1, :)), ' nConv= ', nconvExt
            write (*, '(9x,A,es15.8e2,A,i2)') &
                'csSca= ', sum(sig(2, :)), ' nConv= ', nconvSca
            write (*, '(9x,A,es15.8e2)') 'csAbs= ', sum(sig(3, :))
            write (*, *)
            if (verb > 1) then
                write (*, '(A,A,2(1x,es10.3E2))') &
                    myName, '> Calculation time [calcVTACs] (CPU & real in s): ', &
                    t2, t3
            end if

        end if
        !
    end subroutine calcCs
    !
    !===============================================================
    !
    subroutine calcOAprops(Tmat, sigOA, verb_, rtol_)
        !
        use swav, only: calcWignerBigD, calcVTACsAxial, calcVTACs
        !
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables:
        complex(8), intent(in) :: Tmat(:, :)
        real(8), intent(inout) :: sigOA(:, 0:)
        integer, intent(out), optional :: verb_
        real(8), intent(in), optional :: rtol_
        ! Local variables
        character(*), parameter :: myName = 'calcOAprops'
        integer :: i, j, l, p, q, lmax, nmax, pmax, verb
        integer :: nconv(2)
        real(8) :: pref, hostK, rdum(4), rtol ! debug1, debug2, debug3, debug4
        complex(8), dimension(size(Tmat, 1), size(Tmat, 2)) :: Mat ! lmax by lmax
        logical :: split
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        verb = 0
        if (present(verb_)) verb = verb_
        !
        if (size(Tmat, 1) /= size(Tmat, 2)) then
            write (*, '(A,A)') myName, '> ERROR: Supplied Tmat not square!'
            STOP
        end if
        lmax = size(Tmat, 1)
        pmax = lmax/2
        nmax = int(sqrt(dble(pmax)))
        if (2*nmax*(nmax + 2) /= size(Tmat, 1)) then
            write (*, '(A,A)') myName, '> ERROR: size(Tmat,1) /= 2*nmax*(nmax+2)'
            STOP
        end if
        !
        rtol = 0
        if (present(rtol_)) then
            if (rtol_ >= 1) then
                write (*, '(A,A)') 'ERROR: superunitary rtol= ', rtol_
                STOP
            elseif (rtol_ >= 0) then
                rtol = rtol_
            else
                write (*, '(A,A)') 'ERROR: negative rtol= ', rtol_
                STOP
            end if
        end if
        !
        if (size(sigOA, 1) /= 6) then
            write (*, '(A,A)') myName, '> ERROR: size(sigOA,1) /= 6'
            if (present(verb_)) verb_ = -3
            !return
            STOP
        elseif (size(sigOA, 2) /= 1 .and. size(sigOA, 2) /= nmax + 1) then
            write (*, '(A,A)') myName, '> ERROR: size(sigOA,2) is neither 1 nor nmax+1'
            write (*, '(A,A,i2)') myname, '> nmax= ', nmax
            write (*, '(A,A,i2)') myname, '> size(sigOA,2)= ', size(sigOA, 2)
            if (present(verb_)) verb_ = -3
            !return
            STOP
        end if
        !
        if (size(sigOA, 2) == nmax + 1) then
            split = .true.
        else
            split = .false.
        end if
        !
        hostK = sigOA(1, 0)
        pref = tpi/hostK**2
        !

      !   call dumpTmat(tmat=Tmat, &
      !   filename= 'collective_parity.txt', &
      !    eps_med=1.0d0, &
      !    lambda=1.0d0, &
      !     verb_=verb)


        Mat = matmul(conjg(transpose(Tmat)), Tmat)

      !   call dumpTmat(tmat=Mat, &
      !   filename= 'collective_TT.txt', &
      !    eps_med=1.0d0, &
      !    lambda=1.0d0, &
      !     verb_=verb)

        sigOA = 0
        nconv = 0
        do i = 1, nmax
            rdum = 0
            do q = 0, 1
                do j = -i, i
                    p = i*(i + 1) + j
                    l = pmax*q + p
                    !
                    ! Extinction
                    if (nconv(1) == 0) rdum(1) = rdum(1) + realpart(Tmat(l, l))
                    ! Scattering
                    if (nconv(2) == 0) rdum(2) = rdum(2) + realpart(Mat(l, l))
                    !
                end do
            end do
            !
            quant1: do j = 1, 2
                if (nconv(j) == 0) then
                    sigOA(j, 0) = sigOA(j, 0) + rdum(j)
                    if (split) then
                        sigOA(j, i) = rdum(j)
                    end if
                    if (abs(sigOA(j, 0)) > 0 .and. abs(rdum(j)) > 0) then
                        if (rdum(j)/sigOA(j, 0) < rtol) nconv(j) = i
                    end if
                end if
            end do quant1
            !
        end do
        sigOA(1, :) = -pref*sigOA(1, :)
        sigOA(2, :) = pref*sigOA(2, :)
        sigOA(3, :) = sigOA(1, :) - sigOA(2, :)
        !
        if (verb > 0) then
            write (*, '(A,A,es15.8e2,A,i2)') &
                myName, '> <csExt>= ', sigOA(1, 0), " nConv= ", nconv(1)
            write (*, '(13x,A,es15.8e2,A,i2)') &
                '<csSca>= ', sigOA(2, 0), " nConv= ", nconv(2)
            write (*, '(13x,A,es15.8e2)') '<csAbs>= ', sigOA(3, 0)
            write (*, *)
        end if
        !
        ! Transform Tmat from "parity" (M-N) basis to "helicity" (L-R) basis,
        ! following SuryadharmaR18, and store in Mat
        !
        Mat(1:pmax, 1:pmax) = 0.5d0*( &
                              Tmat(1:pmax, 1:pmax) + &
                              Tmat(1:pmax, pmax + 1:2*pmax) + &
                              Tmat(pmax + 1:2*pmax, 1:pmax) + &
                              Tmat(pmax + 1:2*pmax, pmax + 1:2*pmax)) ! This is the LL block
        !
        Mat(1:pmax, pmax + 1:2*pmax) = 0.5d0*( &
                                       Tmat(1:pmax, 1:pmax) - &
                                       Tmat(1:pmax, pmax + 1:2*pmax) + &
                                       Tmat(pmax + 1:2*pmax, 1:pmax) - &
                                       Tmat(pmax + 1:2*pmax, pmax + 1:2*pmax)) ! This is the LR block
        !
        Mat(pmax + 1:2*pmax, 1:pmax) = 0.5d0*( &
                                       Tmat(1:pmax, 1:pmax) + &
                                       Tmat(1:pmax, pmax + 1:2*pmax) - &
                                       Tmat(pmax + 1:2*pmax, 1:pmax) - &
                                       Tmat(pmax + 1:2*pmax, pmax + 1:2*pmax)) ! This is the RL block
        !
        Mat(pmax + 1:2*pmax, pmax + 1:2*pmax) = 0.5d0*( &
                                                Tmat(1:pmax, 1:pmax) - &
                                                Tmat(1:pmax, pmax + 1:2*pmax) - &
                                                Tmat(pmax + 1:2*pmax, 1:pmax) + &
                                                Tmat(pmax + 1:2*pmax, pmax + 1:2*pmax)) ! This is the LL block
        !
        ! Use formulae in Table 2 of SuryadharmaR18 with S=1, 4 -> 2, and minus
        ! in front of extinction.
        !
        ! Because it is more challenging than average cross-sections, circular dichroism will be summed
        ! over all n <= nmax.
      !   if (nconv(1) == 0) nconv(1) = nmax
      !   if (nconv(2) == 0) nconv(2) = nmax
        
         ! call dumpTmat(tmat=Mat, &
         !               filename= 'collective_helicity.txt', &
         !                eps_med=1.0d0, &
         !                lambda=1.0d0, &
         !                 verb_=verb)

        do i = 1, nmax ! order
            !
            !write(*,*) 'nconv=', nconv
            !
            rdum = 0
            !
            ! Extinction CD
            if (i <= nmax) then
                do j = -i, i  ! degree
                    p = i*(i + 1) + j
                    rdum(1) = rdum(1) + & ! CD
                              realpart(Mat(p, p)) - realpart(Mat(pmax + p, pmax + p))
                end do
                sigOA(4, 0) = sigOA(4, 0) + rdum(1) ! CD
                if (split) then
                    sigOA(4, i) = rdum(1)
                end if
            end if
            !
            ! Scattering CD
            if (i <= nmax) then
                do j = -i, i  ! degree
                    p = i*(i + 1) + j

                    rdum(2) = rdum(2) + & ! CD
                              dot_product(Mat(1:pmax, p), Mat(1:pmax, p)) - & !trace(TLL)
                              dot_product(Mat(1:pmax, pmax + p), Mat(1:pmax, pmax + p)) - & !trace(TLR)
                              dot_product(Mat(pmax + 1:2*pmax, pmax + p), & !trace(TRR)
                                          Mat(pmax + 1:2*pmax, pmax + p)) + &
                              dot_product(Mat(pmax + 1:2*pmax, p), Mat(pmax + 1:2*pmax, p)) !trace(TRL)
                end do
                sigOA(5, 0) = sigOA(5, 0) + rdum(2) ! CD
                if (split) then
                    sigOA(5, i) = rdum(2)
                end if
            end if
            !
        end do

        sigOA(4, :) = -pref*sigOA(4, :)
        sigOA(5, :) = pref*sigOA(5, :)
        sigOA(6, :) = sigOA(4, :) - sigOA(5, :)
        !
        sigOA(4:6, :) = -2*sigOA(4:6, :) ! Extra factor to change CD "convention"
        !
        if (verb > 0) then
            write (*, '(13x,A,1x,es15.8e2)') '<cdExt>=', sigOA(4, 0)
            write (*, '(13x,A,1x,es15.8e2)') '<cdSca>=', sigOA(5, 0)
            write (*, '(13x,A,1x,es15.8e2)') '<cdAbs>=', sigOA(6, 0)
            write (*, *)
            !
        end if
        !
    end subroutine calcOAprops
    !

    subroutine alphaTensor(T, Alpha)

        ! Conversion of l<=3 spherical multipoles into cartesian multipoles
        ! Describing Meta-Atoms Using the Exact Higher-Order Polarizability Tensors
        ! ACS Photonics 2020, 7, 1153−1162
        ! https://dx.doi.org/10.1021/acsphotonics.9b01776
        ! alpha = c_n * c_np / (im*k^3) * V * T * U^-1
        ! ignoring im*k^3 prefactor here

        ! use linalg, only: invSqrMat

        complex(8), intent(in) :: T(:, :)
        complex(8), dimension(30, 30) :: Alpha(:, :) ! (3 + 5 + 7) x 2 -> 30x30 matrix

        integer :: lmax, pmax

        complex(8), dimension(3, 3) :: V1, iU1
        complex(8), dimension(5, 5) :: V2, iU2
        complex(8), dimension(7, 7) :: V3, iU3
        real(8), parameter :: sqrt2 = sqrt(2.0d0), sqrt2o3 = sqrt(2.0d0/3.0d0), &
                              sqrt5 = sqrt(5.0d0), sqrt15 = sqrt(15.0d0), &
                              isqrt5 = 1.0/sqrt(5.0d0), isqrt15 = 1.0/sqrt(15.0d0)
        ! complex(8), parameter :: im = (0.0, 1.0)
        real(8), parameter :: c_1 = sqrt(6.0*pi), c_2 = sqrt(20.0*pi), c_3 = sqrt(105.0*pi/2.0)

        ! V1 = reshape((/ 1.0, -imu, 0.0, &
        !                 0.0, 0.0, sqrt2, &
        !                 -1.0, -imu, 0.0 /), shape(V1))

        ! V2 = reshape((/ -imu, 0.0, 0.0, 1.0, -1.0, &
        !                 0.0, 1.0, -imu, 0.0, 0.0, &
        !                 0.0, 0.0, 0.0, -sqrt2o3, -sqrt2o3, &
        !                 0.0, 1.0, -imu, 0.0, 0.0, &
        !                imu, 0.0, 0.0, 1.0, -1.0 /), shape(V2))

        ! V3 = reshape((/ -1.0,0.0, 0.0, -imu, 0.0, 0.0, 0.0, &
        !                 0.0, 0.0, 0.0, 0.0, sqrt2o3, -sqrt2o3, -sqrt2o3*imu, &
        !                 -isqrt15, 4*isqrt15, -4*imu/isqrt15, imu/isqrt15, 0.0, 0.0, 0.0, &
        !                 0.0, 0.0, 0.0, 0.0, -2*isqrt5, -2*isqrt5, 0.0, &
        !                 isqrt15, -4*isqrt15, -4*imu*isqrt15, imu*isqrt15, 0.0, 0.0, 0.0, &
        !                 0.0, 0.0, 0.0, 0.0, sqrt2o3, -sqrt2o3, imu*sqrt2o3, &
        !                 1.0, 0.0, 0.0, -imu, 0.0, 0.0, 0.0/), shape(V3))

        V1 = reshape((/1.0 + 0*imu, -imu, 0*imu, &
                       0*imu, 0*imu, sqrt2 + 0*imu, &
                       -1.0 + 0*imu, -imu, 0*imu/), shape(V1))

        V2 = reshape((/-imu, 0*imu, 0*imu, 1.0 + 0*imu, -1.0 + 0*imu, &
                       0*imu, 1.0 + 0*imu, -imu, 0*imu, 0*imu, &
                       0*imu, 0*imu, 0*imu, -sqrt2o3 + 0*imu, -sqrt2o3 + 0*imu, &
                       0*imu, 1.0 + 0*imu, -imu, 0*imu, 0*imu, &
                       imu, 0*imu, 0*imu, 1.0 + 0*imu, -1.0 + 0*imu/), shape(V2))

        V3 = reshape((/-1.0 + 0*imu, 0*imu, 0*imu, -imu, 0*imu, 0*imu, 0*imu, &
                       0*imu, 0*imu, 0*imu, 0*imu, sqrt2o3 + 0*imu, -sqrt2o3 + 0*imu, -sqrt2o3 + 0*imu*imu, &
                       -isqrt15 + 0*imu, 4*isqrt15 + 0*imu, -4*imu/isqrt15, imu/isqrt15, 0*imu, 0*imu, 0*imu, &
                       0*imu, 0*imu, 0*imu, 0*imu, -2*isqrt5 + 0*imu, -2*isqrt5 + 0*imu, 0*imu, &
                       isqrt15 + 0*imu, -4*isqrt15 + 0*imu, -4*imu*isqrt15, imu*isqrt15, 0*imu, 0*imu, 0*imu, &
                       0*imu, 0*imu, 0*imu, 0*imu, sqrt2o3 + 0*imu, -sqrt2o3 + 0*imu, imu*sqrt2o3, &
                       1.0 + 0*imu, 0*imu, 0*imu, -imu, 0*imu, 0*imu, 0*imu/), shape(V3))

        V1 = 1.0/sqrt2*V1
        V2 = 1.0/2.0*V2
        V3 = 1.0/sqrt(8.0)*V3

! U = (V^-1)^h
! and we need U^-1
! it seems this actually simplifies to U^-1 == V^h

! U1 = V1
! U2 = V2
! U3 = V3
        ! call invSqrMat( &
        !    A = U1, &
        !    trans_ = .true., &
        !    verb_ = .true.)
        ! U1 = conj(U1)

        iU1 = conjg(transpose(V1))
        iU2 = conjg(transpose(V2))
        iU3 = conjg(transpose(V3))

        lmax = size(T, 1)
        pmax = lmax/2
               
       ! ! ! E-E block
       ! Alpha(1:3,1:3)   = c_1 * c_1 * matmul(V1 , matmul(T(1:3,1:3)  , iU1)) ! e1-e1
       ! Alpha(1:3,4:8)   = c_1 * c_2 * matmul(V1 , matmul(T(1:3,4:8)  , iU2)) ! e1-e2
       ! Alpha(1:3,9:15)  = c_1 * c_3 * matmul(V1 , matmul(T(1:3,9:15) , iU3)) ! e1-e3
       ! Alpha(4:8,1:3)   = c_2 * c_1 * matmul(V2 , matmul(T(4:8,1:3)  , iU1)) ! e2-e1
       ! Alpha(4:8,4:8)   = c_2 * c_2 * matmul(V2 , matmul(T(4:8,4:8)  , iU2)) ! e2-e2
       ! Alpha(4:8,9:15)  = c_2 * c_3 * matmul(V2 , matmul(T(4:8,9:15) , iU3)) ! e2-e3
       ! Alpha(9:15,1:3)  = c_3 * c_1 * matmul(V3 , matmul(T(9:15,1:3) , iU1)) ! e3-e1
       ! Alpha(9:15,4:8)  = c_3 * c_2 * matmul(V3 , matmul(T(9:15,4:8) , iU2)) ! e3-e2
       ! Alpha(9:15,9:15) = c_3 * c_3 * matmul(V3 , matmul(T(9:15,9:15), iU3)) ! e3-e3
       
       ! ! ! ! M-M block
       ! Alpha(15 + 1:15 + 3,  15 + 1:15 + 3)   = c_1 * c_1 * matmul(V1 , matmul(T(pmax + 1:pmax + 3, pmax + 1:pmax + 3)  , iU1  ))! m1-m1
       ! Alpha(15 + 1:15 + 3,  15 + 4:15 + 8)   = c_1 * c_2 * matmul(V1 , matmul(T(pmax + 1:pmax + 3, pmax + 4:pmax + 8)  , iU2  ))! m1-m2
       ! Alpha(15 + 1:15 + 3,  15 + 9:15 + 15)  = c_1 * c_3 * matmul(V1 , matmul(T(pmax + 1:pmax + 3, pmax + 9:pmax + 15) , iU3 ))! m1-m3
       ! Alpha(15 + 4:15 + 8,  15 + 1:15 + 3)   = c_2 * c_1 * matmul(V2 , matmul(T(pmax + 4:pmax + 8, pmax + 1:pmax + 3)  , iU1  ))! m2-m1
       ! Alpha(15 + 4:15 + 8,  15 + 4:15 + 8)   = c_2 * c_2 * matmul(V2 , matmul(T(pmax + 4:pmax + 8, pmax + 4:pmax + 8)  , iU2  ))! m2-m2
       ! Alpha(15 + 4:15 + 8,  15 + 9:15 + 15)  = c_2 * c_3 * matmul(V2 , matmul(T(pmax + 4:pmax + 8, pmax + 9:pmax + 15) , iU3 ))! m2-m3
       ! Alpha(15 + 9:15 + 15, 15 + 1:15 + 3)   = c_3 * c_1 * matmul(V3 , matmul(T(pmax + 9:pmax + 15, pmax + 1:pmax + 3) , iU1  ))! m3-m1
       ! Alpha(15 + 9:15 + 15, 15 + 4:15 + 8)   = c_3 * c_2 * matmul(V3 , matmul(T(pmax + 9:pmax + 15, pmax + 4:pmax + 8) , iU2  ))! m3-m2
       ! Alpha(15 + 9:15 + 15, 15 + 9:15 + 15)  = c_3 * c_3 * matmul(V3 , matmul(T(pmax + 9:pmax + 15, pmax + 9:pmax + 15), iU3 ))! m3-m3
       
       ! ! ! E-M block
       ! Alpha(1:3,  15 + 1:15 + 3)  = c_1 * c_1 * matmul(V1 , matmul(T(1:3,  pmax + 1:pmax + 3) , iU1  ))! e1-m1
       ! Alpha(1:3,  15 + 4:15 + 8)  = c_1 * c_2 * matmul(V1 , matmul(T(1:3,  pmax + 4:pmax + 8) , iU2  ))! e1-m2
       ! Alpha(1:3,  15 + 9:15 + 15) = c_1 * c_3 * matmul(V1 , matmul(T(1:3,  pmax + 9:pmax + 15),  iU3 ))! e1-m3
       ! Alpha(4:8,  15 + 1:15 + 3)  = c_2 * c_1 * matmul(V2 , matmul(T(4:8,  pmax + 1:pmax + 3) , iU1  ))! e2-m1
       ! Alpha(4:8,  15 + 4:15 + 8)  = c_2 * c_2 * matmul(V2 , matmul(T(4:8,  pmax + 4:pmax + 8) , iU2  ))! e2-m2
       ! Alpha(4:8,  15 + 9:15 + 15) = c_2 * c_3 * matmul(V2 , matmul(T(4:8,  pmax + 9:pmax + 15),  iU3 ))! e2-m3
       ! Alpha(9:15, 15 + 1:15 + 3)  = c_3 * c_1 * matmul(V3 , matmul(T(9:15, pmax + 1:pmax + 3) , iU1  ))! e3-m1
       ! Alpha(9:15, 15 + 4:15 + 8)  = c_3 * c_2 * matmul(V3 , matmul(T(9:15, pmax + 4:pmax + 8) , iU2  ))! e3-m2
       ! Alpha(9:15, 15 + 9:15 + 15) = c_3 * c_3 * matmul(V3 , matmul(T(9:15, pmax + 9:pmax + 15),  iU3 ))! e3-m3
       
       ! ! M-E block
       ! Alpha(15 + 1:15 + 3,  1:3)  = c_1 * c_1 * matmul(V1 , matmul(T(pmax + 1:pmax + 3,  1:3) , iU1  ))! m1-e1
       ! Alpha(15 + 1:15 + 3,  4:8)  = c_1 * c_2 * matmul(V1 , matmul(T(pmax + 1:pmax + 3,  4:8) , iU2  ))! m1-e2
       ! Alpha(15 + 1:15 + 3,  9:15) = c_1 * c_3 * matmul(V1 , matmul(T(pmax + 1:pmax + 3,  9:15),  iU3 ))! m1-e3
       ! Alpha(15 + 4:15 + 8,  1:3)  = c_2 * c_1 * matmul(V2 , matmul(T(pmax + 4:pmax + 8,  1:3) , iU1  ))! m2-e1
       ! Alpha(15 + 4:15 + 8,  4:8)  = c_2 * c_2 * matmul(V2 , matmul(T(pmax + 4:pmax + 8,  4:8) , iU2  ))! m2-e2
       ! Alpha(15 + 4:15 + 8,  9:15) = c_2 * c_3 * matmul(V2 , matmul(T(pmax + 4:pmax + 8,  9:15),  iU3 ))! m2-e3
       ! Alpha(15 + 9:15 + 15, 1:3)  = c_3 * c_1 * matmul(V3 , matmul(T(pmax + 9:pmax + 15, 1:3) , iU1  ))! m3-e1
       ! Alpha(15 + 9:15 + 15, 4:8)  = c_3 * c_2 * matmul(V3 , matmul(T(pmax + 9:pmax + 15, 4:8) , iU2  ))! m3-e2
       ! Alpha(15 + 9:15 + 15, 9:15) = c_3 * c_3 * matmul(V3 , matmul(T(pmax + 9:pmax + 15, 9:15),  iU3 ))! m3-e3

        Alpha = 0*imu ! initialise to zeros, in case Nmax(T) < 3

        ! l=1, always
        Alpha(1:3, 1:3) = c_1*c_1*matmul(V1, matmul(T(1:3, 1:3), iU1)) ! e1-e1
        Alpha(15 + 1:15 + 3, 15 + 1:15 + 3) = c_1*c_1*matmul(V1, matmul(T(pmax + 1:pmax + 3, pmax + 1:pmax + 3), iU1))! m1-m1
        Alpha(1:3, 15 + 1:15 + 3) = c_1*c_1*matmul(V1, matmul(T(1:3, pmax + 1:pmax + 3), iU1))! e1-m1
        Alpha(15 + 1:15 + 3, 1:3) = c_1*c_1*matmul(V1, matmul(T(pmax + 1:pmax + 3, 1:3), iU1))! m1-e1

        if (pmax > 3) then ! l=2

            Alpha(1:3, 4:8) = c_1*c_2*matmul(V1, matmul(T(1:3, 4:8), iU2)) ! e1-e2
            Alpha(4:8, 1:3) = c_2*c_1*matmul(V2, matmul(T(4:8, 1:3), iU1)) ! e2-e1
            Alpha(4:8, 4:8) = c_2*c_2*matmul(V2, matmul(T(4:8, 4:8), iU2)) ! e2-e2
            Alpha(15 + 1:15 + 3, 15 + 4:15 + 8) = c_1*c_2*matmul(V1, matmul(T(pmax + 1:pmax + 3, pmax + 4:pmax + 8), iU2))! m1-m2
            Alpha(15 + 4:15 + 8, 15 + 1:15 + 3) = c_2*c_1*matmul(V2, matmul(T(pmax + 4:pmax + 8, pmax + 1:pmax + 3), iU1))! m2-m1
            Alpha(15 + 4:15 + 8, 15 + 4:15 + 8) = c_2*c_2*matmul(V2, matmul(T(pmax + 4:pmax + 8, pmax + 4:pmax + 8), iU2))! m2-m2
            Alpha(1:3, 15 + 4:15 + 8) = c_1*c_2*matmul(V1, matmul(T(1:3, pmax + 4:pmax + 8), iU2))! e1-m2
            Alpha(4:8, 15 + 1:15 + 3) = c_2*c_1*matmul(V2, matmul(T(4:8, pmax + 1:pmax + 3), iU1))! e2-m1
            Alpha(15 + 1:15 + 3, 4:8) = c_1*c_2*matmul(V1, matmul(T(pmax + 1:pmax + 3, 4:8), iU2))! m1-e2
            Alpha(15 + 4:15 + 8, 1:3) = c_2*c_1*matmul(V2, matmul(T(pmax + 4:pmax + 8, 1:3), iU1))! m2-e1
            Alpha(15 + 4:15 + 8, 4:8) = c_2*c_2*matmul(V2, matmul(T(pmax + 4:pmax + 8, 4:8), iU2))! m2-e2

        end if

        if (pmax > 8) then ! l=3

            Alpha(1:3, 9:15) = c_1*c_3*matmul(V1, matmul(T(1:3, 9:15), iU3)) ! e1-e3
            Alpha(4:8, 9:15) = c_2*c_3*matmul(V2, matmul(T(4:8, 9:15), iU3)) ! e2-e3
            Alpha(9:15, 1:3) = c_3*c_1*matmul(V3, matmul(T(9:15, 1:3), iU1)) ! e3-e1
            Alpha(9:15, 4:8) = c_3*c_2*matmul(V3, matmul(T(9:15, 4:8), iU2)) ! e3-e2
            Alpha(9:15, 9:15) = c_3*c_3*matmul(V3, matmul(T(9:15, 9:15), iU3)) ! e3-e3
            Alpha(15 + 1:15 + 3, 15 + 9:15 + 15) = c_1*c_3*matmul(V1, matmul(T(pmax + 1:pmax + 3, pmax + 9:pmax + 15), iU3))! m1-m3
            Alpha(15 + 4:15 + 8, 15 + 9:15 + 15) = c_2*c_3*matmul(V2, matmul(T(pmax + 4:pmax + 8, pmax + 9:pmax + 15), iU3))! m2-m3
            Alpha(15 + 9:15 + 15, 15 + 1:15 + 3) = c_3*c_1*matmul(V3, matmul(T(pmax + 9:pmax + 15, pmax + 1:pmax + 3), iU1))! m3-m1
            Alpha(15 + 9:15 + 15, 15 + 4:15 + 8) = c_3*c_2*matmul(V3, matmul(T(pmax + 9:pmax + 15, pmax + 4:pmax + 8), iU2))! m3-m2
            Alpha(15 + 9:15 + 15, 15 + 9:15 + 15) = c_3*c_3*matmul(V3, matmul(T(pmax + 9:pmax + 15, pmax + 9:pmax + 15), iU3))! m3-m3
            Alpha(1:3, 15 + 9:15 + 15) = c_1*c_3*matmul(V1, matmul(T(1:3, pmax + 9:pmax + 15), iU3))! e1-m3
            Alpha(4:8, 15 + 4:15 + 8) = c_2*c_2*matmul(V2, matmul(T(4:8, pmax + 4:pmax + 8), iU2))! e2-m2
            Alpha(4:8, 15 + 9:15 + 15) = c_2*c_3*matmul(V2, matmul(T(4:8, pmax + 9:pmax + 15), iU3))! e2-m3
            Alpha(9:15, 15 + 1:15 + 3) = c_3*c_1*matmul(V3, matmul(T(9:15, pmax + 1:pmax + 3), iU1))! e3-m1
            Alpha(9:15, 15 + 4:15 + 8) = c_3*c_2*matmul(V3, matmul(T(9:15, pmax + 4:pmax + 8), iU2))! e3-m2
            Alpha(9:15, 15 + 9:15 + 15) = c_3*c_3*matmul(V3, matmul(T(9:15, pmax + 9:pmax + 15), iU3))! e3-m3
            Alpha(15 + 1:15 + 3, 9:15) = c_1*c_3*matmul(V1, matmul(T(pmax + 1:pmax + 3, 9:15), iU3))! m1-e3
            Alpha(15 + 4:15 + 8, 9:15) = c_2*c_3*matmul(V2, matmul(T(pmax + 4:pmax + 8, 9:15), iU3))! m2-e3
            Alpha(15 + 9:15 + 15, 1:3) = c_3*c_1*matmul(V3, matmul(T(pmax + 9:pmax + 15, 1:3), iU1))! m3-e1
            Alpha(15 + 9:15 + 15, 4:8) = c_3*c_2*matmul(V3, matmul(T(pmax + 9:pmax + 15, 4:8), iU2))! m3-e2
            Alpha(15 + 9:15 + 15, 9:15) = c_3*c_3*matmul(V3, matmul(T(pmax + 9:pmax + 15, 9:15), iU3))! m3-e3

        end if

    end subroutine alphaTensor


    subroutine contractTmat(Tin, scatXYZR, rtr, Tout, verb_, mack_)
        !
        use swav, only: xyz2rtp, calcWignerBigD, calcVTACsAxial, calcVTACs
        !
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables:
        real(8), intent(in) :: scatXYZR(:, :)
        complex(8), intent(in) :: Tin(:, :)
        complex(8), intent(inout) :: Tout(:, :)
        complex(8), dimension(30, 30) :: Alpha
        logical, intent(in) :: rtr
        integer, intent(inout), optional :: verb_
        logical, intent(in), optional :: mack_
        ! Local variables
        character(*), parameter :: myName = 'contractTmat'
        integer :: i, j, js, is, pmax2
        integer :: tic, toc, tps, nscat, lmax, nmax, pmax, verb
        real(8) :: rtp(3), rdum, hostK, lambda, hostEps
        real    :: t0, t1, t2, t4, t6, t8
        real(8) :: t3, t5, t7, t9
        complex(8), dimension(size(Tout, 1), size(Tout, 2)) :: Mat, MatJ
        logical :: shift(size(scatXYZR, 2))
        complex(8) :: vtacs(size(Tout, 1), size(Tout, 2), size(scatXYZR, 2)), off(3, 1)
        complex(8), allocatable :: bigD(:, :, :)
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        if (size(Tin, 1) /= size(Tin, 2) .and. &
            size(Tin, 1) /= size(Tin, 2)*size(scatXYZR, 2)) then
            write (*, '(A,A)') myName, '> ERROR: Badly shaped Tin'
            STOP
        elseif (size(scatXYZR, 1) /= 4) then
            write (*, '(A,A)') myName, '> ERROR: size(scatXYZ,1) /= 3'
            STOP
        end if
        nscat = size(scatXYZR, 2)
        lmax = size(Tin, 1)/nscat
        pmax = lmax/2
        nmax = int(sqrt(dble(pmax)))
        if (nscat*2*nmax*(nmax + 2) /= size(Tin, 1)) then
            write (*, '(A,A)') myName, '> ERROR: size(Tin,1) /= nscat*2*nmax*(nmax+2)'
            STOP
        end if
        !
        if (size(Tout, 1) /= size(Tout, 2)) then
            write (*, '(A,A)') myName, '> ERROR: Supplied Tout not square!'
            STOP
        elseif (size(Tout, 1) < lmax) then
            write (*, '(A,A)') myName, '> ERROR: size(Tout,1) < size(Tin,1)/nscat'
            STOP
        end if
        !
        pmax2 = size(Tout, 1)/2 ! used for padding with zeroes when
        j = int(sqrt(dble(pmax2)))
        if (2*j*(j + 2) /= size(Tout, 1)) then
            write (*, '(A,A)') myName, '> ERROR: size(Tout,1) /= 2*n*(n+2)'
            STOP
        end if
        !
        !hostK = realpart(Tout(1,1))
        lambda = realpart(Tout(1, 1))
        hostEps = imagpart(Tout(1, 1))
        hostK = tpi/lambda*sqrt(hostEps)
        !
        verb = 0
        if (present(verb_)) verb = verb_
        !
        shift(1) = .false.
        if (present(mack_)) shift(1) = mack_
        !
        t2 = 0.0d0
        t3 = 0.0d0
        t4 = 0.0d0
        t5 = 0.0d0
        t6 = 0.0d0
        t7 = 0.0d0
        t8 = 0.0d0
        t9 = 0.0d0
        if (size(Tin, 1) /= size(Tin, 2) .or. shift(1)) then ! Mackowski's TJ matrices
            rdum = 0
            Tout = 0
            do j = 1, nscat
                js = (j - 1)*lmax
                off(1:3, 1) = -cmplx(hostK, 0, 8)*scatXYZR(1:3, j)  ! regular
                Mat = 0
                !Mat(1:lmax,1:lmax) = Tin(js+1:js+lmax,1:lmax) ! WRONG!!
                Mat(1:pmax, 1:pmax) = Tin(js + 1:js + pmax, 1:pmax)
                Mat(pmax2 + 1:pmax2 + pmax, 1:pmax) = &
                    Tin(js + pmax + 1:js + lmax, 1:pmax)
                Mat(1:pmax, pmax2 + 1:pmax2 + pmax) = &
                    Tin(js + 1:js + pmax, pmax + 1:lmax)
                Mat(pmax2 + 1:pmax2 + pmax, pmax2 + 1:pmax2 + pmax) = &
                    Tin(js + pmax + 1:js + lmax, pmax + 1:lmax)
                call xyz2rtp(scatXYZR(1:3, j), rtp)
                if (rtp(1)*hostK > tiny) then

                    if (verb > 1) then
                        call cpu_time(t0)
                        call system_clock(tic)
                    end if
                    !
                    call offsetTmat( &
                        Tmat=Mat, &
                        off=off, &
                        miet=.false., &
                        rtr=rtr, &
                        right=.false.)

                    if (verb > 1) then
                        call cpu_time(t1)
                        call system_clock(toc, count_rate=tps)
                        t2 = t2 + t1 - t0
                        t3 = t3 + dble(toc - tic)/dble(tps)
                    end if
                    !
                end if
                ! Accumulate the collective
                Tout = Tout + Mat
                !
                ! Start of sanity check >>>
                ! if(rtp(1)*hostK > tiny) then
                !    call calcVTACs( &
                !         r0 = -scatXYZR(1:3,j), &
                !         k = cmplx(hostK,0,8), &
                !         regt = .true., &
                !         vtacs = matJ )
                !    rdum = rdum + &
                !         realpart(sum(transpose(MatJ)*Tin(js+1:js+lmax,1:lmax)))
                ! else
                !    do i=1,lmax
                !       rdum = rdum + realpart(Tin(js+i,i))
                !    enddo
                ! endif
                ! <<< end of sanity check
                !
            end do
            !
        else ! Stout's pairwise TIJ matrices
            !
            allocate (bigD(pmax2, pmax2, nscat))
            ! >>> Could still make the rotaton optional, so that it can be
            ! skipped for z-linear systems, but the speed-up won't be great.
            !
            do j = 1, nscat
                call xyz2rtp(scatXYZR(1:3, j), rtp)
                if (rtp(1)*hostK > tiny) then
                    shift(j) = .true.
                    if (rtr) then
                        if (verb > 1) then
                            call cpu_time(t0)
                            call system_clock(tic)
                        end if
                        call calcWignerBigD( &
                            angles=(/rtp(3), rtp(2), 0.0d0/), &
                            pmax=pmax2, &
                            bigD=bigD(:, :, j))
                        if (verb > 1) then
                            call cpu_time(t1)
                            call system_clock(toc, count_rate=tps)
                            t4 = t4 + t1 - t0
                            t5 = t5 + dble(toc - tic)/dble(tps)
                            call cpu_time(t0)
                            call system_clock(tic)
                        end if
                        !
                        call calcVTACsAxial( &
                            r0=rtp(1), &
                            k=cmplx(hostK, 0, kind(hostK)), &
                            pmax=pmax2, &
                            regt=.true., & ! regular variant
                            flip=.false., & ! parallel to local z axis
                            vtacs=vtacs(:, :, j))
                        !
                        if (verb > 1) then
                            call cpu_time(t1)
                            call system_clock(toc, count_rate=tps)
                            t6 = t6 + t1 - t0
                            t7 = t7 + dble(toc - tic)/dble(tps)
                        end if
                    else
                        if (verb > 1) then
                            call cpu_time(t0)
                            call system_clock(tic)
                        end if
                        call calcVTACs( &
                            r0=scatXYZR(1:3, j), &
                            k=cmplx(hostK, 0, kind(hostK)), &
                            regt=.true., &
                            vtacs=vtacs(:, :, j))
                        if (verb > 1) then
                            call cpu_time(t1)
                            call system_clock(toc, count_rate=tps)
                            t8 = t8 + t1 - t0
                            t9 = t9 + dble(toc - tic)/dble(tps)
                        end if
                    end if
                else
                    shift(j) = .false.
                end if
            end do
            !
            rdum = 0
            Tout = 0
            do i = 1, nscat
                is = (i - 1)*lmax
                !
                ! Can use equation 33 of Mackowski94 paper, which is analogous to
                ! equation 44 of Stout02, and avoid matrix multiplication by
                ! exploitint the fact that the trace of a matrix product is
                ! related to the sum of element-wise products:
                !
                !        trace(matmul(transpose(A),B)) = sum(A*B),
                !
                ! where * is element-wise product.
                !
                !do j=1,lmax
                !   rdum = rdum + realpart(Tin(is+j,is+j))
                !enddo
                !
                Mat = 0
                do j = 1, nscat
                    js = (j - 1)*lmax
                    MatJ = 0
                    !MatJ(1:lmax,1:lmax) = Tin(is+1:is+lmax,js+1:js+lmax) ! WRONG!!
                    MatJ(1:pmax, 1:pmax) = Tin(is + 1:is + pmax, js + 1:js + pmax)
                    MatJ(pmax2 + 1:pmax2 + pmax, 1:pmax) = &
                        Tin(is + pmax + 1:is + lmax, js + 1:js + pmax)
                    MatJ(1:pmax, pmax2 + 1:pmax2 + pmax) = &
                        Tin(is + 1:is + pmax, js + pmax + 1:js + lmax)
                    MatJ(pmax2 + 1:pmax2 + pmax, pmax2 + 1:pmax2 + pmax) = &
                        Tin(is + pmax + 1:is + lmax, js + pmax + 1:js + lmax)
                    !
                    if (shift(j)) then
                        if (verb > 1) then
                            call cpu_time(t0)
                            call system_clock(tic)
                        end if
                        call offsetTmat( &
                            Tmat=MatJ, &
                            off=vtacs(:, :, j), &
                            miet=.false., &
                            rtr=rtr, &
                            right=.true., &
                            bigD_=bigD(:, :, j), &
                            useD_=.true.)

                        if (verb > 1) then
                            call cpu_time(t1)
                            call system_clock(toc, count_rate=tps)
                            t2 = t2 + t1 - t0
                            t3 = t3 + dble(toc - tic)/dble(tps)
                        end if

                    end if
                    !
                    Mat = Mat + MatJ
                    !
                    ! if(j > i) then ! Sanity check
                    !    call calcVTACs( &
                    !         r0 = scatXYZR(1:3,j) - scatXYZR(1:3,i), &
                    !         k = cmplx(hostK,0,8), &
                    !         regt = .true., &
                    !         vtacs = matJ )
                    !    rdum = rdum + &
                    !         realpart(sum(transpose(MatJ)*Tin(is+1:is+lmax,js+1:js+lmax)))
                    !    matJ = conjg(matJ) ! Assumed real hostK
                    !    rdum = rdum + &
                    !         realpart(sum(MatJ*Tin(js+1:js+lmax,is+1:is+lmax)))
                    ! endif
                    !
                end do
                !
                ! assuming hostK is real!!!!
                if (shift(i)) then
                    if (verb > 1) then
                        call cpu_time(t0)
                        call system_clock(tic)
                    end if

                    call offsetTmat( &
                        Tmat=Mat, &
                        off=conjg(transpose(vtacs(:, :, i))), &
                        miet=.false., &
                        rtr=rtr, &
                        right=.false., &
                        bigD_=bigD(:, :, i), &
                        useD_=.true.)

                    if (verb > 1) then
                        call cpu_time(t1)
                        call system_clock(toc, count_rate=tps)
                        t2 = t2 + t1 - t0
                        t3 = t3 + dble(toc - tic)/dble(tps)
                    end if
                end if
                !
                Tout = Tout + Mat
                !
            end do
            !
            if (allocated(bigD)) deallocate (bigD)
            !
        end if
        !
        !if(verb > 1) then
        !   call cpu_time(t2)
        !   call system_clock(toc,count_rate=tps)
        !   write(*,'(A,A,2(1x,es10.3e2))') &
        !        myName,'> Calculation time [contraction] (CPU & real in s): ', &
        !        t2-t1, dble(toc-tic)/dble(tps)
        !   call cpu_time(t1)
        !   call system_clock(tic)
        !endif
        !
        if (dumpTmatCol_G) then
            call dumpTmat(tmat=Tout, &
                          filename=tfilename_G, &
                          eps_med=hostEps, &
                          lambda=lambda, &
                          verb_=verb)

            ! also dump the cartesian multipoles for l<=3
            call alphaTensor(Tout, Alpha)
            call dumpTmat(tmat=Alpha, &
                          filename=afilename_G, &
                          eps_med=hostEps, &
                          lambda=lambda, &
                          verb_=verb)

        end if
        call diagnoseTmat(Tmat=Tout, &
                          mode_=1, &
                          verb_=verb)
        !
        if (verb > 1) then
            write (*, '(A,A,2(1x,es10.3e2))') &
                myName, '> Calculation time [offsetTmat] (CPU & real in s): ', &
                t2, t3
            write (*, '(A,A,2(1x,es10.3e2))') &
                myName, '> Calculation time [calcWignerBigD] (CPU & real in s): ', &
                t4, t5
            write (*, '(A,A,2(1x,es10.3e2))') &
                myName, '> Calculation time [calcVTACsAxial] (CPU & real in s): ', &
                t6, t7
            write (*, '(A,A,2(1x,es10.3e2))') &
                myName, '> Calculation time [calcVTACs] (CPU & real in s): ', &
                t8, t9
            write (*, *)
        end if
        !
    end subroutine contractTmat
    !
    subroutine diagnoseTmat(Tmat, mode_, verb_)
        !
        ! Determine the value of n <= nmax when Tr(Re{Tcol}) converges to rtol_G.
        ! If mode_ > 0, also test for the general symmetry, which applies to all
        ! T-matrices. (See equation 5.34 on pg. 121 of Mishchenko's book).
        !
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables:
        complex(8), intent(inout) :: Tmat(:, :)
        integer, intent(in), optional :: mode_, verb_
        ! Local variables
        character(*), parameter :: myName = 'diagnoseTmat'
        integer :: q, n, m, l1, l2, qp, np, mp, lp1, lp2, pmax, nmax, nconv, ij(2), verb, mode
        real(8) :: rdum, rdum2, maxAbs, aveAbs
        real(8), allocatable :: SymErr(:, :), RelSymErr(:, :)
        complex(8) :: tr, trn
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        if (size(Tmat, 1) /= size(Tmat, 2)) then
            write (*, '(A,A)') myName, '> ERROR: size(Tmat,1) /= size(Tmat,2)'
            STOP
        end if
        pmax = size(Tmat, 1)/2
        nmax = int(sqrt(dble(pmax)))
        if (2*nmax*(nmax + 2) /= size(Tmat, 1)) then
            write (*, '(A,A)') myName, '> ERROR: 2*nmax*(nmax+2) /= size(Tmat,1)'
            STOP
        end if
        !
        mode = 0
        if (present(mode_)) mode = mode_
        !
        verb = 0
        if (present(verb_)) verb = verb_
        !
        if (mode > 0) then
            allocate (SymErr(nmax, nmax), RelSymErr(nmax, nmax))
            SymErr = 0; RelSymErr = 0
        end if
        !
        aveAbs = sum(abs(Tmat))/size(Tmat)
        ij = maxloc(abs(Tmat))
        maxAbs = abs(Tmat(ij(1), ij(2)))
        !
        if (verb > 1) then
            write (*, '(A,A,1x,es9.2E3,1x,es11.5E2,A,2(1x,i7))') myName, &
                '> Ave&Max(|Tmat|)=', aveAbs, maxAbs, ' MaxLoc=', ij
            write (*, '(A,A,3(1x,es12.5E2))') myName, &
                '> For i=1:3, |Tmat(i,i)|=', (abs(Tmat(n, n)), n=1, 3)
            write (*, '(A,A,3(1x,es12.5E2))') myName, &
                '> For i=pmax+1:pmax+3, ~=', (abs(Tmat(n + pmax, n + pmax)), n=1, 3)
        end if
        !
        tr = 0
        nconv = 0
        nloop: do n = 1, nmax
            trn = 0
            do q = 1, 2
                ij(1) = (q - 1)*pmax + n*(n + 1)
                do m = -n, n
                    !
                    l1 = ij(1) + m
                    l2 = ij(1) - m
                    !
                    if (mode > 0) then
                        do qp = 1, 2
                            do np = 1, nmax
                                ij(2) = (qp - 1)*pmax + np*(np + 1)
                                do mp = 1, np
                                    lp1 = ij(2) + mp
                                    lp2 = ij(2) - mp
                                    !
                                    rdum = max(abs(Tmat(l1, lp1)), abs(Tmat(l2, lp2)))
                                    if (rdum > 0) then
                                        rdum2 = abs(Tmat(l1, lp1) - (-1)**(m + mp)*Tmat(lp2, l2))
                                        if (rdum2 > symErr(n, np)) symErr(n, np) = rdum2
                                        rdum = rdum2/rdum
                                        if (rdum > RelSymErr(n, np)) RelSymErr(n, np) = rdum
                                    end if
                                    !
                                end do
                            end do
                        end do
                    end if
                    !
                    if (nconv == 0) trn = trn + Tmat(l1, l1)
                    !
                end do
            end do
            if (nconv == 0) then
                tr = tr + trn
                rdum = abs(tr)
                rdum2 = abs(trn)
                if (rdum > 0 .and. rdum2 > 0) then
                    if (rdum2/rdum < rtol_G) nconv = n
                end if
            end if
        end do nloop
        !
        if (nconv == 0) then
            write (*, '(/,A,A,es7.1E2,A,i2,/)') myname, &
                '> WARNING: Tr(Tmat) not converged to ', &
                rtol_G, ' for n <= nmax= ', nmax
        elseif (verb > 1) then
            write (*, '(A,A,es7.1E2,A,i2)') myname, &
                '> Tr(Tmat) converged to ', rtol_G, ' for n= ', nconv
        end if
        !
        if (mode > 0) then
            !
            ij = maxloc(RelSymErr)
            if (verb > 1) write (*, '(A,A,es9.3E2,A,2(1x,i2))') myname, &
                '> Max(|RelSymErr(Tmat)|)= ', &
                RelSymErr(ij(1), ij(2)), ' for n,n''= ', ij
            !
            ! if(RelSymErr(ij(1),ij(2)) > rtol_G) then
            !    write(*,'(A,A,es7.1E2)') myName, &
            !         '> WARNING: Max(|RelSymErr(Tmat)|) > tol= ', rtol_G
            ! endif
            !
            symErr = symErr/maxAbs ! scale by maximum element magnitude
            ij = maxloc(SymErr)
            if (verb > 1) write (*, '(A,A,es9.3E2,A,2(1x,i2))') myname, &
                '> Max(|SymErr(Tmat)|)/Max(|Tmat|)= ', &
                symErr(ij(1), ij(2)), ' for n,n''= ', ij
            !
            if (symErr(ij(1), ij(2)) > rtol_G) then
                write (*, '(/,A,A,es7.1E2,/)') myName, &
                    '> WARNING: Max(|SymErr(Tmat)|)/Max(|Tmat|)= ', &
                    symErr(ij(1), ij(2))
            end if
            !
            deallocate (SymErr, RelSymErr)
            !
        end if
        !
    end subroutine diagnoseTmat
    !
    subroutine calcOaStout(TIJ, scatXYZR, sigOA, cdOA_, verb_, jAbsOA, escat, ehost, vtacsJK)
        ! ============================================================
        ! Calculate the orientationally averaged extinction and
        ! scattering cross-sections defined in equations 44 and 47
        ! of Stout02. The absorption cross-section is then deduced
        ! from the difference (by conservation of energy).
        ! ============================================================
        !
        use swav, only: calcVTACs, calcAbsMat
        !
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables:
        ! Passed variables:
        complex(8), intent(in) :: TIJ(:, :)
        real(8), intent(in) :: scatXYZR(:, :)
        real(8), intent(inout) :: sigOA(3)
        real(8), intent(out), optional :: cdOA_, jAbsOA(0:size(scatXYZR, 2))
        integer, intent(in), optional :: verb_
        complex(8), intent(in), optional :: escat(:)
        real(8), intent(in), optional :: ehost
        complex(8), intent(out), optional :: vtacsJK(size(TIJ, 1), size(TIJ, 2))
        ! Local variables
        character(*), parameter :: myName = 'calcOaStout'
        integer :: i, j, k, l, js, ks, is, ls, nscat
        integer :: lmax, pmax, nmax, tic, toc, tps, verb, im, ip
        real(8) :: pref, cdOA, jsigOA(0:size(scatXYZR, 2), 2), hostK
        real(8) :: t3, t5
        real    :: t0, t1, t2, t4
        complex(8), dimension(size(TIJ, 1)/size(scatXYZR, 2), &
                              size(TIJ, 2)/size(scatXYZR, 2)) :: mat, mat2
        complex(8) :: absMat(size(TIJ, 1), size(TIJ, 2))
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        if (size(TIJ, 1) /= size(TIJ, 2)) then
            write (*, '(A,A)') myName, '> ERROR: size(TIJ,1) /= size(TIJ,2)'
            STOP
        elseif (size(scatXYZR, 1) /= 4) then
            write (*, '(A,A)') myName, '> ERROR: size(scatXYZR,1) /= 4'
            STOP
        end if
        !
        verb = 0
        if (present(verb_)) verb = verb_
        !
        nscat = size(scatXYZR, 2)
        lmax = size(TIJ, 1)/nscat
        pmax = lmax/2
        nmax = int(sqrt(dble(pmax)))
        absMat = 0
        hostK = sigOA(1)
        if (verb > 0) write (*, '(A,A,4(5x,A,4x))') myname, '> J', &
            'csExtOA', 'csScaOA', 'csAbsOA', 'cdExtOA'
        !
        t2 = 0.0d0
        t3 = 0.0d0
        t4 = 0.0d0
        t5 = 0.0d0
        ! Precompute all the regular translation matrices
        do j = 1, nscat
            js = (j - 1)*lmax
            do k = j + 1, nscat
                ks = (k - 1)*lmax

                if (verb > 1) then
                    call cpu_time(t0)
                    call system_clock(tic)
                end if
                call calcVTACs( &
                    r0=scatXYZR(1:3, j) - scatXYZR(1:3, k), &
                    k=cmplx(hostK, 0, kind(hostK)), &
                    regt=.true., &
                    vtacs=vtacsJK(js + 1:js + lmax, ks + 1:ks + lmax))

                vtacsJK(ks + 1:ks + lmax, js + 1:js + lmax) = &
                    conjg(transpose(vtacsJK(js + 1:js + lmax, ks + 1:ks + lmax)))

                if (verb > 1) then
                    call cpu_time(t1)
                    call system_clock(toc, count_rate=tps)
                    t2 = t2 + t1 - t0
                    t3 = t3 + dble(toc - tic)/dble(tps)
                end if
            end do
            vtacsJK(js + 1:js + lmax, js + 1:js + lmax) = 0
            do k = 1, lmax
                vtacsJK(js + k, js + k) = 1
            end do
        end do
        !
        ! Precompute absorption matrix gamma (Stout 2002, eq:49)
        do j = 1, nscat
            js = (j - 1)*lmax

            if (verb > 1) then
                call cpu_time(t0)
                call system_clock(tic)
            end if
            call calcAbsMat( &
                Xi=cmplx(hostK*scatXYZR(4, j), 0, kind(hostK)), &
                ro=sqrt(escat(j)/ehost), &
                mat=absMat(js + 1:js + lmax, js + 1:js + lmax))
            if (verb > 1) then
                call cpu_time(t1)
                call system_clock(toc, count_rate=tps)
                t4 = t4 + t1 - t0
                t5 = t5 + dble(toc - tic)/dble(tps)
            end if

        end do
        pref = tpi/(hostK*hostK)
        jsigOA(0, :) = 0   !the 1st line of this matrix contains sum of the ext. (1st column) and sca. (2nd column)
        ! by all particles and
        !the remaining rows contain ext. and sca. by each particle.
        jAbsOA = 0
        cdOA = 0
        do j = 1, nscat

            js = (j - 1)*lmax
            jsigOA(j, :) = 0
            !
            do k = 1, nscat
                ks = (k - 1)*lmax
                !
                ! First process extinction cross-section
                if (k == j) then
                    mat = TIJ(js + 1:js + lmax, ks + 1:ks + lmax)
                else
                    mat = matmul(TIJ(js + 1:js + lmax, ks + 1:ks + lmax), &
                                 vtacsJK(ks + 1:ks + lmax, js + 1:js + lmax))
                end if
                do im = 1, lmax ! Compute the trace
                    jsigOA(j, 1) = jsigOA(j, 1) + realpart(mat(im, im))
                end do
                do im = 1, pmax
                    cdOA = cdOA + realpart(mat(pmax + im, im) + mat(im, pmax + im)) !CD extinction
                end do
                !
                ! Now process scattering cross-section
                !mat(:,:) = 0
                ! mat2(:,:) = 0
                !do l=1,nscat
                !   ls=(l-1)*lmax
                !   if(l==k) then
                !      mat = mat + TIJ(js+1:js+lmax, ls+1:ls+lmax)
                !   else
                !     mat = mat + matmul(TIJ(js+1:js+lmax, ls+1:ls+lmax), &
                !          vtacsJK(ls+1:ls+lmax, ks+1:ks+lmax))
                !   endif
                ! Wil need mat for absorption cross-section!!!
                !   if(l==j) then
                !     mat2 = mat2 + TIJ(ls+1:ls+lmax, ks+1:ks+lmax)
                !  else
                !     mat2 = mat2 + matmul(vtacsJK(js+1:js+lmax, ls+1:ls+lmax), &
                !          TIJ(ls+1:ls+lmax, ks+1:ks+lmax) ) !mat2 is actually transpose of the real mat2.
                !this is because of the last equation which calculates trace.
                !  endif
                !enddo
                ! Left-multiply mat2 by mat and compute the trace
                !do n=1,lmax
                !   do l=1,lmax
                !      jsigOA(j,2)=jsigOA(j,2)+realpart(conjg(mat(n,l))*mat2(n,l))
                !   enddo
                !enddo
                !****************D Code *************************************
                !****************A code *************************************

                do l = 1, nscat   !***
                    ! Now process absorption cross-section
                    ls = (l - 1)*lmax
                    mat = matmul(conjg(transpose(TIJ(js + 1:js + lmax, ks + 1:ks + lmax))), &
                                 absMat(js + 1:js + lmax, js + 1:js + lmax))

                    if (l == k) then
                        mat2 = TIJ(js + 1:js + lmax, ls + 1:ls + lmax)
                    else
                        mat2 = matmul(TIJ(js + 1:js + lmax, ls + 1:ls + lmax), &
                                      vtacsJK(ls + 1:ls + lmax, ks + 1:ks + lmax))
                    end if
                    mat2 = matmul(mat, mat2)
                    do im = 1, lmax
                        jAbsOA(j) = jAbsOA(j) + realpart(mat2(im, im))
                    end do
                    ! Now process scattering cross-section
                    do i = 1, nscat   !++

                        is = (i - 1)*lmax
                        if (j == k) then
                            mat = conjg(transpose(TIJ(js + 1:js + lmax, ls + 1:ls + lmax)))
                        else
                            mat = matmul(conjg(transpose(TIJ(js + 1:js + lmax, ls + 1:ls + lmax))), &
                                         vtacsJK(js + 1:js + lmax, ks + 1:ks + lmax))
                        end if
                        if (i == l) then
                            mat2 = TIJ(ks + 1:ks + lmax, is + 1:is + lmax)
                        else
                            mat2 = matmul(TIJ(ks + 1:ks + lmax, is + 1:is + lmax), &
                                          vtacsJK(is + 1:is + lmax, ls + 1:ls + lmax))
                        end if

                        do im = 1, lmax
                            do ip = 1, lmax
                                jsigOA(j, 2) = jsigOA(j, 2) + realpart(mat(im, ip)*mat2(im, ip))
                            end do
                        end do

                    end do   !++
                end do  !***
                !**********************************************************************
                !
            end do
            !
            jsigOA(j, 1) = -pref*jsigOA(j, 1) ! Extinction
            jsigOA(j, 2) = pref*jsigOA(j, 2) ! Scattering

            !jAbsOA(j)= jsigOA(j,1)- jsigOA(j,2)  !Absorption
            jAbsOA(j) = pref*jAbsOA(j)  ! Absorption
            jAbsOA(0) = jAbsOA(0) + jAbsOA(j)
            if (verb > 0) write (*, '(13x,i3,3(1x,es15.8E2))') j, jsigOA(j, 1:2), jAbsOA(j)
            jsigOA(0, 1:2) = jsigOA(0, 1:2) + jsigOA(j, 1:2)

        end do
        !
        cdOA = 2*pref*cdOA !
        !
        if (verb > 0) then
            write (*, '(13x,67("-"))')
            write (*, '(13x,A,4(1x,es15.8E2))') 'Sum', &
                jsigOA(0, 1:2), jsigOA(0, 1) - jsigOA(0, 2), cdOA
            write (*, *)
        end if
        !
        sigOA(1:2) = jsigOA(0, 1:2)
        sigOA(3) = jsigOA(0, 1) - jsigOA(0, 2)
        if (present(cdOA_)) cdOA_ = cdOA
        !
        if (verb > 1) then

            write (*, '(A,A,2(1x,es10.3e2))') &
                myName, '> Calculation time [calcVTACs] (CPU & real in s): ', &
                t2, t3
            write (*, '(A,A,2(1x,es10.3e2))') &
                myName, '> Calculation time [calcAbsMat] (CPU & real in s): ', &
                t4, t5
            write (*, *)
        end if
        !
    end subroutine calcOaStout
    !
    !===============================================================
    !
    subroutine applyRotTranzRotOnMat(vtacs, mat, bigdOP, rightOP)
        !
        ! ============================================================
        ! INPUT:
        ! ------
        ! vtacs(2*pmax,2*pmax) - Axial VTACs with (m,n,q) indexing [COMPLEX]
        ! bigdOP(pmax,pmax) - Optional Wigner D-functions for rotation [COMPLEX]
        ! rightOP - option for applying from right [LOGICAL]
        !
        ! IN/OUTPUT:
        ! ------
        ! mat(2*pmax,2*pmax) - matrix [COMPLEX]
        ! ============================================================
        !
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables
        complex(8), intent(in) :: vtacs(:, :) ! (m,n,q) indexing
        complex(8), intent(inout) :: mat(size(vtacs, 1), size(vtacs, 2))
        complex(8), intent(in), optional :: bigdOP(size(vtacs, 1)/2, size(vtacs, 2)/2)
        logical, intent(in), optional :: rightOP
        ! Local variables
        character(*), parameter :: myName = 'applyRotTranzRotOnMat'
        logical :: right, doRot
        integer :: pmax, nmax, n, m, q, l, p, lo, hi, qq!, tic,toc,tps
        complex(8), dimension(size(mat, 1), size(mat, 2)) :: lmat, lvtacs
        complex(8) :: lbigD(size(mat, 1)/2, size(mat, 2)/2)
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        if (present(bigdOP)) then
            doRot = .true.
        else
            doRot = .false.
        end if
        if (present(rightOP)) then
            right = rightOP
        else
            right = .false.
        end if
        !
        ! Check the passed matrix and infer pmax with nmax
        if (size(mat, 1) /= size(mat, 2)) then
            write (*, '(A,A)') myname, '> ERROR: Passed matrix not square'
            STOP
        else
            pmax = size(mat, 1)/2
            nmax = int(sqrt(dble(pmax)))
            n = 2*nmax*(nmax + 2)
            if (n /= size(mat, 1)) then
                if (n < size(mat, 1)) then
                    write (*, '(A,A)') myname, '> WARNING: 2*nmax*(nmax+2) < size(mat,1)'
                else
                    write (*, '(A,A)') myname, '> ERROR: 2*nmax*(nmax+2) > size(vtacs,1)'
                    STOP
                end if
            end if
        end if
        !
        ! Left multiplication by Mat is preferable, because rearranging Mat's
        ! columns is faster than rearranging the rows, when re-indexing from
        ! (q,n,m) to (m,n,q) convention and back.
        if (right) then ! Do Mat*(RotTransIRot)
            lvtacs = vtacs
            if (doRot) lbigD = bigdOP
        else ! Do (RotTransIRot)*Mat = transpose(Mat)*transpose(RotTransIRot)
            mat = transpose(mat)
            lvtacs = transpose(vtacs)
            if (doRot) lbigD = conjg(bigdOP)
        end if
        !
        ! >>>  First rotate >>>
        !call system_clock(tic,tps)
        if (doRot) then
            do q = 0, 1
                qq = q*pmax
                do n = 1, nmax
                    lo = n*n
                    hi = lo + 2*n
                    mat(:, lo + qq:hi + qq) = matmul(mat(:, lo + qq:hi + qq), lbigD(lo:hi, lo:hi))
                end do
            end do
        end if
        !call system_clock(toc,tps)
        !write(*,'(A,A,es15.8E2)') myName,'> 1st rotation t= ',dble(toc-tic)/real(tps)
        ! <<< End of rotation <<<
        !
        ! >>> Then translate axially >>>
        !call system_clock(tic,tps)
        ! Rearrange the columns of mat from (q,n,m) to (m,n,q) indexing
        l = 0
        do m = -nmax, nmax
            do n = max(1, abs(m)), nmax
                do q = 0, 1
                    qq = q*pmax
                    l = l + 1
                    p = n*(n + 1) + m
                    lmat(:, l) = mat(:, p + qq)
                end do
            end do
        end do
        ! With lvtacs and lmat now both index by the (m,n,q) convention,
        ! multiply and exploit the block-diagonal form of lvtacs.
        lo = 1
        do m = -nmax, nmax
            hi = lo + 2*(nmax - max(1, abs(m)) + 1) - 1
            lmat(:, lo:hi) = matmul(lmat(:, lo:hi), lvtacs(lo:hi, lo:hi))
            lo = hi + 1
        end do
        ! Rearrange the columns of lmat from (m,n,q) back to (q,n,m) indexing
        l = 0
        do m = -nmax, nmax
            do n = max(1, abs(m)), nmax
                do q = 0, 1
                    qq = q*pmax
                    l = l + 1
                    p = n*(n + 1) + m
                    mat(:, p + qq) = lmat(:, l)
                end do
            end do
        end do
        !call system_clock(toc,tps)
        !write(*,'(A,A,es15.8E2)') myName,'> Axial translation t= ',dble(toc-tic)/real(tps)
        ! <<< End of axial translation <<<
        !
        ! >>> Finally, rotate back >>>
        !call system_clock(tic,tps)
        if (doRot) then
            do q = 0, 1
                qq = q*pmax
                do n = 1, nmax
                    lo = n*n
                    hi = lo + 2*n
                    mat(:, lo + qq:hi + qq) = matmul(mat(:, lo + qq:hi + qq), &
                                                     conjg(transpose(lbigD(lo:hi, lo:hi))))
                end do
            end do
        end if
        !call system_clock(toc,tps)
        !write(*,'(A,A,es15.8E2)') myName,'> 2nd rotation t= ',dble(toc-tic)/real(tps)
        ! <<< End of rotation <<<
        !
        if (.not. right) then
            mat = transpose(mat)
        end if
        !
    end subroutine applyRotTranzRotOnMat
    !
    !===============================================================
    !

    subroutine calcField(r, geometry, ipwVec, ipwE0, scaCJ, &
                         intCJreg_, intCJirr_, scatK_, reE, imE, reB, imB, &
                         reE_sca, imE_sca, reB_sca, imB_sca, p_label, dipoles, verb_)
        !
        use swav, only: calcVSWs
        !
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables
        real(8), intent(in)  :: r(:, :), geometry(:, :), ipwVec(3)
        complex(8), intent(in) :: scaCJ(:), ipwE0(3)
        complex(8), intent(in), optional :: intCJirr_(size(scaCJ), 4), &
                                            intCJreg_(size(scaCJ), 4), scatK_(size(geometry, 2), 4)
        real(8), intent(out) :: reE(size(r, 1), size(r, 2)), imE(size(r, 1), size(r, 2))
        real(8), intent(out), optional :: reE_sca(size(r, 1), size(r, 2)), imE_sca(size(r, 1), size(r, 2))
        real(8), intent(out), optional :: reB_sca(size(r, 1), size(r, 2)), imB_sca(size(r, 1), size(r, 2))
        real(8), intent(out), optional :: reB(size(r, 1), size(r, 2)), imB(size(r, 1), size(r, 2))
        real(8), intent(in), optional :: dipoles(:, :)
        integer, intent(inout) :: p_label(:, :)
        integer, intent(in), optional :: verb_
        ! Local variables
        character(*), parameter :: myName = 'calcField'
        integer :: i, j, k, js, npts, nscat, lmax, idepth
        integer :: tic, toc, tps, verb
        real :: t1, t2
        real(8) :: rj(3), xyz(3), rdum, scatRot(3, 3, size(geometry, 2)), r2, hostK
        complex(8) :: vsw(3, size(scaCJ)/size(geometry, 2)), E(3), E_sca(3), B_sca(3)
        complex(8) :: vswb(3, size(scaCJ)/size(geometry, 2)), B(3)

        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------

        if (size(r, 1) /= 3) then
            write (*, '(A,A)') myName, '> ERROR: size(r,1) /= 3'
            STOP
        elseif (size(geometry, 1) /= 8) then
            write (*, '(A,A)') myName, '> ERROR: size(geometry,1) /= 8'
            STOP
        end if
        npts = size(r, 2)
        nscat = size(geometry, 2)
        lmax = size(scaCJ)/nscat
        hostK = sqrt(dot_product(ipwVec, ipwVec))
        i = int(sqrt(dble(lmax/2)))
        if (2*i*(i + 2) /= lmax) then
            write (*, '(A,A)') myName, '> ERROR: Incommensurate lmax= ', lmax
            STOP
        end if
        !
        ! First compute the rotation matrices for spheroidal scatterers.
        ! The matrices are used for checking if a point is inside a spheroid.
        do j = 1, nscat
            if (geometry(8, j) > tiny1) then
                scatRot(:, :, j) = matmul(RotMatY(-geometry(6, j)), & ! \theta -> \beta
                                          RotMatZ(-geometry(7, j)))  ! \phi -> \gamma
            else
                scatRot(:, :, j) = 0
                do i = 1, 3
                    scatRot(i, i, j) = 1
                end do
            end if
        end do
        verb = 0
        if (present(verb_)) verb = verb_
        call cpu_time(t1)
        call system_clock(tic)
        !if(verb > 1) write(*,'(A,A,i8)') myName,'> Number of grid-points: ',npts
        !
        ! NOTE: scaling ~npts*nscat, not just ~npts.
        points: do i = 1, npts
            ! this routine's output should be verbosity dependent
            ! if(verb > 1) write(*,*)'r=', r(1:3,i)
            ! Check if point is inside a scatterer
            idepth = 0

            jscat1: do j = 1, nscat
                js = (j - 1)*lmax
                r2 = geometry(4, j)*geometry(4, j)
                rj(1:3) = r(1:3, i) - geometry(1:3, j)
                rdum = dot_product(rj, rj)
                if (rdum < r2) then ! Inside smallest circumscribing sphere
                    if (geometry(8, j) > tiny1) then ! check for spheroid
                        ! check if inside the "focal" sphere
                        ! if(geometry(8,j) >= 1) then
                        !    ! prolate spheroid, rugby ball
                        !    if(rdum < r2*(1.0d0 - 1.0d0/geometry(8,j)**2)) idepth = 1
                        ! else
                        !    ! oblate spheroid, pumpkin
                        !    if(rdum < r2*(1 - geometry(8,j)**2)) idepth = 1
                        ! endif
                        if (idepth /= 1) then
                            xyz(:) = matmul(scatRot(:, :, j), rj(:))
                            xyz(1:2) = xyz(1:2)*geometry(8, j) ! scale x,y by aspect ratio
                            rdum = dot_product(xyz, xyz)
                            if (rdum < r2) then
                                idepth = 1
                                p_label(i, 1) = j
                                p_label(i, 2) = 0
                            end if
                        end if
                    else ! inside sphere
                        idepth = -1
                        p_label(i, 1) = j
                        shells: do k = 1, abs(nint(geometry(8, j))) ! check shells
                            if (rdum < geometry(4 + k, j)**2) then
                                p_label(i, 2) = idepth
                                idepth = idepth - 1

                            else
                                exit shells
                            end if
                        end do shells
                    end if
                    if (abs(idepth) > 0) exit jscat1
                end if
            end do jscat1
            !
            if (abs(idepth) > 0) then
                E = 0
                E_sca = 0
                B = 0
                B_sca = 0
                if (present(intCJreg_) .and. present(intCJirr_) .and. present(scatK_) .and. idepth < 0) then
                    ! write(*,*)'idepth1=',idepth ! too verbose
                    ! Re-use rj, j, and js from jscat1 loop
                    call calcVSWs( &
                        r=rj, &
                        k=scatK_(j, abs(idepth)), &
                        pmax=lmax/2, &
                        regt=.true., &
                        cart=.true., &
                        waves=vsw, &
                        wavesB=vswb)
                    E = E + matmul(vsw, intCJreg_(js + 1:js + lmax, abs(idepth)))
                    B = B + matmul(vswb, intCJreg_(js + 1:js + lmax, abs(idepth)))
                    E_sca(1:3) = 0
                    B_sca(1:3) = 0
                    if (abs(idepth) /= 1 + abs(nint(geometry(8, j)))) then
                        ! write(*,*)'idepth2=',idepth ! too verbose
                        call calcVSWs( &
                            r=rj, &
                            k=scatK_(j, abs(idepth)), &
                            pmax=lmax/2, &
                            regt=.false., &
                            cart=.true., &
                            waves=vsw, &
                            wavesB=vswb)

                        E = E + matmul(vsw, intCJirr_(js + 1:js + lmax, abs(idepth)))
                        B = B + matmul(vswb, intCJirr_(js + 1:js + lmax, abs(idepth)))
                        E_sca(1:3) = 0
                        B_sca(1:3) = 0
                    end if
                end if
            else ! the point is out of the particle
                E = ipwE0*exp(cmplx(0, dot_product(ipwVec, r(:, i)), kind(ipwVec)))  !incident field
                E_sca = 0
                B_sca = 0
                B(1) = ipwVec(2)*ipwE0(3) - ipwVec(3)*ipwE0(2)   !B=(k x E)/omega, I omit omega
                B(2) = ipwVec(3)*ipwE0(1) - ipwVec(1)*ipwE0(3)
                B(3) = ipwVec(1)*ipwE0(2) - ipwVec(2)*ipwE0(1)
                B = B*exp(cmplx(0, dot_product(ipwVec, r(:, i)), kind(ipwVec)))
                jscat2: do j = 1, nscat
                    js = (j - 1)*lmax
                    rj(1:3) = r(1:3, i) - geometry(1:3, j)
                    call calcVSWs( &
                        r=rj, &
                        k=cmplx(hostK, 0, kind(hostK)), &
                        pmax=lmax/2, &
                        regt=.false., &
                        cart=.true., &
                        waves=vsw, &
                        wavesB=vswb)
                    E(:) = E(:) + matmul(vsw, scaCJ(js + 1:js + lmax))
                    E_sca(:) = E_sca(:) + matmul(vsw, scaCJ(js + 1:js + lmax))
                    B(:) = B(:) + matmul(vswb, scaCJ(js + 1:js + lmax))  !Note: I omitted 1/omega in B, I add this in mapNF
                    B_sca(:) = B_sca(:) + matmul(vswb, scaCJ(js + 1:js + lmax))
                end do jscat2

            end if

            do j = 1, 3
                reE(j, i) = realpart(E(j))
                imE(j, i) = imagpart(E(j))
            end do

            if (present(reB)) then
                do j = 1, 3
                    reB(j, i) = realpart(B(j))
                    imB(j, i) = imagpart(B(j))
                end do
            end if
            if (present(reE_sca)) then
                do j = 1, 3
                    reE_sca(j, i) = realpart(E_sca(j))
                    imE_sca(j, i) = imagpart(E_sca(j))
                end do

            end if
            if (present(reB_sca)) then
                do j = 1, 3
                    reB_sca(j, i) = realpart(B_sca(j))
                    imB_sca(j, i) = imagpart(B_sca(j))
                end do

            end if

        end do points

        call cpu_time(t2)
        call system_clock(toc, count_rate=tps)
        if (verb > 1) then
            write (*, '(A,A,2(1x,es10.3e2))') &
                myName, '> Calculation time (CPU & real in s): ', &
                t2 - t1, dble(toc - tic)/dble(tps)
            write (*, *)
        end if
        !
    end subroutine calcField
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    subroutine dumpTmat(tmat, filename, lambda, eps_med, tol_, verb_)
        !
        ! Routine for dumping the collective T-matrix to file
        ! in the format: s, s', n, n', m, m', T_re, T_im
        !
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables:
        complex(8), intent(in), dimension(:, :) :: tmat
        character(*), intent(in) :: filename
        real(8), intent(in) :: lambda, eps_med
        real(8), intent(in), optional :: tol_
        integer, intent(in), optional :: verb_
        ! Local variables:
        character(*), parameter :: myName = 'dumpTmat'
        integer :: s, n, m, sp, np, mp, nmax, i, ip, verb, lunit
        real(8) :: tol
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        if (size(tmat, 1) /= size(tmat, 2)) then
            write (*, '(A,A)') myname, '> ERROR: Matrix tmat not square'
            STOP
        end if
        !
        verb = 0
        if (present(verb_)) verb = verb_
        !
        nmax = int(sqrt(dble(size(tmat, 1)/2)))
        if (2*nmax*(nmax + 2) /= size(tmat, 1)) then
            write (*, '(A,A)') myname, '> ERROR: Incommensurate nmax!'
            STOP
        end if
        !
        tol = 0
        if (present(tol_)) tol = tol_*maxval(abs(tmat))
        !
        inquire (file=filename, number=lunit)
        if (lunit == -1) then
            lunit = tunit
            open (lunit, file=filename, status='replace')
            write (lunit, '(A)') '# s sp n np m mp Tr Ti'
        elseif (lunit /= tunit) then
            write (*, '(A,A,A,A)') myname, '> ERROR: File ', trim(filename), &
                ' connected to wrong unit'
            STOP
        end if
        !
        write (lunit, '(A,f8.2,A,i9,A,es15.8)') '# lambda= ', lambda, &
            ' nelements= ', size(tmat, 1)**2, ' eps_med= ', eps_med
        !
        ip = 0
        do sp = 1, 2
            do np = 1, nmax
                do mp = -np, np
                    ip = ip + 1
                    i = 0
                    do s = 1, 2
                        do n = 1, nmax
                            do m = -n, n
                                i = i + 1
                                !
                                ! Ideally would like to use tolerance, but that
                                ! means we don't kow nelements in advance... so
                                ! print entire matrix (of known size) for now.
                                !
                                !if( abs(tmat(i,ip)) > tol ) then
                                write (lunit, &
                                       '(2(1x,i1),2(1x,i2),2(1x,i3),2(1x,es22.15))') &
                                    s, sp, n, np, m, mp, &
                                    realpart(tmat(i, ip)), imagpart(tmat(i, ip))
                                !endif
                                !
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !
    end subroutine dumpTmat
    !
    subroutine dumpMatrix(mat, ofile, tolOP, verb_)
        !
        !============================================================
        ! Print matrix M, with i, i', Re(M_ii'), Im(M_ii') on each
        ! printed line. The first line contains i_{max} and i'_{max}
        ! in the first two fields, and then an optional tolerance,
        ! 'tolOP', to discard small values.
        !------------------------------------------------------------
        ! For a scatterer-centred T-matrix, i and i' are composite
        ! indices, each spanning the usual multipolatity indices
        ! n <= n_{max} and m(n), index q=1,2 for M_{nm} and N_{nm}
        ! waves, and also the scatterer index j <= j_{max}.
        !
        ! More precisely:
        !
        ! i(j,k) = (j-1)*l_{max} + l <= i_{max}(j_{max},l_{max})
        ! l(q,p) = (q-1)*p_{max} + p <= l_{max}(q_{max},p_{max})
        ! p(n,m) = n*(n+1)+m <= p_{max}(n_{max}) = n_{max}*(n_{max}+2)
        !
        ! so 1 <= i <= i_{max} = j_{max}*l_{max}
        !                      = 2*j_{max}*p_{max}
        !                      = 2*j_{max}*n_{max}*(n_{max}+2).
        !
        ! To unpack i need only j_{max} and n_{max}, then use:
        !
        ! j = j_{max}*(k-1)/i_{max} + 1, note integer division,
        ! l = k - (j-1)*i_{max}/j_{max},
        ! q = 2*(l-1)/l_{max} + 1, again note integer division,
        ! p = l - (q-1)*l_{max}/2,
        ! n = floor(sqrt(p)) or int(sqrt(p)),
        ! m = p - n*(n+1).
        !
        ! Note that for j_{max}=1 we have i=l, corresponding to a
        ! single T-matrix centred about the origin.
        !============================================================
        !
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables:
        complex(8), intent(in), dimension(:, :) :: mat
        character(*), intent(in), optional :: ofile
        real(8), intent(in), optional :: tolOP
        integer, intent(in), optional :: verb_
        ! Local variables:
        character(*), parameter :: myName = 'dumpMatrix'
        integer, parameter :: lunit = 777
        integer :: i, j, imax, jmax, nzeroes, verb
        real(8) :: tol, t0, t1
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        verb = 0
        if (present(verb_)) verb = verb_
        !
        if (verb > 1) call cpu_time(t0)
        !
        imax = ubound(mat, 1)
        jmax = ubound(mat, 2)
        if (present(ofile)) then
            open (lunit, file=ofile, status='replace')
            write (lunit, '(2(1x,i6))') imax, jmax
        else
            write (*, '(2(1x,i6))') imax, jmax
        end if
        nzeroes = 0
        if (present(tolOP)) then
            tol = tolOP*maxval(abs(mat))
        else
            tol = 0.0d0
        end if
        do i = 1, imax
            do j = 1, jmax
                if (abs(mat(i, j)) <= tol) then
                    nzeroes = nzeroes + 1
                else
                    if (present(ofile)) then
                        write (lunit, '(2(1x,i6),2(1x,es23.15E3))') &
                            i, j, mat(i, j)
                    else
                        write (*, '(2(1x,i6),2(1x,es23.15E3))') &
                            i, j, mat(i, j)
                    end if
                end if
            end do
        end do
        if (present(ofile)) close (lunit)
        if (verb > 1) then
            call cpu_time(t1)
            write (*, '(A,A,i12,A,es10.3E2,A)') &
                myname, '> Wrote ', imax*jmax - nzeroes, &
                ' non-zero values in ', t1 - t0, ' seconds.'
        end if
        !
    end subroutine dumpMatrix
    !
    subroutine offsetTmat(Tmat, off, miet, rtr, right, bigD_, useD_, balJI_)
        !
        ! Offset the supplied <Tmat> by <off>, which can be either a square matrix
        ! of VTACs or a (complex!) displacement vector kr(3) from which VTACs will
        ! be generated. Regular or irregular VTACs will be generated depending on
        ! whether kr(3) is pureley real or purely imaginary, respectively, and the
        ! wavenumber for medium will be treated as a real number. If logical miet
        ! is true, Tmat will be treated as diagonal. If logical rtr is true, then
        ! offsetting will be based on factorised translation. If logical right is
        ! true, then offsetting will done by post-multiplying Tmat (from right).
        !
        ! NB: balJI_ triggers balancing of the VTACs and the T-matrix individually,
        ! before the offsetting, but currently works only without factorised trans-
        ! lation.
        !
        use swav, only: calcVTACs, calcWignerBigD, calcVTACsAxial, xyz2rtp
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables:
        complex(8), intent(inout) :: Tmat(:, :)
        complex(8), intent(in) :: off(:, :)
        logical, intent(in) :: miet, rtr, right
        complex(8), intent(in), optional :: bigD_(size(Tmat, 1)/2, size(Tmat, 2)/2)
        logical, intent(in), optional :: useD_
        integer, intent(in), optional :: balJI_(2)
        ! Local variables:
        character(*), parameter :: myName = 'offsetTmat'
        logical :: regular, applyRot, flip, balance
        integer :: lmax, pmax, nmax, regtest(3), i, j, n, m
        real(8) :: rdum, xyz(3), rtp(3), k
        complex(8), allocatable, dimension(:, :) :: vtacs, bigD
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        !
        ! First check that Tmat is square
        if (size(Tmat, 1) /= size(Tmat, 2)) then
            write (*, '(A,A)') myName, '> ERROR: Tmat not square'
            STOP
        else
            lmax = size(Tmat, 1)
            pmax = lmax/2
            nmax = int(sqrt(dble(pmax)))
            if (2*nmax*(nmax + 2) /= lmax) then
                write (*, '(A,A)') myName, '> ERROR: 2*nmax*(nmax+2) /= lmax'
                STOP
            end if
        end if
        !
        !write(*,*) 'TEST: rtr,right=', rtr,right
        if (verb_G > 1) then
            rdum = 0
            nloop1: do n = 1, nmax
                do i = 1, 2
                    j = (i - 1)*pmax + n*(n + 1)
                    do m = -n, n
                        k = abs(Tmat(j + m, j + m))
                        rdum = rdum + k
                        if (k > 0 .and. rdum > 0) then
                            if (k/rdum < rtol_G) then
                                exit nloop1
                            end if
                        end if
                    end do
                end do
            end do nloop1
            write (*, '(A,A,es8.1E3)', advance='NO') &
                myName, '> Initial Max(|Tmat|)= ', maxval(abs(Tmat))
            if (n <= nmax) then
                write (*, '(A,es8.1E3,A,i2)') ' Tr(|Tmat|)= ', rdum, ' nconv= ', n
            else
                write (*, '(A,es8.1E3,A,i2)') ' Tr(|Tmat|)= ', rdum, &
                    ' unconverged!'
            end if
        end if
        !
        balance = .false.
        if (present(balJI_)) then
            if (balJI_(1) > 0 .and. balJI_(2) > 0) then
                if (rtr) then
                    write (*, '(A,A)') myName, &
                        '> ERROR: Cannot balance with factorised translation!'
                    STOP
                end if
                balance = .true.
                if (right) then
                    i = balJI_(1)
                else
                    i = balJI_(2)
                end if
                call balanceMatJI( &
                    j=i, jregt=.false., &
                    Mat=Tmat, &
                    i=i, iregt=.true.)
            end if
        end if
        !
        ! Process the offset matrix <off> and do the offsetting
        if (size(off, 1) == 3 .and. size(off, 2) == 1) then
            !
            ! Generate VTACs from <off>, which we now process to infer "regularity"
            allocate (vtacs(lmax, lmax))
            k = sqrt(dot_product(off(:, 1), off(:, 1)))
            if (k < 1.0d-99) return ! negligible offset magnitude
            do i = 1, 3
                rdum = abs(off(i, 1))
                if (rdum/k > tiny) then
                    if (abs(imagpart(off(i, 1)))/rdum < tiny) then
                        regtest(i) = 1
                        xyz(i) = realpart(off(i, 1))/k
                    elseif (abs(realpart(off(i, 1)))/rdum < tiny) then
                        regtest(i) = -1
                        xyz(i) = imagpart(off(i, 1))/k
                    else
                        write (*, '(A,A,i1)') &
                            myName, '> ERROR: Non-zero Re and Im parts for i= ', i
                        STOP
                    end if
                else
                    xyz(i) = 0
                    regtest(i) = 0
                end if
                do j = 1, i - 1
                    if (regtest(i)*regtest(j) == -1) then
                        write (*, '(A,A)') &
                            myName, '> ERROR: Re and Im components present!'
                        STOP
                    end if
                end do
            end do
            i = sum(regtest)
            if (i > 0) then
                regular = .true.
            elseif (i < 0) then
                regular = .false.
            else
                write (*, '(A,A)') myName, '> ERROR: Failed to infer regularity!'
                STOP
            end if
            !
            ! Compute and apply the rotation and axial translation matrices
            if (rtr) then
                !
                call xyz2rtp(xyz, rtp)
                !
                if ((abs(rtp(2)) < tiny1 .or. abs(rtp(2) - pi) < tiny1) .and. &
                    (abs(rtp(3)) < tiny1 .or. abs(rtp(3) - tpi) < tiny1)) then
                    applyRot = .false.
                    ! In the absence of rotation, must check the sign of xyz(3)
                    ! and assign the axial vtacs_mN and vtacs_Nm correctly, to
                    ! make sure the intended axial translation is consistent
                    ! with the z-axis direction.
                    if (xyz(3) < 0.0d0) then
                        flip = .true.
                    else
                        flip = .false.
                    end if
                else
                    applyRot = .true.
                    ! With rotation, the sign of dr(3) is irrelevan since the
                    ! rotation will always align the local axes so that dr points
                    ! in the +ve z-direction, so no need to flip!
                    flip = .false.
                    allocate (bigD(pmax, pmax))
                    call calcWignerBigD( &
                        angles=(/rtp(3), rtp(2), 0.0d0/), &
                        pmax=pmax, &
                        bigD=bigD)
                end if
                !
                call calcVTACsAxial( &
                    r0=rtp(1), &
                    k=cmplx(k, 0, kind(k)), &
                    pmax=pmax, &
                    regt=regular, &
                    flip=flip, &
                    vtacs=vtacs)
                !
                if (applyRot) then
                    call applyRotTranzRotOnMat( &
                        vtacs=vtacs, &
                        mat=Tmat, &
                        bigdOP=bigD, & ! Include rotation matrices
                        rightOP=right)
                else
                    call applyRotTranzRotOnMat( &
                        vtacs=vtacs, &
                        mat=Tmat, &
                        rightOP=right)
                end if
                !
                if (allocated(bigD)) deallocate (bigD)
                !
            else
                !
                !write(*,*) 'TEST2'
                !
                !i=1
                !if(right) i=-1
                !
                call calcVTACs( &
                    r0=xyz, &
                    k=cmplx(k, 0, kind(k)), &
                    regt=regular, &
                    vtacs=vtacs)
                !
                if (balance) call balanceMatJI( &
                    j=balJI_(1), jregt=.true., &
                    Mat=vtacs, &
                    i=balJI_(2), iregt=.false.)
                !
                if (miet) then
                    if (right) then
                        do j = 1, lmax ! Right-multiply diagonal T-matrix
                            Tmat(j, :) = Tmat(j, j)*vtacs(j, :)
                        end do
                    else
                        do i = 1, lmax ! Left-multiply diagonal T-matrix
                            Tmat(:, i) = vtacs(:, i)*Tmat(i, i)
                        end do
                    end if
                else
                    if (right) then
                        Tmat = matmul(Tmat, vtacs)
                    else
                        Tmat = matmul(vtacs, Tmat)
                    end if
                end if
                !
            end if
            !
            if (allocated(vtacs)) deallocate (vtacs)
            !
        elseif (size(off, 1) == lmax .and. size(off, 2) == lmax) then
            !
            ! VTACs supplied in off
            if (rtr) then
                !
                if (present(bigD_) .and. present(useD_)) then ! not f2py-safe
                    applyRot = useD_
                else
                    applyRot = .false.
                end if
                !
                if (applyRot) then
                    call applyRotTranzRotOnMat( &
                        vtacs=off, &
                        mat=Tmat, &
                        bigdOP=bigD_, & ! Include rotation matrices
                        rightOP=right)
                else
                    call applyRotTranzRotOnMat( &
                        vtacs=off, &
                        mat=Tmat, &
                        rightOP=right)
                end if
                !
            else
                !
                allocate (vtacs(lmax, lmax))
                vtacs = off ! Recall: off has intent in
                !
                if (balance) call balanceMatJI( &
                    j=balJI_(1), jregt=.true., &
                    Mat=vtacs, &
                    i=balJI_(2), iregt=.false.)
                !
                if (miet) then
                    if (right) then
                        do j = 1, lmax
                            do i = 1, lmax ! Right-multiply diagonal T-matrix
                                Tmat(i, j) = Tmat(i, i)*vtacs(i, j) ! was off(i,j)
                            end do
                        end do
                    else
                        do i = 1, lmax ! Left-multiply diagonal T-matrix
                            Tmat(:, i) = vtacs(:, i)*Tmat(i, i) ! was off(:,i)
                        end do
                    end if
                else
                    if (right) then
                        Tmat = matmul(Tmat, vtacs) ! was off
                    else
                        Tmat = matmul(vtacs, Tmat) ! was off
                    end if
                end if
                !
                deallocate (vtacs)
                !
            end if
            !
        else
            write (*, '(A,A)') myName, '> ERROR: Badly shaped off'
        end if
        !
        !call diagnoseTmat(Tmat, verb_= 1)
        !
        if (verb_G > 1) then
            rdum = 0
            nloop2: do n = 1, nmax
                do i = 1, 2
                    j = (i - 1)*pmax + n*(n + 1)
                    do m = -n, n
                        k = abs(Tmat(j + m, j + m))
                        rdum = rdum + k
                        if (k > 0 .and. rdum > 0) then
                            if (k/rdum < rtol_G) then
                                exit nloop2
                            end if
                        end if
                    end do
                end do
            end do nloop2
            write (*, '(A,A,es8.1E3)', advance='NO') &
                myName, '>   Final Max(|Tmat|)= ', maxval(abs(Tmat))
            if (n <= nmax) then
                write (*, '(A,es8.1E3,A,i2)') ' Tr(|Tmat|)= ', rdum, ' nconv= ', n
            else
                write (*, '(A,es8.1E3,A,i2)') ' Tr(|Tmat|)= ', rdum, &
                    ' unconverged!'
            end if
            ! SymTest >>>
            ! do n=1,nmax
            !    do i=1,2
            !       j = (i-1)*pmax + n*(n+1) ! p
            !       rdum = maxval(abs(Tmat(j-n:j+n,j-n:j+n)))*rtol_G
            !       do m=1,n
            !          if(abs(Tmat(j+m,j+m)-Tmat(j-m,j-m)) > rdum) then
            !             write(*,'(A,A,es7.1E2,A,3(1x,i2))') &
            !                  myname,'> WARNING: |RelSymErr| > tol= ',rtol_G,&
            !                  ' for q,n,|m|=',i,n,m
            !             return
            !          endif
            !       enddo
            !    enddo
            ! enddo
            ! <<< SymTest
            !
        end if
        !
    end subroutine offsetTmat
    !
    !function RotMatX(ang) result(rotMat)
    !  real(8), intent(in) :: ang
    !  real(8) :: rotMat(3,3),c,s
    !  c = cos(ang); s = sin(ang)
    ! rotMat(1,1) = 1; rotMat(1,2) = 0; rotMat(1,3) = 0
    ! rotMat(2,1) = 0; rotMat(2,2) = c; rotMat(2,3) =-s
    ! rotMat(3,1) = 0; rotMat(3,2) = s; rotMat(3,3) = c
    ! return
    !end function RotMatX
    !
    function RotMatY(ang) result(rotMat)
        real(8), intent(in) :: ang
        real(8) :: rotMat(3, 3), c, s
        c = cos(ang); s = sin(ang)
        rotMat(1, 1) = c; rotMat(1, 2) = 0; rotMat(1, 3) = s
        rotMat(2, 1) = 0; rotMat(2, 2) = 1; rotMat(2, 3) = 0
        rotMat(3, 1) = -s; rotMat(3, 2) = 0; rotMat(3, 3) = c
        return
    end function RotMatY
    !
    function RotMatZ(ang) result(rotMat)
        real(8), intent(in) :: ang
        real(8) :: rotMat(3, 3), c, s
        c = cos(ang); s = sin(ang)
        rotMat(1, 1) = c; rotMat(1, 2) = -s; rotMat(1, 3) = 0
        rotMat(2, 1) = s; rotMat(2, 2) = c; rotMat(2, 3) = 0
        rotMat(3, 1) = 0; rotMat(3, 2) = 0; rotMat(3, 3) = 1
        return
    end function RotMatZ
    !
    subroutine readTmatFile(Tmat, filename, unit, wavelen,HDF5_in_, verb_)
        
        use HDFfive, only: h5_rd_file, h5_rd_vec
        !
        complex(8), intent(inout) :: Tmat(:, :)
        character(*), intent(in) :: filename
        integer, intent(in) :: unit
        real(8), intent(in) :: wavelen
        integer, intent(in), optional :: verb_
        logical, intent(in), optional :: HDF5_in_
        !
        character(*), parameter :: myName = 'readTmatFile'
        character(64) :: comment
        logical :: yes
        integer :: i, j, eof, m, mp, n, np, q, qp, k, ij(1)
        integer :: lmax, pmax, nmax, nelements, verb
        real(8) :: x, y
        complex(8), allocatable :: Tmat_h5(:, :, :)
        real(8), allocatable :: ang_wave_num(:), wavelen_h5(:), ldum(:), mdum(:)  
        integer, allocatable :: pdum(:), sdum(:), qdum(:)     
        !
        verb = 0
        if (present(verb_)) verb = verb_
        !
        if (verb > 1) write (*, '(A,A,A)', advance='NO') &
            myName, '> Filename ', trim(filename)
        !
        lmax = size(Tmat, 1)
        print *, 'lmax'
        print *, lmax
        pmax = lmax/2
        nmax = int(sqrt(dble(pmax)))
        if (size(Tmat, 2) /= lmax) then
            write (*, '(A,A)') myname, '> ERROR: Passed Tmat not square!'
            STOP
        elseif (lmax /= 2*nmax*(nmax + 2)) then
            write (*, '(A,A)') myname, '> ERROR: lmax /= 2*nmax*(nmax+2)'
            STOP
        end if
        !
        Tmat = 0
        !
        inquire (file=trim(filename), exist=yes)
        if (.not. yes) then
            write (*, '(/,A,A)') myname, '> ERROR: Missing file ', trim(filename)
            STOP
        end if
        inquire (file=trim(filename), number=i)
        if (i == -1) then
            open (unit, file=trim(filename), status='old')
            read (unit, *, iostat=eof) ! comment with column labels
        elseif (i /= unit) then
            write (*, '(A,A,i5)') myName, '> ERROR: File connected to unit i /= ', unit
            STOP
        end if
        !
        read (unit, '(A)', iostat=eof) comment ! comment line
        !
        if (.not. HDF5_in_)then
              i = index(comment, '=')  ! look for 1st equal
              comment = adjustl(comment(i + 1:))
              read (comment, *, iostat=eof) x ! will read just the wavelength
              if (eof /= 0) then
                 write (*, '(/,A,A)') myName, &
                   '> ERROR: Failed to extract wavelength from comment'
                 STOP
              elseif (abs(x - wavelen) > tiny1) then
                  write (*, '(/,A,A)') myName, '> ERROR: Wavelengh mismatch'
                  STOP
              elseif (verb > 1) then
                  write (*, '(A,es15.8E2)') '  wavelength=', x
              end if
        !
              i = index(comment, '=')  ! look for 2nd equal
              comment = adjustl(comment(i + 1:))
              read (comment, *, iostat=eof) nelements ! read the expected element count
              if (eof /= 0) then
                  write (*, '(A,A)') myName, &
                  '> ERROR: Failed to extract element count from comment'
                  STOP
              end if
        !
              kloop: do k = 1, nelements
                  read (unit, *, iostat=eof) q, qp, n, np, m, mp, x, y
                  if (eof /= 0) then
                      write (*, '(A,A,A)') myName, '> ERROR: Premature end of file ', &
                          trim(filename)
                      STOP
                  elseif (q < 1 .or. qp < 1) then
                      write (*, '(A,A)') myName, '> ERROR: q < 1, wrong column order?'
                      STOP
                  end if
                  if (n > nmax .or. np > nmax) cycle kloop
                  i = (q - 1)*pmax + n*(n + 1) + m
                  j = (qp - 1)*pmax + np*(np + 1) + mp
                  Tmat(i, j) = cmplx(x, y, kind(Tmat))
                  end do kloop
          else
                 print *, filename
           	 call h5_rd_vec(filename,'/','vacuum_wavelength', wavelen_h5)
            	!  if (allocated(wavelen_h5)) deallocate(wavelen_h5)
            	!  allocate(wavelen_h5(size(ang_wave_num)))
            	!  wavelen_h5 = tpi/ang_wave_num
            	 ij=minloc(abs(wavelen_h5-wavelen), 1)           	          	  
          	 !----------------------------------------------------
          	 call h5_rd_vec(filename, '/modes','l', ldum)
                 call h5_rd_vec(filename, '/modes','m', mdum)    
                 allocate(pdum(size(ldum)), sdum(size(ldum)), qdum(size(ldum)) )   
                   pdum = int(ldum)*(int(ldum) + 1) + int(mdum)
                   sdum(1:size(sdum):2)=2
                   sdum(2:size(sdum):2)=1
                   qdum=(sdum-1)*maxval(pdum)+pdum 
                   if (allocated(Tmat_h5)) deallocate(Tmat_h5)  	      
            	   call h5_rd_file(filename, '/','tmatrix', Tmat_h5)
            	   !------------New: making Tmatrix compatible with TERMS -------------------
            	  ! allocate(tmat(size(Tmat_h5,1), size(Tmat_h5,2), size(Tmat_h5,3)))
            	   Tmat =0
            	   do i=1,size(Tmat,1)
            	   	do j=1, size(Tmat,1)
            		   Tmat(qdum(j),qdum(i))=Tmat_h5(i,j,ij(1))
            		end do
            	   end do 
          	 !----------------------------------------------------   
          end if
        !
    end subroutine readTmatFile
    !
    subroutine parseInc(inc, inc_dirn, inc_ampl, verb_)
        !
        real(8), intent(in) :: inc(:)
        real(8), intent(out) :: inc_dirn(3)
        complex(8), intent(out), optional :: inc_ampl(3)
        integer, intent(in), optional :: verb_
        !
        character(*), parameter :: myName = 'parseInc'
        real(8), parameter :: dirn(3) = (/0, 0, 1/)
        real(8), parameter :: isqrt2 = 1/sqrt(2.0d0)
        complex(8), parameter :: ampl(3, 4) = &
                                 !
                                 reshape((/ &
                                         !
                                         ! Linear polarization along x
                                         !
                                         cmplx(1, 0, kind(isqrt2)), &
                                         cmplx(0, 0, kind(isqrt2)), &
                                         cmplx(0, 0, kind(isqrt2)), &
                                         !
                                         ! Linear polarization along y
                                         !
                                         cmplx(0, 0, kind(isqrt2)), &
                                         cmplx(1, 0, kind(isqrt2)), &
                                         cmplx(0, 0, kind(isqrt2)), &
                                         !
                                         ! Right circular polarization (R) convention of the
                                         ! receiver (looking at beam coming to you, rotating
                                         ! clockwise). Taking the conjugate will produce the
                                         ! left circular polarisation (L).
                                         !
                                         cmplx(isqrt2, 0, kind(isqrt2)), &
                                         cmplx(0, -isqrt2, kind(isqrt2)), &
                                         cmplx(0, 0, kind(isqrt2)), &
                                         !
                                         ! Left circular polarization (L)
                                         cmplx(isqrt2, 0, kind(isqrt2)), &
                                         cmplx(0, isqrt2, kind(isqrt2)), &
                                         cmplx(0, 0, kind(isqrt2)) &
                                         !
                                         ! indices: components, polarisation id.
                                         /), (/3, 4/))
        !
        integer :: pol_id, verb
        real(8) :: angles(3) ! Euler angles for ZY'Z' rotation
        !
        verb = 0
        if (present(verb_)) verb = verb_
        !
        angles = 0
        !
        if (size(inc) > 2) angles(1:2) = inc(2:3)

        inc_dirn(1:3) = matmul(rotZYZmat(angles(1:3)), dirn(1:3))

        if (present(inc_ampl)) then
            pol_id = nint(inc(1))

            if (pol_id < 1 .or. pol_id > 4) then
                write (*, '(A,A)') myname, '> ERROR: Bad pol_id= ', pol_id
                STOP
            end if
            if (size(inc) > 3) angles(3) = inc(4)
            inc_ampl(1:3) = matmul(rotZYZmat(angles(1:3)), &
                                   ampl(1:3, pol_id))

            if (verb > 1) then
                write (*, '(A,A,3(9x,A))') myName, '>', 'dirn', 'Re(ampl)', 'Im(ampl)'
                write (*, '(10x,A,3(1x,es15.8))') &
                    'x:', inc_dirn(1), realpart(inc_ampl(1)), imagpart(inc_ampl(1))
                write (*, '(10x,A,3(1x,es15.8))') &
                    'y:', inc_dirn(2), realpart(inc_ampl(2)), imagpart(inc_ampl(2))
                write (*, '(10x,A,3(1x,es15.8))') &
                    'z:', inc_dirn(3), realpart(inc_ampl(3)), imagpart(inc_ampl(3))
            end if
        else
            if (verb > 1) then
                write (*, '(A,A)') 'parseSca >', 'dirn'
                write (*, '(10x,A,1(1x,es15.8))') 'x:', inc_dirn(1)
                write (*, '(10x,A,1(1x,es15.8))') 'y:', inc_dirn(2)
                write (*, '(10x,A,3(1x,es15.8))') 'z:', inc_dirn(3)

            end if
        end if
        !
    end subroutine parseInc

!--------------------------------------------------------------------------------------------
    !
    function rotZYZmat(angles) result(mat)
        !
        ! Euler angles (alpha,beta,gamma) for ZY'Z' rotation
        real(8), intent(in) :: angles(3)
        real(8) :: c(3), s(3), c1c2, s1s3, c3s1, mat(3, 3)
        !
        ! matrix R(alpha,beta,gamma) = Rz(alpha)Ry(beta)Rz(gamma)
        ! represents extrinsic rotations about the global, space-fixed
        ! axes z,y, and z by angles gamma, beta, and alpha, respectively
        ! (in that order!); or intrinsic rotations about the local
        ! moving axes z, y', z' by angles alpha, beta, gamma (in that
        ! order!). Confusing...
        !
        c = cos(angles)
        s = sin(angles)
        !
        c1c2 = c(1)*c(2)
        s1s3 = s(1)*s(3)
        c3s1 = c(3)*s(1)
        !
        mat(1, 1) = c1c2*c(3) - s1s3
        mat(2, 1) = c(1)*s(3) + c(2)*c3s1
        mat(3, 1) = -c(3)*s(2)
        !
        mat(1, 2) = -c3s1 - c1c2*s(3)
        mat(2, 2) = c(1)*c(3) - c(2)*s1s3
        mat(3, 2) = s(2)*s(3)
        !
        mat(1, 3) = c(1)*s(2)
        mat(2, 3) = s(1)*s(2)
        mat(3, 3) = c(2)
        !
        return
        !
    end function rotZYZmat

    !-------------------------------------------------------------------------------------------------------

    subroutine calcStokesScaVec(sca_angles, inc2, ncut, wavelen, ehost, geometry, scheme, &
                                tfiles_, escat_, nselect_, noRTR_, verb_, &
                                StokesPhaseMat, StokesScaVec, diff_sca, HDF5_in)

        use sphmsv, only: calcScatMat, calcStokesIncVec, calcStokesPhaseMat

        !-----------------------------------------------------------------------------
    !!This subroutine calculates the Stokes Scattering vector based on the formula 2.102
        !of Mishchenko's book on page 51.
        !---------------------------------------------------
        ! Start of variable declarations
        !---------------------------------------------------
        ! Passed variables

        integer, intent(in) :: scheme, ncut(3)
        real(8), intent(in) :: wavelen(:)
        real(8), intent(in) :: ehost(size(wavelen))
        real(8), intent(in) :: geometry(:, :)
        real(8), intent(in)     :: inc2(:, :), sca_angles(:, :)
        complex(8), intent(in), optional :: &
            escat_(size(geometry, 2), 4, size(wavelen))
        ! NOTE: f2py fails for character(*) argument, so chose 256 chars
        character*256, intent(in), optional :: tfiles_(size(geometry, 2))
        integer, intent(in), optional :: nselect_(2, 2, 2, size(geometry, 2)), verb_
        logical, intent(in), optional :: noRTR_
        real(8), intent(inout) :: StokesScaVec(:, :, :, :), diff_sca(:, :, :, :)
        real(8), intent(inout) :: StokesPhaseMat(:, :, :, :)
	logical, intent(in) :: HDF5_in(:)
        ! Local variables
        character(*), parameter :: myName = 'calcStokesScaVec'

        integer :: nmax, lmax, nscat, nwav
        character*256 :: tfiles(size(geometry, 2))
        logical :: scatMiet(size(geometry, 2)), ldum
        logical ::  noRTR
        real(8) :: hostK
        integer :: iwav, verb, j, js, i, is, ierr, isca, idum
        integer :: nselect(2, 2, 2, size(geometry, 2)), units(size(geometry, 2))
        complex(8), allocatable :: TIJ(:, :), Tcol(:, :), cJ(:, :, :), TIJ_(:, :)

        real(8)  :: StokesinVec(4, size(inc2, 2))
        !real(8), dimension(4,4) :: StokesPhaseMat
        complex(8), dimension(2, 2) :: ScatMat
        complex(8), allocatable :: ipwAmpl(:, :)
        real(8), allocatable  :: ipwDirn(:, :), spwDirn(:, :)
        real :: t0, t1
        integer :: tic, toc, tps
        !------------------------------------------------
        nscat = size(geometry, 2)
        nmax = ncut(1)
        lmax = 2*nmax*(nmax + 2)
        nwav = size(wavelen)
        rtol_G = 10.0d0**ncut(3)
        !
        !---------------------------------------------------
        ! Start processing arguments
        !---------------------------------------------------
        allocate (ipwDirn(3, size(inc2, 2)), ipwAmpl(3, size(inc2, 2)))
        allocate (spwDirn(3, size(sca_angles, 2)))

        ! if (size(wavelen) > 1) then
        !   allocate(StokesScaVec(4,1,size(wavelen),size(sca_angles,2)))
        !  allocate(StokesPhaseMat(4,4,size(wavelen),size(sca_angles,2)))
        ! else
        !   allocate(StokesScaVec(4,1,size(inc2,2),size(sca_angles,2)))
        !  allocate(StokesPhaseMat(4,4,size(inc2,2),size(sca_angles,2)))
        ! endif

        if (size(geometry, 1) /= 8) then
            write (*, *) myName, '> ERROR: size(geometry,1) /= 8'
            STOP
        end if
        !
        tfiles = ''; units = 0
        if (present(tfiles_)) then
            do j = 1, nscat
                if (len_trim(tfiles_(j)) > 0) then
                    tfiles(j) = tfiles_(j)
                    units(j) = 1000 + j
                    iloop: do i = 1, j - 1
                        ! Check for filename duplicates and, if present,
                        ! store scatterer ID with the same T-matrix file.
                        ! Store -ve of the ID to distinguish from actual
                        ! file handles in "units".
                        if (tfiles(j) == tfiles(i)) then
                            units(j) = -i
                            exit iloop
                        end if
                    end do iloop
                end if
            end do
        end if
        !
        if (maxval(geometry(8, :)) >= tiny1) then
            ! non-Mie scatterer(s) present
            if (maxval(units) <= 0) then
                write (*, *) myName, &
                    '> ERROR: Missing tfiles for non-Mie scatterers'
                STOP
            end if
        end if
        !
        if (minval(geometry(8, :)) < tiny1) then
            ! Mie scatterer(s) present
            ldum = .false.
            if (present(escat_)) then
                if (maxval(abs(escat_)) > tiny2) then
                    ldum = .true. ! non-zero escat_ present
                end if
            end if
            if (.not. ldum) then
                write (*, *) myName, '> ERROR: Missing non-zero escat_'
                STOP
            end if
        end if
        !
        noRTR = .false.
        if (present(noRTR_)) noRTR = noRTR_
        !
        j = ncut(2)
        j = 2*j*(j + 2)
        allocate (TIJ(nscat*lmax, nscat*lmax), Tcol(j, j))
        allocate (TIJ_(nscat*lmax, nscat*lmax))
        !TIJ=0
        !
        nselect = 0 ! f2py will default to zero
        if (present(nselect_)) nselect = nselect_
        !
        verb = 0
        if (present(verb_)) verb = verb_

        !---------------------------------------------------------------
        ! Done processing/initialising arguments. Start calculations.
        !---------------------------------------------------------------
        allocate (cJ(nscat*lmax, 2, 1))
        waves: do iwav = 1, nwav
            !
            if (verb > 0) then
                write (*, *)
                write (*, '(A,A,f8.2,A)') myname, '> ===== Wavelength: ', &
                    wavelen(iwav), ' (nm) ======================'
                write (*, *)
            end if
            if (verb > 1) then
                call cpu_time(t0)
                call system_clock(tic)
            end if
            !
            TIJ = 0
            TIJ_ = 0
            !
            do j = 1, nscat
                js = (j - 1)*lmax
                if (geometry(8, j) < tiny1) then ! Mie scatterer
                    scatMiet(j) = .true.
                    do i = 1, 1 - nint(geometry(8, j))
                        TIJ(js + i, js + i) = escat_(j, i, iwav)
                    end do
                else
                    scatMiet(j) = .false.
                    if (units(j) > 0) then
                        call readTmatFile( &
                            Tmat=TIJ(js + 1:js + lmax, js + 1:js + lmax), &
                            filename=trim(tfiles(j)), &
                            unit=units(j), &
                            wavelen=wavelen(iwav), &
                            HDF5_in_ = HDF5_in(j))
                    elseif (units(j) < 0) then
                        i = abs(units(j))
                        is = (i - 1)*lmax
                        TIJ(js + 1:js + lmax, js + 1:js + lmax) = &
                            TIJ(is + 1:is + lmax, is + 1:is + lmax)
                    else
                        write (*, '(A,A,i5)') &
                            myName, '> ERROR: units(j)=0 for j=', j
                        STOP
                    end if
                end if
            end do

            !if(seek_fi) then
            !write(*,*)'size(inc2,2)',size(inc2,2)
            !write(*,*)'size(sca_angles,2)',size(sca_angles,2)
            do idum = 1, size(inc2, 2)
                TIJ_ = TIJ
                cJ = 0
                cJ(1:4, 1, 1) = inc2(:, idum)
                Tcol = 0
                call solve( &
                    wavelen=wavelen(iwav), &
                    ehost=ehost(iwav), &
                    geometry=geometry, &
                    TIJ=TIJ_, &
                    scheme_=scheme, &
                    nselect_=nselect, &
                    verb_=verb, &
                    cJ_=cJ, &
                    noRTR_=noRTR, &
                    ierr_=ierr)

                ! wavenumber k_host
                hostK = tpi/wavelen(iwav)*sqrt(ehost(iwav))
                Tcol(1, 1) = cmplx(wavelen(iwav), ehost(iwav), kind(Tcol))
                j = lmax*nscat
                if (scheme == 3) then
                    j = lmax
                    ldum = .true. ! need this for single scatterer!
                else
                    ldum = .false.
                end if

                call contractTmat( &
                    Tin=TIJ_(:, 1:j), &
                    scatXYZR=geometry(1:4, :), &
                    Tout=Tcol, &
                    rtr=.not. noRTR, &
                    verb_=verb, &
                    mack_=ldum)

                call parseInc( &
                    inc=inc2(:, idum), &
                    inc_dirn=ipwDirn(:, idum), &
                    inc_ampl=ipwAmpl(:, idum), &
                    verb_=verb_)
                call calcStokesIncVec(ehost(iwav), ipwDirn(:, idum), ipwAmpl(:, idum), StokesinVec(:, idum), verb_)

                do isca = 1, size(sca_angles, 2)
                    !write(*,*)'sca_angles(:,isca)',sca_angles(:,isca)
                    call parseInc( &
                        inc=(/0.0d0, sca_angles(:, isca)/), &
                        inc_dirn=spwDirn(:, isca), &
                        verb_=verb_)

                    call calcScatMat( &
                        Tcol, &
                        hostK, &  !k=cmplx(hostK,0,kind(hostK)), &
                        spwDirn(:, isca), &
                        ipwDirn(:, idum), &
                        ScatMat, &
                        verb_)
                    if (iwav > 1) then
                        call calcStokesPhaseMat(ScatMat, StokesPhaseMat(:, :, iwav, isca), verb_)
                    StokesScaVec(1:4, 1, iwav, isca) = (0.5d0)*sqrt(ehost(iwav)*eps0/mu0)*matmul(StokesPhaseMat(:, :, iwav, isca), &
                                                                                                     StokesinVec(:, idum))
                        diff_sca(1, 1, iwav, isca) = StokesScaVec(1, 1, iwav, isca)/StokesinVec(1, idum)
                    else
                        call calcStokesPhaseMat(ScatMat, StokesPhaseMat(:, :, idum, isca), verb_)
                    StokesScaVec(1:4, 1, idum, isca) = (0.5d0)*sqrt(ehost(iwav)*eps0/mu0)*matmul(StokesPhaseMat(:, :, idum, isca), &
                                                                                                     StokesinVec(:, idum))
                        diff_sca(1, 1, idum, isca) = StokesScaVec(1, 1, idum, isca)/StokesinVec(1, idum)

                    end if

                end do
            end do !idum
            if (verb > 1) then
                call cpu_time(t1)
                call system_clock(toc, count_rate=tps)
                write (*, '(A,A,2(1x,es10.3E2))') &
                    myName, '> Calculation time (CPU & real in s): ', &
                    t1 - t0, dble(toc - tic)/dble(tps)
            end if
        end do waves

        if (verb > 0) write (*, '(A,A,/)') myname, '> Done!'

        if (allocated(TIJ)) deallocate (TIJ)
        if (allocated(Tcol)) deallocate (Tcol)
        if (allocated(cJ)) deallocate (cJ)
        !if(allocated(inc2)) deallocate(inc2)
    end subroutine calcStokesScaVec
    !----------------------------------------------------------------

    subroutine calcLDOC(Ef, Bf, N_OpC, verb_)

        !-----------------------------------------------------------------
        ! This subroutine calculates the optical chirality and normalized optical
        ! chirality relative to the oc of RHCP light.
        ! OC of RHCP is :omega*eps0*|E|^2 /(2*c) , c:speed of light,
        ! OC= -omega*eps0/2 *Im(E*.B)
        !
        ! note: omega is removed with the omega in B's formula
        !-----------------------------------------------------------------
        ! real(8), intent(in) :: hostk
        real(8), intent(in) :: Ef(:, :, :), Bf(:, :, :)
        real(8), intent(out) :: N_OpC(1, 1, size(Ef, 3))
        integer, intent(in), optional :: verb_
        ! Local variables
        character(*), parameter :: myName = 'calcLDOC'
        integer :: ngrid, idum
        real(8), parameter :: eps0 = 8.8541878128d-12

        !real(8), parameter :: c= 299792458
        real :: t0, t1
        integer :: tic, toc, tps, verb
        complex(8) :: e(3), b(3)
        !----------------------------------------------------------------
        verb = 0
        if (present(verb_)) verb = verb_

        if (verb > 1) then
            call cpu_time(t0)
            call system_clock(tic)
        end if
        ngrid = size(Ef, 3)
        do idum = 1, ngrid
            e = cmplx(Ef(:, 2, idum), Ef(:, 3, idum), kind(Ef))
            b = cmplx(Bf(:, 2, idum), Bf(:, 3, idum), kind(Bf))
            N_OpC(1, 1, idum) = dble(imagpart(DOT_PRODUCT(e, b))) !dot prodduct for complex vector is dot(dconjg(a),b)

        end do
        ! OpC_=dble(-0.5d0*eps0*N_OpC)
        N_OpC = dble(-N_OpC/hostK_G)
        if (verb > 1) then
            call cpu_time(t1)
            call system_clock(toc, count_rate=tps)
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time (CPU & real in s): ', &
                t1 - t0, dble(toc - tic)/dble(tps)

        end if
        if (verb > 0) write (*, '(A,A,/)') myname, '> Done!'

    end subroutine calcLDOC
!--------------------------------------------------------------------------------------------------------

!******************************************************************************
!Average of OC2
!******************************************************************************

    subroutine calcOaLDOC(pol_type, r, geometry, TIJ, Or_OC, lambda, &
                          ehost, p_label, scatK_, verb_)

        !------------------------------------------------------------------
        ! This subroutine calculates the orientation average of the local degree of chirality
        ! based on Phys. Rev. B 103, 115405 (2021) http://dx.doi.org/10.1103/PhysRevB.103.115405
        ! notations follow the TERMS user guide
        !============================================================
        use swav, only: calcVSWs, calcVTACs, xyz2rtp, calcVTrtp2xyz, calcLamMat
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables:
        integer, intent(in) :: pol_type
        real(8), intent(in), dimension(:, :) :: r
        !real(8), intent(in)  ::  hostK
        integer, intent(in), optional :: verb_
        complex(8), intent(in) :: TIJ(:, :)
        complex(8), intent(in), optional ::  scatK_(:, :)
        real(8), intent(in) :: geometry(:, :), ehost, lambda
        real(8), intent(out), optional :: Or_OC(size(r, 2), 10)
        integer, intent(inout) :: p_label(:, :)

        ! Local variables:
        character(*), parameter :: myName = 'calcOaLDOC'
        integer :: i, j, k, l, np, is, js, ks, ls, pmax, nmax, nscat, npts, lmax, n, idepth, ndum
        integer :: tic, toc, tps, verb, ip
        real    :: t0, t1, t2, t4, t6
        real(8) :: rp(3), t3, t5, t7
        real(8) :: xyz(3), r2, rj(3), rnp(3), rdum, scatRot(3, 3, size(geometry, 2))
        complex(8) ::  B0R, D0R, B0L, D0L, tr
        complex(8), allocatable :: t11(:, :), t22(:, :), t33(:, :), t44(:, :)

        complex(8), allocatable :: vtacsKL(:, :), si_reg(:, :, :), si_ireg(:, :, :)
        complex(8), allocatable :: ZR_reg(:, :, :), ZR_ireg(:, :, :)
        complex(8), allocatable :: ZL_reg(:, :, :), ZL_ireg(:, :, :), vtacsKL_h(:, :)
        complex(8), allocatable :: TIJ_h(:, :)
        real(8), allocatable :: scatXYZR(:, :)
        real(8), parameter :: isqrt2 = 1/sqrt(2.0d0)
        logical :: RCP = .true., LP = .true.
        complex(8) :: lam_mat(size(TIJ, 1), size(TIJ, 2))
        complex(8), allocatable :: VR(:, :, :, :), UR(:, :, :, :), VL(:, :, :, :), UL(:, :, :, :), &
                                   VpUR(:, :, :, :), VmUR(:, :, :, :), &
                                   VpUL(:, :, :, :), VmUL(:, :, :, :)
        complex(8), parameter :: Alpha = cmplx(1.0d0, 0.0d0, kind(isqrt2)), &
                                 beta = cmplx(0.0d0, 0.0d0, kind(isqrt2))

        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------

        nscat = size(geometry, 2)
        npts = size(r, 2)

        Or_OC = 0
        if (size(TIJ, 1) /= size(TIJ, 2)) then
            write (*, '(A,A)') myName, '> ERROR: TIJ not square'
            STOP
        else
            lmax = size(TIJ, 1)/nscat
            pmax = lmax/2
            nmax = int(sqrt(dble(pmax)))
            if (2*nmax*(nmax + 2) /= lmax) then
                write (*, '(A,A)') myName, '> ERROR: 2*nmax*(nmax+2) /= lmax'
                STOP
            end if
        end if
        verb = 0
        if (present(verb_)) verb = verb_

        if (pol_type == 3) then
            LP = .false.
        elseif (pol_type == 4) then
            LP = .false.
            RCP = .false.
        end if
        allocate (scatXYZR(4, nscat), vtacsKL(size(TIJ, 1), size(TIJ, 2)))
        allocate (si_reg(3, lmax, nscat), si_ireg(3, lmax, nscat))
        allocate (ZR_reg(3, pmax, nscat), ZR_ireg(3, pmax, nscat))
        allocate (ZL_reg(3, pmax, nscat), ZL_ireg(3, pmax, nscat))
        allocate (VR(3, pmax, nscat, nscat), UR(3, pmax, nscat, nscat), &
                  VL(3, pmax, nscat, nscat), UL(3, pmax, nscat, nscat), &
                  VpUR(3, pmax, nscat, nscat), VmUR(3, pmax, nscat, nscat), &
                  VpUL(3, pmax, nscat, nscat), VmUL(3, pmax, nscat, nscat))
        allocate (TIJ_h(size(TIJ, 1), size(TIJ, 2)), &
                  vtacsKL_h(size(TIJ, 1), size(TIJ, 2)))
        scatXYZR(1:4, :) = geometry(1:4, :)
        vtacsKL = 0
        vtacsKL_h = 0
        do j = 1, nscat
            js = (j - 1)*lmax
            do l = 1, nscat
                ls = (l - 1)*lmax
                TIJ_h(js + 1:js + pmax, ls + 1:ls + pmax) =  &  !TIJ_hLL According to Eq. 23 of PRB paper
                 & TIJ(js + 1:js + pmax, ls + 1:ls + pmax) + TIJ(js + pmax + 1:js + lmax, ls + 1:ls + pmax) + &
                 & TIJ(js + 1:js + pmax, ls + 1 + pmax:ls + lmax) + TIJ(js + 1 + pmax:js + lmax, ls + 1 + pmax:ls + lmax)

                TIJ_h(js + 1:js + pmax, ls + 1 + pmax:ls + lmax) =  & !TIJ_hLR
                 & TIJ(js + 1:js + pmax, ls + 1:ls + pmax) + TIJ(js + pmax + 1:js + lmax, ls + 1:ls + pmax) - &
                 & TIJ(js + 1:js + pmax, ls + 1 + pmax:ls + lmax) - TIJ(js + 1 + pmax:js + lmax, ls + 1 + pmax:ls + lmax)

                TIJ_h(js + pmax + 1:js + lmax, ls + 1:ls + pmax) =  & !TIJ_hRL
                 & TIJ(js + 1:js + pmax, ls + 1:ls + pmax) - TIJ(js + pmax + 1:js + lmax, ls + 1:ls + pmax) + &
                 & TIJ(js + 1:js + pmax, ls + 1 + pmax:ls + lmax) - TIJ(js + 1 + pmax:js + lmax, ls + 1 + pmax:ls + lmax)

                TIJ_h(js + pmax + 1:js + lmax, ls + pmax + 1:ls + lmax) =  &  !TIJ_hRR
                 & TIJ(js + 1:js + pmax, ls + 1:ls + pmax) - TIJ(js + pmax + 1:js + lmax, ls + 1:ls + pmax) - &
                 & TIJ(js + 1:js + pmax, ls + 1 + pmax:ls + lmax) + TIJ(js + 1 + pmax:js + lmax, ls + 1 + pmax:ls + lmax)
                !-------------------------------------
            end do
        end do
        TIJ_h = 0.5d0*TIJ_h

        !compute translation
        if (verb > 1) then
            call cpu_time(t0)
            call system_clock(tic)
        end if

        do k = 1, nscat
            ks = (k - 1)*lmax
            do l = k + 1, nscat
                ls = (l - 1)*lmax
                call calcVTACs( &
                    r0=scatXYZR(1:3, k) - scatXYZR(1:3, l), &
                    k=cmplx(hostK_G, 0, kind(hostK_G)), &
                    regt=.true., &
                    vtacs=vtacsKL(ks + 1:ks + lmax, ls + 1:ls + lmax))

                vtacsKL(ls + 1:ls + lmax, ks + 1:ks + lmax) = &  !It is not needed to calculate here! I just calculate its helicity
                    dconjg(transpose(vtacsKL(ks + 1:ks + lmax, ls + 1:ls + lmax)))

                !-------------calculation of VtacsKL_h :just diagonal blocks are needed for D0R and D0L
                VtacsKL_h(ks + 1:ks + pmax, ls + 1:ls + pmax) = &  !VtacsKL_hLL According to Eq. 23 of PRB paper
                     & vtacsKL(ks + 1:ks + pmax, ls + 1:ls + pmax) + vtacsKL(ks + pmax + 1:ks + lmax, ls + 1:ls + pmax) + &
                    & vtacsKL(ks + 1:ks + pmax, ls + 1 + pmax:ls + lmax) + vtacsKL(ks + 1 + pmax:ks + lmax, ls + 1 + pmax:ls + lmax)

                VtacsKL_h(ks + pmax + 1:ks + lmax, ls + pmax + 1:ls + lmax) =  &  !VtacsKL_hRR
                    & vtacsKL(ks + 1:ks + pmax, ls + 1:ls + pmax) - vtacsKL(ks + pmax + 1:ks + lmax, ls + 1:ls + pmax) - &
                    & vtacsKL(ks + 1:ks + pmax, ls + 1 + pmax:ls + lmax) + vtacsKL(ks + 1 + pmax:ks + lmax, ls + 1 + pmax:ls + lmax)

                !-------------
                ! vtacsKL_h(ls+1:ls+pmax, ks+1:ks+pmax) = &
                !   & vtacsKL(ls+1:ls+pmax, ks+1:ks+pmax)+ vtacsKL(ls+pmax+1:ls+lmax, ks+1:ks+pmax)+ &
                !   & vtacsKL(ls+1+pmax:ls+lmax, ks+1:ks+pmax) + vtacsKL(ls+1+pmax:ls+lmax, ks+1+pmax:ks+lmax)

                vtacsKL_h(ls + 1:ls + pmax, ks + 1:ks + pmax) = dconjg(transpose( &     !just because the corresponding off
                  VtacsKL_h(ks + 1:ks + pmax, ls + 1:ls + pmax)))                  !diagonal blocks of vtacsKL are conjg(transpose))

                ! vtacsKL_h(ls+pmax+1:ls+lmax, ks+pmax+1:ks+lmax) = & !VtacsKL_hRR
                ! & vtacsKL(ls+1:ls+pmax, ks+1:ks+pmax) - vtacsKL(ls+1:ls+pmax, ks+pmax+1:ks+lmax)- &
                !    & vtacsKL(ls+1+pmax:ls+lmax, ks+1:ks+pmax)+ vtacsKL(ls+1+pmax:ls+lmax, ks+1+pmax:ks+lmax)

                vtacsKL_h(ls + pmax + 1:ls + lmax, ks + pmax + 1:ks + lmax) = dconjg(transpose( &
                                                                       VtacsKL_h(ks + pmax + 1:ks + lmax, ls + pmax + 1:ls + lmax)))
                !---------------------------------------------------

            end do
            vtacsKL(ks + 1:ks + lmax, ks + 1:ks + lmax) = 0   !diagonal block of vtacsKL
            do l = 1, lmax
                vtacsKL(ks + l, ks + l) = 1
            end do
            !------------------------------
            VtacsKL_h(ks + 1:ks + lmax, ks + 1:ks + lmax) = &  !According to Eq. 23 of PRB paper, The diagonal block of vtacsKL are Identity matrix
                2.0d0*vtacsKL(ks + 1:ks + lmax, ks + 1:ks + lmax) !and if one put it in Eq. 23, arrive at this formula.
            !----------------------------
        end do

        VtacsKL_h = 0.5d0*VtacsKL_h

        if (verb > 1) then
            call cpu_time(t1)
            call system_clock(toc, count_rate=tps)
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [calcVTACs] (CPU & real in s): ', &
                t1 - t0, dble(toc - tic)/dble(tps)
        end if
        !-------------------------------------
        allocate (t11(pmax, 3), &
                  t22(pmax, pmax), t33(3, pmax), t44(pmax, pmax))
        t2 = 0.0d0 ! timing for calcVSWs
        t3 = 0.0d0 ! timing for calcVSWs
        t4 = 0.0d0 ! timing for regular calcVSWs
        t5 = 0.0d0 ! timing for regular calcVSWs
        t6 = 0.0d0 ! timing for irregular calcVSWs
        t7 = 0.0d0 ! timing for irregular calcVSWs
        points: do n = 1, npts

            rp = r(1:3, n)
            !-----------------------------check for points outside particles
            idepth = 0
            np = 0
            np_scat1: do np = 1, nscat
                r2 = geometry(4, np)*geometry(4, np)
                rnp(1:3) = rp - geometry(1:3, np)
                rdum = dot_product(rnp, rnp)

                if (rdum < r2) then   ! Inside smallest circumscribing sphere

                    if (geometry(8, np) > tiny1) then ! check for spheroid
                        ! First compute the rotation matrices for spheroidal scatterers.
                        ! The matrices are used for checking if a point is inside a spheroid.
                        scatRot(:, :, np) = matmul(RotMatY(-geometry(6, np)), & ! \theta -> \beta
                                                   RotMatZ(-geometry(7, np)))  ! \phi -> \gamma
                        xyz(:) = matmul(scatRot(:, :, np), rnp(:))
                        xyz(1:2) = xyz(1:2)*geometry(8, np) ! scale x,y by aspect ratio
                        rdum = dot_product(xyz, xyz)
                        if (rdum < r2) then
                            idepth = 1
                            exit np_scat1
                        end if

                    else ! inside sphere
                        idepth = -1
                        p_label(n, 1) = np
                        ! shells: do k=1,abs(nint(geometry(8,j))) ! check shells
                        !       if(rdum < geometry(4+k,j)**2) then
                        !           p_label(n,2)=idepth
                        !          idepth = idepth - 1

                        !       else
                        !          exit shells
                        !       endif
                        !  enddo shells
                        exit np_scat1
                    end if
                end if
            end do np_scat1

!------------------------------------------------------------

            if (idepth == 0) then

                ! restart times for calcVSWs
                if (verb > 1) then
                    call cpu_time(t0)
                    call system_clock(tic)
                end if

                call calcVSWs( &
                    r=rp, &
                    k=cmplx(hostK_G, 0, kind(hostK_G)), &
                    pmax=lmax/2, &
                    regt=.true., &
                    cart=.true., &
                    waves=si_reg(:, :, 1))

                ! sum times for calcVSWs
                if (verb > 1) then
                    call cpu_time(t1)
                    call system_clock(toc, count_rate=tps)
                    t2 = t2 + t1 - t0
                    t3 = t3 + dble(toc - tic)/dble(tps)
                end if

                do j = 1, nscat
                    rj(1:3) = rp - geometry(1:3, j)
                    ! restart times for calcVSWs
                    if (verb > 1) then
                        call cpu_time(t0)
                        call system_clock(tic)
                    end if

                    call calcVSWs( &
                        r=rj, &
                        k=cmplx(hostK_G, 0, kind(hostK_G)), &
                        pmax=lmax/2, &
                        regt=.false., &
                        cart=.true., &
                        waves=si_ireg(:, :, j))

                    ! sum times for irr calcVSWs
                    if (verb > 1) then
                        call cpu_time(t1)
                        call system_clock(toc, count_rate=tps)
                        t4 = t4 + t1 - t0
                        t5 = t5 + dble(toc - tic)/dble(tps)
                    end if

                    ZR_ireg(:, 1:pmax, j) = isqrt2*(si_ireg(:, 1:pmax, j) - si_ireg(:, pmax + 1:lmax, j)) !Eq. 22a in PRB paper
                    ZL_ireg(:, 1:pmax, j) = isqrt2*(si_ireg(:, 1:pmax, j) + si_ireg(:, pmax + 1:lmax, j)) !Eq. 22b in PRB paper

                    ! restart times for reg calcVSWs
                    if (verb > 1) then
                        call cpu_time(t0)
                        call system_clock(tic)
                    end if

                    call calcVSWs( &
                        r=rj, &
                        k=cmplx(hostK_G, 0, kind(hostK_G)), &
                        pmax=lmax/2, &
                        regt=.true., &
                        cart=.true., &
                        waves=si_reg(:, :, j))

                    ! sum times for reg calcVSWs
                    if (verb > 1) then
                        call cpu_time(t1)
                        call system_clock(toc, count_rate=tps)
                        t6 = t6 + t1 - t0
                        t7 = t7 + dble(toc - tic)/dble(tps)
                    end if

                    ZR_reg(:, 1:pmax, j) = isqrt2*(si_reg(:, 1:pmax, j) - si_reg(:, pmax + 1:lmax, j)) !Eq. 22a in PRB paper
                    ZL_reg(:, 1:pmax, j) = isqrt2*(si_reg(:, 1:pmax, j) + si_reg(:, pmax + 1:lmax, j)) !Eq. 22b in PRB paper

                end do

!---------------------------------------------------
                B0R = 0
                B0L = 0

                do j = 1, nscat
                    js = (j - 1)*lmax
                    do l = 1, nscat
                        ls = (l - 1)*lmax

                        if (RCP .OR. (LP)) then
                            ! U_R = Z_R * T_RR aka TIJ_h(js+pmax+1:js+lmax,ls+pmax+1:ls+lmax)
                            call ZGEMM('N', 'N', 3, pmax, pmax, Alpha, &
                                       ZR_ireg(:, :, j), 3, TIJ_h(js + pmax + 1:js + lmax, ls + pmax + 1:ls + lmax), pmax, &
                                       beta, UR(:, :, j, l), 3)
                            !  Z^h [-Ujl]
                            call calcTrace('C', 'N', ZR_reg(:, :, l), -UR(:, :, j, l), tr)

                            B0R = B0R + (2.0d0)*tr ! Re(B0+C0) = -2Tr(Z^h U)

                            ! V_R = Z_L * TIJ_hLR
                            call ZGEMM('N', 'N', 3, pmax, pmax, Alpha, &  !use it later for D0R
                                       ZL_ireg(:, :, j), 3, TIJ_h(js + 1:js + pmax, ls + 1 + pmax:ls + lmax), pmax, &
                                       beta, VR(:, :, j, l), 3)

                        end if
                        if ((.NOT. RCP) .OR. (LP)) then
                            ! U_L = Z_L * T_LL
                            call ZGEMM('N', 'N', 3, pmax, pmax, Alpha, &
                                       ZL_ireg(:, :, j), 3, TIJ_h(js + 1:js + pmax, ls + 1:ls + pmax), pmax, &
                                       beta, UL(:, :, j, l), 3)
                            ! Z^h [Ujl]
                            call calcTrace('C', 'N', ZL_reg(:, :, l), UL(:, :, j, l), tr)

                            B0L = B0L + (2.0d0)*tr ! Re(B0+C0) = 2Tr(Z^h U)

                            ! V_L = Z_R * TIJ_hRL
                            call ZGEMM('N', 'N', 3, pmax, pmax, Alpha, &   !use it later for D0L
                                       ZR_ireg(:, :, j), 3, TIJ_h(js + pmax + 1:js + lmax, ls + 1:ls + pmax), pmax, &
                                       beta, VL(:, :, j, l), 3)
                        end if
                    end do
                end do
                !write(*,*)'B0R', B0R
                !write(*,*)'B0L', B0L

!---------------end calculation of B0+C0
!---------------Start writing D0:-----------------------------

                D0R = 0
                D0L = 0

                do j = 1, nscat
                    do l = 1, nscat
                        ls = (l - 1)*lmax
                        do i = 1, nscat
                            do k = 1, nscat
                                ks = (k - 1)*lmax

                                if (RCP .OR. LP) then

                                    t33 = VR(:, :, i, k) - UR(:, :, i, k)  !3*pmax
                                    t11 = VR(:, :, j, l) + UR(:, :, j, l)  !pmax*3

                                    ! [U^h + V^h] [V - U] = t11^h * t33
                                    call ZGEMM('C', 'N', pmax, pmax, 3, Alpha, &
                                               t11, 3, t33, 3, beta, t22, pmax)

                                    ! O_RR [U^h + V^h] [V - U] = VtacsKL_h22 * t22
                                    call calcTrace('N', 'N', vtacsKL_h(ks + pmax + 1:ks + lmax, ls + pmax + 1:ls + lmax), &
                                                   t22, tr)

                                    D0R = D0R + tr

                                end if

                                if ((.NOT. RCP) .OR. LP) then

                                    t33 = UL(:, :, i, k) - VL(:, :, i, k)  !3*pmax
                                    t11 = UL(:, :, j, l) + VL(:, :, j, l)  !pmax*3

                                    ! [U^h + V^h] [U - V] = t11^h * t33
                                    call ZGEMM('C', 'N', pmax, pmax, 3, Alpha, &
                                               t11, 3, t33, 3, beta, t22, pmax)

                                    ! O_LL [U^h + V^h] [U - V] = VtacsKL_h11 * t22
                                    call calcTrace('N', 'N', vtacsKL_h(ks + 1:ks + pmax, ls + 1:ls + pmax), &
                                                   t22, tr)

                                    D0L = D0L + tr

                                end if

                            end do
                        end do
                    end do
                end do
                !write(*,*)'D0R', D0R
                !write(*,*)'D0L', D0L
!---------------------------------------------------------

                if (RCP .OR. LP) then

                    Or_OC(n, 1) = -1.0d0 + 2.0d0*tpi*(realpart(B0R + D0R)) !Normalized value
                    Or_OC(n, 2) = 0.5d0*eps0*hostK_G*(-1.0d0 + 2.0d0*tpi*realpart(B0R + D0R)) !Total value
                    Or_OC(n, 3) = -1.0d0
                    Or_OC(n, 4) = 2.0d0*tpi*realpart(B0R)
                    Or_OC(n, 5) = 2.0d0*tpi*realpart(D0R)
                end if

                if (.not. RCP .AND. (.not. LP)) then
                    Or_OC(n, 1) = 1.0d0 + 2.0d0*tpi*(realpart(B0L + D0L)) !Normalized value
                    Or_OC(n, 2) = 0.5d0*eps0*hostK_G*(1.0d0 + 2.0d0*tpi*realpart(B0L + D0L)) !Total value
                    Or_OC(n, 3) = 1.0d0
                    Or_OC(n, 4) = 2.0d0*tpi*realpart(B0L)
                    Or_OC(n, 5) = 2.0d0*tpi*realpart(D0L)
                else
                    Or_OC(n, 6) = 1.0d0 + 2.0d0*tpi*(realpart(B0L + D0L)) !Normalized value
                    Or_OC(n, 7) = 0.5d0*eps0*hostK_G*(1.0d0 + 2.0d0*tpi*realpart(B0L + D0L)) !Total value
                    Or_OC(n, 8) = 1.0d0
                    Or_OC(n, 9) = 2.0d0*tpi*realpart(B0L)
                    Or_OC(n, 10) = 2.0d0*tpi*realpart(D0L)
                end if

            elseif (idepth == -1) then !inside sphere
                if (verb > 1) then
                    call cpu_time(t0)
                    call system_clock(tic)
                end if
                !do j=1,nscat
                j = p_label(n, 1)
                rj(1:3) = rp - geometry(1:3, j)

                call calcVSWs( &
                    r=rj, &
                    k=scatK_(j, abs(idepth)), &
                    pmax=lmax/2, &
                    regt=.true., &
                    cart=.true., &
                    waves=si_reg(:, :, j))

                ! enddo
                if (verb > 1) then   !Note:check timing!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    call cpu_time(t1)
                    call system_clock(toc, count_rate=tps)
                    t2 = t2 + t1 - t0
                    t3 = t3 + dble(toc - tic)/dble(tps)
                end if
                ! Precompute Lambda matrix (Stout 2002, eq:53)
                if (verb > 1) then
                    call cpu_time(t0)
                    call system_clock(tic)
                end if
                !do j=1,nscat
                js = (j - 1)*lmax
                call calcLamMat( &
                    Xi=cmplx(hostK_G*scatXYZR(4, j), 0, kind(hostK_G)), &
                    ro=scatK_(j, abs(idepth))/hostK_G, &
                    mat=lam_mat(js + 1:js + lmax, js + 1:js + lmax))

                !enddo
                if (verb > 1) then
                    call cpu_time(t1)
                    call system_clock(toc, count_rate=tps)
                    t4 = t4 + t1 - t0
                    t5 = t5 + dble(toc - tic)/dble(tps)
                end if
                if (RCP .OR. LP) then
                    do k = 1, nscat   !calculate VpUR =(U+V) & VmUR = (V-U)
                        ks = (k - 1)*lmax
                        !B = t22, A=t44

                        call ZGEMM('N', 'N', pmax, pmax, pmax, Alpha, &
                                lam_mat(js + 1:js + pmax, js + 1:js + pmax), pmax, &
                                (TIJ_h(js + 1:js + pmax, ks + 1 + pmax:ks + lmax) + &
                          & TIJ_h(js + pmax + 1:js + lmax, ks + pmax + 1:ks + lmax)), &
                           pmax, beta, t22, pmax)

                        ! t22 = matmul(lam_mat(js + 1:js + pmax, js + 1:js + pmax), &
                        !    TIJ_h(js + 1:js + pmax, ks + 1 + pmax:ks + lmax)+ &
                        !    & TIJ_h(js + pmax + 1:js + lmax, ks + pmax + 1:ks + lmax))

                        call ZGEMM('N', 'N', pmax, pmax, pmax, Alpha, &
                                lam_mat(js + 1 + pmax:js + lmax, js + 1 + pmax:js + lmax), pmax, &
                                (TIJ_h(js + 1:js + pmax, ks + 1 + pmax:ks + lmax) - &
                          & TIJ_h(js + pmax + 1:js + lmax, ks + pmax + 1:ks + lmax)), &
                          pmax, beta, t44, pmax)

                        ! t44 = matmul(lam_mat(js + 1 +pmax :js + lmax, js + 1 + pmax:js + lmax), &
                        !   TIJ_h(js + 1:js + pmax, ks + 1 + pmax:ks + lmax)- &
                        !   & TIJ_h(js + pmax + 1:js + lmax, ks + pmax + 1:ks + lmax))

                        VpUR(:, :, j, k) = sqrt(2.0d0)*(matmul(si_reg(:, 1:pmax, j), t22) + &
                                                        matmul(si_reg(:, pmax + 1:lmax, j), t44))

                        VmUR(:, :, j, k) = sqrt(2.0d0)*(matmul(si_reg(:, 1:pmax, j), t44) + &
                                                        matmul(si_reg(:, pmax + 1:lmax, j), t22))
                    end do
                end if
                if ((.NOT. RCP) .OR. LP) then
                    do k = 1, nscat
                        ks = (k - 1)*lmax
                        !B = t22, A=t44

                        call ZGEMM('N', 'N', pmax, pmax, pmax, Alpha, &
                                lam_mat(js + 1:js + pmax, js + 1:js + pmax), pmax, &
                                (TIJ_h(js + pmax + 1:js + lmax, ks + 1:ks + pmax) + &
                          & TIJ_h(js + 1:js + pmax, ks + 1:ks + pmax)), &
                           pmax, beta, t22, pmax)

                        ! t22 = matmul(lam_mat(js + 1:js + pmax, js + 1:js + pmax), &
                        !    TIJ_h(js + pmax +1:js + lmax, ks + 1 :ks + pmax)+ &
                        !    & TIJ_h(js + pmax + 1:js + lmax, ks + 1:ks + pmax))

                        call ZGEMM('N', 'N', pmax, pmax, pmax, Alpha, &
                                lam_mat(js + 1 + pmax:js + lmax, js + 1 + pmax:js + lmax), pmax, &
                                (TIJ_h(js + 1:js + pmax, ks + 1:ks + pmax) - &
                          & TIJ_h(js + pmax + 1:js + lmax, ks + 1:ks + pmax)), &
                          pmax, beta, t44, pmax)

                        ! t44 = matmul(lam_mat(js + 1 +pmax :js + lmax, js + 1 + pmax:js + lmax), &
                        !   TIJ_h(js + 1:js + pmax, ks + 1 :ks + pmax)- &
                        !   & TIJ_h(js + pmax + 1:js + lmax, ks + 1:ks + pmax))

                        VpUL(:, :, j, k) = sqrt(2.0d0)*(matmul(si_reg(:, 1:pmax, j), t22) + &
                                                        matmul(si_reg(:, pmax + 1:lmax, j), t44))

                        VmUL(:, :, j, k) = sqrt(2.0d0)*(matmul(si_reg(:, 1:pmax, j), t44) + &
                                                        matmul(si_reg(:, pmax + 1:lmax, j), t22))

                    end do

                end if

                B0R = 0
                B0L = 0
                do k = 1, nscat
                    ks = (k - 1)*lmax
                    do l = 1, nscat
                        ls = (l - 1)*lmax
                        if (RCP .OR. LP) then

                            call ZGEMM('C', 'N', pmax, pmax, 3, Alpha, &
                                       VpUR(:, :, j, l), 3, VmUR(:, :, j, k), 3, &
                                       beta, t22, pmax)

                            !t22 = matmul(dconjg(transpose(VpUR(:,:,j,l))), VmU(:,:,j,k))
                            call calcTrace('N', 'N', vtacsKL_h(ks + pmax + 1:ks + lmax, ls + pmax + 1:ls + lmax), &
                                           t22, tr)

                            B0R = B0R + tr
                        end if
                        if ((.NOT. RCP) .OR. LP) then
                            call ZGEMM('C', 'N', pmax, pmax, 3, Alpha, &
                                       VpUL(:, :, j, l), 3, VmUL(:, :, j, k), 3, &
                                       beta, t22, pmax)

                            !t22 = matmul(dconjg(transpose(VpUR(:,:,j,l))), VmU(:,:,j,k))
                            call calcTrace('N', 'N', vtacsKL_h(ks + 1:ks + pmax, ls + 1:ls + pmax), &
                                           t22, tr)

                            B0L = B0L + tr
                        end if
                    end do
                end do
                if (RCP .OR. LP) then
                    Or_OC(n, 1) = 2.0d0*tpi*scatK_(j, abs(idepth))*realpart(B0R)/hostK_G !Normalized value
                    Or_OC(n, 2) = tpi*eps0*scatK_(j, abs(idepth))*realpart(B0R) !Total value
                end if

                if (.not. RCP .AND. (.not. LP)) then
                    Or_OC(n, 1) = 2.0d0*tpi*scatK_(j, abs(idepth))*realpart(B0L)/hostK_G !Normalized value
                    Or_OC(n, 2) = tpi*eps0*scatK_(j, abs(idepth))*realpart(B0L) !Total value
                else
                    Or_OC(n, 6) = 2.0d0*tpi*scatK_(j, abs(idepth))*realpart(B0L)/hostK_G !Normalized value
                    Or_OC(n, 7) = tpi*eps0*scatK_(j, abs(idepth))*realpart(B0L) !Total value
                end if

!-------------------------------------------------

            end if
            if (verb > 1 .AND. (idepth .NE. 0)) write (*, '(i6,12x,A)') n, '-in '

        end do points
        ! report timings
        if (verb > 1) then
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [calcVSWs] (CPU & real in s): ', &
                t2, t3

            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [irreg calcVSWs] (CPU & real in s): ', &
                t4, t5

            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [reg calcVSWs] (CPU & real in s): ', &
                t6, t7

        end if

        write (*, *) ! empty line
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (verb > 0) write (*, '(A,A,/)') myname, '> Done!'

    end subroutine calcOaLDOC
!------------------------------------------------------------------------------
    subroutine calcOaExtField(r, geometry, TIJ, oaEB2, verb_, lambda, ehost, scatK_, p_label)

        !------------------------------------------------------------------
        ! This subroutine calculates the orientation average of the total external electric
        ! field intensity based on formulae 16 of Stout 2008
        !============================================================
        !note : by changing cart to true, I think we don't need prr any more!, but I write it according
        !paper
        use swav, only: calcVSWs, calcVTACs, xyz2rtp, calcVTrtp2xyz, calcLamMat
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables:
        real(8), intent(in), dimension(:, :) :: r
        !real(8), intent(in)  :: hostK
        integer, intent(in), optional :: verb_
        complex(8), intent(in) :: TIJ(:, :)
        complex(8), intent(in), optional ::  scatK_(:, :)
        real(8), intent(in) :: geometry(:, :), ehost, lambda
        real(8), intent(out), optional :: oaEB2(size(r, 2), 6)
        integer, intent(inout) :: p_label(:, :)
        ! Local variables:
        character(*), parameter :: myName = 'calcOaExtField'

        integer :: i, j, k, l, np, is, js, ks, ls, pmax, nmax, nscat, npts, lmax, n, idepth
        integer :: tic, toc, tps, verb
        real :: t0, t1, t2, t4
        real(8), parameter :: tpi = 2.0d0*pi
        real(8) :: A0, rp(3), t3, t5
        real(8) :: xyz(3), r2, rj(3), rnp(3), rdum, scatRot(3, 3, size(geometry, 2))
        complex(8) :: B0E, B0B
        complex(8), allocatable :: t11E(:, :), t11B(:, :), t22(:, :), t33(:, :), t44(:, :)
        real(8), allocatable ::  C0E, C0B
        complex(8), allocatable :: vtacsKL(:, :), si_reg(:, :, :), si_ireg(:, :, :)
        complex(8), allocatable :: psi_reg(:, :, :), psi_ireg(:, :, :)
        real(8), allocatable :: scatXYZR(:, :)
        complex(8) :: lam_mat(size(TIJ, 1), size(TIJ, 2)), tr
        real(8), parameter :: isqrt2 = 1/sqrt(2.0d0)
        complex(8), allocatable :: VE(:, :, :, :), VB(:, :, :, :)
        complex(8), parameter :: Alpha = cmplx(1.0d0, 0.0d0, kind(isqrt2)), &
                                 beta = cmplx(0.0d0, 0.0d0, kind(isqrt2))
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------

        nscat = size(geometry, 2)
        npts = size(r, 2)
        lam_mat = 0
        oaEB2 = 0
        if (size(TIJ, 1) /= size(TIJ, 2)) then
            write (*, '(A,A)') myName, '> ERROR: TIJ not square'
            STOP
        else
            lmax = size(TIJ, 1)/nscat
            pmax = lmax/2
            nmax = int(sqrt(dble(pmax)))
            if (2*nmax*(nmax + 2) /= lmax) then
                write (*, '(A,A)') myName, '> ERROR: 2*nmax*(nmax+2) /= lmax'
                STOP
            end if
        end if
        verb = 0
        if (present(verb_)) verb = verb_
!
        allocate (scatXYZR(4, nscat), vtacsKL(size(TIJ, 1), size(TIJ, 2)))
        allocate (si_reg(3, lmax, nscat), si_ireg(3, lmax, nscat))
        allocate (psi_reg(3, lmax, nscat), psi_ireg(3, lmax, nscat))
        scatXYZR(1:4, :) = geometry(1:4, :)
        vtacsKL = 0

        !compute translation
        if (verb > 1) then
            call cpu_time(t0)
            call system_clock(tic)
        end if
        do k = 1, nscat
            ks = (k - 1)*lmax
            do l = k + 1, nscat
                ls = (l - 1)*lmax
                call calcVTACs( &
                    r0=scatXYZR(1:3, k) - scatXYZR(1:3, l), &
                    k=cmplx(hostK_G, 0, kind(hostK_G)), &
                    regt=.true., &
                    vtacs=vtacsKL(ks + 1:ks + lmax, ls + 1:ls + lmax))

                vtacsKL(ls + 1:ls + lmax, ks + 1:ks + lmax) = &
                    conjg(transpose(vtacsKL(ks + 1:ks + lmax, ls + 1:ls + lmax)))
            end do
            vtacsKL(ks + 1:ks + lmax, ks + 1:ks + lmax) = 0
            do l = 1, lmax
                vtacsKL(ks + l, ks + l) = 1
            end do
        end do
        if (verb > 1) then
            call cpu_time(t1)
            call system_clock(toc, count_rate=tps)
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [calcVTACs] (CPU & real in s): ', &
                t1 - t0, dble(toc - tic)/dble(tps)
        end if
        !-------------------------------------

        allocate (t11E(lmax, lmax), t11B(lmax, lmax), t44(lmax, lmax), &
                  t22(lmax, lmax), t33(lmax, lmax), &
                  VE(3, lmax, nscat, nscat), VB(3, lmax, nscat, nscat))

        A0 = 1/tpi; 
        !if (verb > 1) then

        !   write(*,'(A,A)') myName, '> Calculation time [point in/out shape] (CPU & real in s): '
        !end
        !if (verb > 2) then
        !    write(*, '(A,A)') myName, '> ===================='
        !   write(*,'(A,A)') myName, '> Calculation time [point in/out shape] (CPU & real in s): oaE2 & oaB2*(omega)^2'
        ! endif
        t2 = 0.0d0 ! timing for calcVSWs
        t3 = 0.0d0 ! timing for calcVSWs
        points: do n = 1, npts

            rp = r(1:3, n)
            !-----------------------------check for points outside particles
            idepth = 0
            np = 0
            np_scat1: do np = 1, nscat

                r2 = geometry(4, np)*geometry(4, np)
                rnp(1:3) = rp - geometry(1:3, np)
                rdum = dot_product(rnp, rnp)
                ! if (rdum == 0.0d0) cycle points
                if (rdum < r2) then   ! Inside smallest circumscribing sphere

                    if (geometry(8, np) > tiny1) then ! check for spheroid
                        ! First compute the rotation matrices for spheroidal scatterers.
                        ! The matrices are used for checking if a point is inside a spheroid.
                        scatRot(:, :, np) = matmul(RotMatY(-geometry(6, np)), & ! \theta -> \beta
                                                   RotMatZ(-geometry(7, np)))  ! \phi -> \gamma
                        xyz(:) = matmul(scatRot(:, :, np), rnp(:))
                        xyz(1:2) = xyz(1:2)*geometry(8, np) ! scale x,y by aspect ratio
                        rdum = dot_product(xyz, xyz)
                        if (rdum < r2) then
                            idepth = 1
                            exit np_scat1
                        end if

                    else ! inside sphere
                        idepth = -1
                        p_label(n, 1) = np
                        exit np_scat1
                        ! shells: do k=1,abs(nint(geometry(8,j))) ! check shells
                        !       if(rdum < geometry(4+k,j)**2) then
                        !           p_label(n,2)=idepth
                        !          idepth = idepth - 1

                        !       else
                        !          exit shells
                        !       endif
                        !  enddo shells
                    end if
                    !  if(abs(idepth) > 0) exit np_scat1
                end if
            end do np_scat1

!------------------------------------------------------------
            if (idepth == 0) then
                !write(*,*)'rp, idepth', rp, idepth
                if (verb > 1) then
                    call cpu_time(t0)
                    call system_clock(tic)
                end if
                do j = 1, nscat
                    rj(1:3) = rp - geometry(1:3, j)
                    call calcVSWs( &
                        r=rj, &
                        k=cmplx(hostK_G, 0, kind(hostK_G)), &
                        pmax=lmax/2, &
                        regt=.false., &
                        cart=.true., &
                        waves=si_ireg(:, :, j), &
                        wavesB=psi_ireg(:, :, j))
                    call calcVSWs( &
                        r=rj, &
                        k=cmplx(hostK_G, 0, kind(hostK_G)), &
                        pmax=lmax/2, &
                        regt=.true., &
                        cart=.true., &
                        waves=si_reg(:, :, j), &
                        wavesB=psi_reg(:, :, j))

                end do
                if (verb > 1) then
                    call cpu_time(t1)
                    call system_clock(toc, count_rate=tps)
                    t2 = t2 + t1 - t0
                    t3 = t3 + dble(toc - tic)/dble(tps)
                end if
!----------------------------B0-----------------------
                B0E = 0
                B0B = 0
                do j = 1, nscat
                    js = (j - 1)*lmax
                    do l = 1, nscat
                        ls = (l - 1)*lmax
                        call ZGEMM('N', 'N', 3, lmax, lmax, Alpha, &
                                   si_ireg(:, :, j), 3, TIJ(js + 1:js + lmax, ls + 1:ls + lmax), lmax, &
                                   beta, VE(:, :, j, l), 3)

                        call calcTrace('C', 'N', si_reg(:, :, l), VE(:, :, j, l), tr)

                        B0E = B0E + tr
!
                        call ZGEMM('N', 'N', 3, lmax, lmax, Alpha, &
                                   psi_ireg(:, :, j), 3, TIJ(js + 1:js + lmax, ls + 1:ls + lmax), lmax, &
                                   beta, VB(:, :, j, l), 3)

                        call calcTrace('C', 'N', psi_reg(:, :, l), VB(:, :, j, l), tr)

                        B0B = B0B + tr
                    end do
                end do
                !------------------------------------end calculation of B0
!!!!!!!!!!!! Start writing C0:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                C0E = 0
                C0B = 0
                do j = 1, nscat
                    js = (j - 1)*lmax
                    do l = 1, nscat
                        ls = (l - 1)*lmax
                        do i = 1, nscat
                            is = (i - 1)*lmax
                            do k = 1, nscat
                                ks = (k - 1)*lmax

                                call ZGEMM('C', 'N', lmax, lmax, 3, Alpha, &
                                           VE(:, :, i, k), 3, VE(:, :, j, l), 3, &
                                           beta, t33, lmax)

                                call calcTrace('N', 'N', vtacsKL(ls + 1:ls + lmax, ks + 1:ks + lmax), &
                                               t33, tr)

                                C0E = C0E + realpart(tr)

                                !-----------------
                                call ZGEMM('C', 'N', lmax, lmax, 3, Alpha, &
                                           VB(:, :, i, k), 3, VB(:, :, j, l), 3, &
                                           beta, t33, lmax)

                                call calcTrace('N', 'N', vtacsKL(ls + 1:ls + lmax, ks + 1:ks + lmax), &
                                               t33, tr)

                                C0B = C0B + realpart(tr)

                            end do
                        end do
                    end do
                end do
!---------------------
                oaEB2(n, 1) = 1 + tpi*(2*realpart(B0E) + C0E)
                oaEB2(n, 2) = tpi*2*realpart(B0E)
                oaEB2(n, 3) = tpi*C0E
                oaEB2(n, 4) = (hostK_G**2 + tpi*(2*realpart(B0B) + C0B))
                oaEB2(n, 5) = tpi*2*realpart(B0B)
                oaEB2(n, 6) = tpi*C0B
                oaEB2(n, 4:6) = ((lambda/(tpi*sp_light))**2)*oaEB2(n, 4:6)
                !Note: (1/omega)*orAveB-->total orAveB2

            elseif (idepth == -1) then !inside sphere
                if (verb > 1) then
                    call cpu_time(t0)
                    call system_clock(tic)
                end if
                !do j=1,nscat
                j = p_label(n, 1)
                rj(1:3) = rp - geometry(1:3, j)

                call calcVSWs( &
                    r=rj, &
                    k=scatK_(j, abs(idepth)), &
                    pmax=lmax/2, &
                    regt=.true., &
                    cart=.true., &
                    waves=si_reg(:, :, j), &
                    wavesB=psi_reg(:, :, j))

                ! enddo
                if (verb > 1) then
                    call cpu_time(t1)
                    call system_clock(toc, count_rate=tps)
                    t2 = t2 + t1 - t0
                    t3 = t3 + dble(toc - tic)/dble(tps)
                end if
                ! Precompute Lambda matrix (Stout 2002, eq:53)
                if (verb > 1) then
                    call cpu_time(t0)
                    call system_clock(tic)
                end if
                !do j=1,nscat
                js = (j - 1)*lmax

                call calcLamMat( &
                    Xi=cmplx(hostK_G*scatXYZR(4, j), 0, kind(hostK_G)), &
                    ro=scatK_(j, abs(idepth))/hostK_G, &
                    mat=lam_mat(js + 1:js + lmax, js + 1:js + lmax))
                !enddo

                if (verb > 1) then
                    call cpu_time(t1)
                    call system_clock(toc, count_rate=tps)
                    t4 = t4 + t1 - t0
                    t5 = t5 + dble(toc - tic)/dble(tps)
                end if
!---------------------------------------------------
                B0E = 0
                B0B = 0
                do k = 1, nscat
                    ks = (k - 1)*lmax
                    do l = 1, nscat
                        ls = (l - 1)*lmax

                        call ZGEMM('N', 'N', lmax, lmax, lmax, Alpha, & !a common term between E^2 & B^2
                                   lam_mat(js + 1:js + lmax, js + 1:js + lmax), &
                                   lmax, TIJ(js + 1:js + lmax, ks + 1:ks + lmax), &
                                   lmax, beta, t22, lmax)

                        call ZGEMM('N', 'N', lmax, lmax, lmax, Alpha, & !a common term between E^2 & B^2
                                   lam_mat(js + 1:js + lmax, js + 1:js + lmax), &
                                   lmax, TIJ(js + 1:js + lmax, ls + 1:ls + lmax), &
                                   lmax, beta, t33, lmax)

                        call ZGEMM('N', 'N', lmax, lmax, lmax, Alpha, & !a common term between E^2 & B^2
                                   t33, lmax, vtacsKL(ls + 1:ls + lmax, ks + 1:ks + lmax), &
                                   lmax, beta, t44, lmax)

                        call ZGEMM('C', 'N', lmax, lmax, 3, Alpha, &  !for E^2
                                   si_reg(:, :, j), 3, si_reg(:, :, j), 3, &
                                   beta, t11E, lmax)

                        call ZGEMM('C', 'N', lmax, lmax, 3, Alpha, &  !for B^2
                                   psi_reg(:, :, j), 3, psi_reg(:, :, j), 3, &
                                   beta, t11B, lmax)

                        call ZGEMM('C', 'N', lmax, lmax, lmax, Alpha, & !for E^2
                                   t22, lmax, t11E, lmax, beta, t33, lmax)

                        call calcTrace('N', 'N', t33, t44, tr)
                        B0E = B0E + tr
                        !------------------------------
                        call ZGEMM('C', 'N', lmax, lmax, lmax, Alpha, & !for B^2
                                   t22, lmax, t11B, lmax, beta, t33, lmax)

                        call calcTrace('N', 'N', t33, t44, tr)
                        B0B = B0B + tr
                    end do
                end do

                oaEB2(n, 1) = tpi*realpart(B0E)
                oaEB2(n, 4) = (1/tpi)*(lambda/sp_light)**2*realpart(B0B)
            end if
            if (verb > 1 .AND. (idepth .NE. 0)) write (*, '(i6,12x,A)') n, '-in '

        end do points
        if (verb > 1) then
            write (*, '(A,A,2(1x,es10.3E2))') &
                myName, '> Calculation time [calcVSWs] (CPU & real in s): ', &
                t2, t3
        end if
        write (*, *)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (verb > 0) write (*, '(A,A,/)') myname, '> Done!'
    end subroutine calcOaExtField

!--------------------------------------------------------------------------
    subroutine calcTrace(TRANSA, TRANSB, A, B, tr)

        ! This subroutine calculate the trace of (A*B)
        !---------------------------------------------------
        ! Start of variable declarations.
        !---------------------------------------------------
        ! Passed variables:
        Character, intent(in) :: TRANSA, TRANSB
        complex(8), intent(in) :: A(:, :), B(:, :)
        complex(8), intent(out) :: tr
        ! Local variables:
        character(*), parameter :: myName = 'calcTrace'
        integer :: is, js
        !---------------------------------------------------
        ! End of variable declarations. Directives start now
        !---------------------------------------------------
        tr = 0.0
        if (TRANSA == 'C' .AND. TRANSB == 'N') then ! Tr(A^h B) = sum(a*_ji b_ji)
            do is = 1, size(A, 1)
                do js = 1, size(A, 2)
                    tr = tr + dconjg(A(is, js))*B(is, js)

                end do
            end do
        elseif (TRANSA == 'N' .AND. TRANSB == 'N') then ! Tr(AB) = sum(a_ij b_ji)
            do is = 1, size(A, 1)
                do js = 1, size(A, 2)
                    tr = tr + A(is, js)*B(js, is)

                end do
            end do
        elseif (TRANSA == 'N' .AND. TRANSB == 'C') then ! Tr(AB^h) = sum(a_ij b*_ij)
            do is = 1, size(A, 1)
                do js = 1, size(A, 2)
                    tr = tr + A(is, js)*dconjg(B(is, js))
                end do
            end do
        elseif (TRANSA == 'C' .AND. TRANSB == 'C') then ! Tr(A^h B^h) = sum(a*_ij b*_ji)
            do is = 1, size(A, 1)
                do js = 1, size(A, 2)
                    tr = tr + dconjg(A(is, js)*B(js, is))
                end do
            end do
        end if

    end subroutine calcTrace

end module multiscat
