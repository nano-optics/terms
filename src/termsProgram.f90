program termsProgram
   !
   use multiscat, only: mapNF, spectrumFF, calcStokesScaVec
   use HDFfive, only: h5_crtgrp, h5_wrt2file, h5_wrtvec2file
   !
   implicit none
   !
   character(*), parameter :: myname = 'termsProgram'
   real(8), parameter :: pi = acos(-1.0d0)
   real(8), parameter :: tpi = 2*pi
   !
   ! --- Parameters to be read from the input file: ----------------
   integer :: nscat = 0, scheme = 3, verb = 1, mode = 2
   integer :: ndipoles = 0

   integer :: ncut(3) = (/8, 8, -8/) ! ncut(3) := [n1, n2, tol]
   ! n1 is nmax [for irregular/initial offsetting when staging]
   ! n2 >= n1 [for regular offsetting when contracting]
   character(len=3), allocatable :: labels(:, :)
   character(len=18), allocatable :: string(:)
   real(8) :: ehost_const = 1
   real(8), allocatable :: geometry(:, :), dipoles(:, :)
   integer :: Epower = 2, pol_type = 1
   real(8), allocatable :: ehost(:)
   real(8), allocatable :: wavelen(:)
   real(8), allocatable :: incidences(:, :), sca_angles(:, :)
   real(8), allocatable :: inc(:, :)
   real(8) :: gridLB(3) = (/0.0d0, 0.0d0, 0.0d0/)
   real(8) :: gridUB(3) = (/0.0d0, 0.0d0, 0.0d0/)
   integer :: gridNbins(3) = (/0, 0, 0/), nGridPoints = 1
   logical :: noRTR = .false., HDF5_out = .false. !, HDF5_in= .true.
   logical :: nf = .false., PWinc = .false., Ef = .true.
   logical :: dump_E = .true., dump_B = .false., dump_C = .false., split_absOA = .false.
   ! --- Declare in/output arguments for the solve function
   integer, allocatable :: nselect(:, :, :, :)
   complex(8), allocatable :: escat(:, :, :)
   real(8), allocatable :: work2(:, :, :), work1(:, :, :, :), work3(:, :, :, :, :), work5(:, :)
   real(8), allocatable :: work6(:, :), work7(:, :)
   real(8), allocatable :: work4(:, :, :, :, :), N_OC(:, :, :, :, :)
   real(8), allocatable :: orAvextEB_int(:, :, :), oa_ldoc(:, :, :)
   ! --- Declare local variables
   logical :: MiePresent = .false., incfile = .false., Scafile = .false., ldum
   logical :: dump_oaE2 = .false., dump_oaB2 = .false., dump_oaLdoc = .false.
   integer :: tic, toc, tps, idum, jdum, kdum, mdum, ndum, odum, pdum
   real :: t1, t2
   real(8) :: rdum(3), Ave_N_OC
   real(8), allocatable :: StokesScaVec(:, :, :, :), StokesPhaseMat(:, :, :, :)
   real(8), allocatable :: diff_Sca_CS(:, :, :, :), jsig_abs_oa(:, :, :)
   integer, allocatable ::  p_label(:, :)
   character(len=32) :: str, str1
   character(len=64) :: filename, pfilename, filename1, filename2
   character(len=64) :: hdf5_filename, fname
   character(len=256), allocatable :: dfuns(:)
   character(len=256), allocatable :: tfiles(:)
   character(len=16), dimension(18), parameter :: filenames = (/ &
                                                  'csExtOA', 'csScaOA', 'csAbsOA', &
                                                  'cdExtOA', 'cdScaOA', 'cdAbsOA', &
                                                  'csExt1X', 'csSca1X', 'csAbs1X', &
                                                  'csExt2Y', 'csSca2Y', 'csAbs2Y', &
                                                  'csExt3R', 'csSca3R', 'csAbs3R', &
                                                  'csExt4L', 'csSca4L', 'csAbs4L'/)
   ! NOTE: gfortran complains when filenames are of different
   !       length, so added trailing blanks in last three.

   CHARACTER(LEN=32) ::  groupname, dsetname
   CHARACTER(LEN=256) :: attribute
   CHARACTER(LEN=32), allocatable :: SUBGRPs(:)!, SUB_SUBGRPs(:)
   ! ----------------------------------------------------------------
   call cpu_time(t1)
   call system_clock(tic)
   !
   ! Parse inputfile and prepare for calculation
   if (iargc() > 0) then
      call getarg(1, filename)
   else
      filename = 'inputfile'
   end if
   call readInputFile(inputfile=filename)
   !

   if (mode == 2) then ! spectrum for far-field properties
      !
      if (scheme > 0) then
         allocate (work2(6, 0:ncut(2), size(wavelen)))
         allocate (jsig_abs_oa(nscat + 1, 1, size(wavelen)))
         jsig_abs_oa = 0
      else
         allocate (work2(1, 1, 1)) ! not seeking sig_oa
         allocate (jsig_abs_oa(1, 1, 1))
         jsig_abs_oa = 0
      end if
      !
      if (size(incidences, 2) > 1) then
         ! 4th index spans incidences
         allocate (work1(3, 4, size(wavelen), size(incidences, 2)))
      else
         ! 4th index spans multipole order
         allocate (work1(3, 4, size(wavelen), ncut(2)))
      end if
      work1 = 0
      allocate (work3(4, nscat, 4, size(wavelen), size(incidences, 2)))
      work3 = 0
      do ndum = 1, size(incidences, 2)
         work1(1:3, 1, 1, ndum) = incidences(1:3, ndum)
         !work1(1,2,1,ndum) = pol_type
      end do
      if (size(incidences, 2) == 1) work1(1:3, 1, 1, ncut(2)) = work1(1:3, 1, 1, 1)

      call spectrumFF( &
         ncut=ncut, &
         wavelen=wavelen, &
         ehost=ehost, &
         geometry=geometry, &
         scheme=scheme, & !
         sig_oa_=work2, & ! main output
         sig_=work1, & !
         sig_abs_=work3, &
         tfiles_=tfiles, &
         escat_=escat, &
         nselect_=nselect, &
         noRTR_=noRTR, &
         verb_=verb, &
         jsig_abs_oa=jsig_abs_oa)
      !
      ! ------------------------------------------------------
      ! Deal with (multiple) single-orientation cross-sections
      ! ------------------------------------------------------

      if (.not. HDF5_out) then

         cs_loop: do tps = 1, 3
            jones_loop: do toc = 1, 4

               write (filename, '(A,A)') &
                  trim(filenames(6 + (toc - 1)*3 + tps)), '.dat'
               open (unit=11, file=filename, status='replace')
               !
               if (size(incidences, 2) > 1) then
                  write (11, '(A)', advance='NO') '# Weights:'
                  do ndum = 1, size(incidences, 2)
                     write (11, '(1x,f12.10)', advance='NO') incidences(4, ndum)
                  end do
                  write (11, *)
                  ! Print comment with column labels
                  write (11, '(A)') '# Wavelength || Average || inc1 | inc2 | ...'
               else
                  write (11, '(A)') '# Wavelength || Total || n=1 | n=2 | ...'
               end if
               !
               ! Now loop over wavelengths
               wave_loop: do kdum = 1, size(wavelen)
                  write (11, '(f12.6)', advance='NO') wavelen(kdum)

                  if (size(incidences, 2) > 1) then
                     ! Calculate the weighted averege over incidences
                     rdum(1) = dot_product(work1(tps, toc, kdum, :), incidences(4, :))
                     write (11, '(1x, ES24.17E2)', advance='NO') rdum(1)
                     do ndum = 1, size(incidences, 2)
                        write (11, '(1x, ES24.17E2)', advance='NO') &
                           work1(tps, toc, kdum, ndum)
                     end do
                     write (11, *) ! end of line
                  else
                     write (11, '(1x, ES24.17E2)', advance='NO') &
                        sum(work1(tps, toc, kdum, :))
                     do ndum = 1, ncut(1)
                        write (11, '(1x, ES24.17E2)', advance='NO') &
                           work1(tps, toc, kdum, ndum)
                     end do
                     write (11, *) ! end of line
                  end if
               end do wave_loop
               !
               close (11)
               !
               ! ------------------------------------------------------
               ! Deal with partial absorptions for Mie scatterers
               ! ------------------------------------------------------
               if (tps == 3) then ! absorption cross-section
                  do idum = 1, nscat ! loop over scatterers
                     if (geometry(8, idum) <= 0) then
                        kdum = abs(nint(geometry(8, idum)))
                        do jdum = 1, 1 + kdum ! loop over shells
                           write (filename, '(A,A,i0.3,A,i1,A)') &
                              filenames(6 + (toc - 1)*3 + tps) (1:7), '_scat', idum, &
                              'coat', kdum - jdum + 1, '.dat'
                           open (unit=11, file=filename, status='replace')
                           write (11, '(A)') &
                              '# Wavelength || Average || inc1 | inc2 | ...'
                           do ndum = 1, size(wavelen) ! loop over wavelengths
                              write (11, '(f12.6)', advance='NO') wavelen(ndum)
                              ! Calculate the weighted averege over incidences
                              rdum(1) = dot_product( &
                                        work3(jdum, idum, toc, ndum, :), &
                                        incidences(4, :))
                              write (11, '(1x, ES24.17E2)', advance='NO') rdum(1)
                              do mdum = 1, size(incidences, 2)
                                 write (11, '(1x, ES24.17E2)', advance='NO') &
                                    work3(jdum, idum, toc, ndum, mdum)
                              end do
                              write (11, *) ! end of line
                           end do
                           close (11)
                        end do
                     end if
                  end do
               end if
               !
            end do jones_loop
         end do cs_loop
         !
         ! -----------------------------------------------------
         ! Now deal with orientationally-averaged cross-sections
         ! -----------------------------------------------------
         !write(*,*)'scheme', scheme
         !write(*,*)'split_absOA', split_absOA
         if (scheme > 0) then
            if ((scheme /= 3) .and. split_absOA) then !@

               do ndum = 1, 4 !size(work2,1)

                  write (filename, '(A,A)') &
                     trim(filenames(ndum)), '.dat'
                  open (unit=11, file=filename, status='replace')
                  write (11, '(A)') '# Wavelength | Total '
                  do idum = 1, size(wavelen)
                     write (11, '(f12.6,1x,ES24.17E2)', advance='NO') &
                        wavelen(idum), work2(ndum, 0, idum)
                     write (11, *)
                  end do
                  close (11)
               end do
            else !@
               do ndum = 1, size(work2, 1)

                  write (filename, '(A,A)') &
                     trim(filenames(ndum)), '.dat'
                  open (unit=11, file=filename, status='replace')

                  ! open(unit=11,file=filenames(ndum),'.dat',status='replace')
                  write (11, '(A)') '# Wavelength | Total | n= 1 | ...'
                  do idum = 1, size(wavelen)
                     write (11, '(f12.6,1x,ES24.17E2)', advance='NO') &
                        wavelen(idum), work2(ndum, 0, idum)
                     do jdum = 1, ncut(1)
                        write (11, '(1x,ES24.17E2)', advance='NO') &
                           work2(ndum, jdum, idum)
                     end do
                     write (11, *)
                  end do
                  close (11)
               end do
            end if !@
            !
            if (split_absOA) then
               write (filename, '(A)') 'csAbsOA_split_particle.dat'
               open (unit=11, file=filename, status='replace')
               write (11, '(A)') &
                  '# Wavelength || Total || absOA(particle=1) | absOA(particle=2) | ...'

               do idum = 1, size(wavelen)
                  write (11, '(f12.6,1x,10000(ES24.17E2))') &
                     wavelen(idum), jsig_abs_oa(:, 1, idum)
               end do
               close (11)
            end if

         end if

      else

         fname = hdf5_filename
         groupname = 'Far-Field '

         allocate (SUBGRPs(3))
         SUBGRPs = (/'/Far-Field/fixed_incidence   ', '/Far-Field/partial_absorption', '/Far-Field/oa_incidence      '/)
         call h5_crtgrp(fname, groupname, SUBGRPs)
         groupname = '/Far-Field'
         dsetname = 'Wavelengths'
         attribute = 'Wavelength(nm)'
         call h5_wrtvec2file(fname, groupname, dsetname, wavelen, attribute)

         !--------------------------------------------------
         dsetname = 'Incidences'
         attribute = 'Euler Angles:(alpha, beta, gamma) and Weights'
         call h5_wrt2file(fname, groupname, dsetname, TRANSPOSE(incidences), attribute)
         if (size(incidences, 2) > 1) then
            !   write(11,'(A)') ' Average || inc1 | inc2 | ...'
            groupname = '/Far-Field/fixed_incidence   '
            ! dsetname='Weights'
            ! call h5_wrtvec2file(fname,groupname, dsetname,incidences(4,:),attr)

            if (allocated(work5)) deallocate (work5)
            allocate (work5(size(wavelen), size(incidences, 2) + 1))
         else
            if (allocated(work5)) deallocate (work5)
            allocate (work5(size(wavelen), ncut(1) + 1))
         end if

         cs_loop2: do tps = 1, 3
         jones_loop2: do toc = 1, 4
            !
            dsetname = trim(filenames(6 + (toc - 1)*3 + tps))
            wave_loop2: do kdum = 1, size(wavelen)

               if (size(incidences, 2) > 1) then
                  ! Calculate the weighted averege over incidences
                  rdum(1) = dot_product(work1(tps, toc, kdum, :), incidences(4, :))
                  work5(kdum, 1) = rdum(1)
                  work5(kdum, 2:size(incidences, 2) + 1) = work1(tps, toc, kdum, :)
                  attribute = 'Average, inc1, inc2, ...'
               else
                  work5(kdum, 1) = sum(work1(tps, toc, kdum, :))
                  work5(kdum, 2:ncut(1) + 1) = work1(tps, toc, kdum, :)
                  attribute = 'Total, n=1, n=2, ...'
               end if
            end do wave_loop2
            groupname = '/Far-Field/fixed_incidence   '

            call h5_wrt2file(fname, groupname, dsetname, work5, attribute)
            ! ------------------------------------------------------
            ! Deal with partial absorptions for Mie scatterers
            ! ------------------------------------------------------
            if (tps == 3) then ! absorption cross-section
               do idum = 1, nscat ! loop over scatterers
                  if (geometry(8, idum) <= 0) then
                     kdum = abs(nint(geometry(8, idum)))

                     do jdum = 1, 1 + kdum ! loop over shells
                        write (str, '(i0.3)') idum
                        write (str1, '(i1)') kdum - jdum + 1
                        dsetname = filenames(6 + (toc - 1)*3 + tps) (1:7)//'_scat'//trim(str)&
                             &//'coat'//trim(str1)
                        !  '# Wavelength || Average || inc1 | inc2 | ...'
                        attribute = 'Average, inc1, inc2, ...'
                        if (allocated(work6)) deallocate (work6)
                        allocate (work6(size(wavelen), size(incidences, 2) + 1))
                        do ndum = 1, size(wavelen) ! loop over wavelengths
                           ! Calculate the weighted averege over incidences
                           rdum(1) = dot_product( &
                                     work3(jdum, idum, toc, ndum, :), &
                                     incidences(4, :))
                           work6(ndum, 1) = rdum(1)
                           work6(ndum, 2:size(incidences, 2) + 1) = work3(jdum, idum, toc, ndum, :)
                        end do
                        groupname = '/Far-Field/partial_absorption'
                        call h5_wrt2file(fname, groupname, dsetname, work6, attribute)

                     end do

                  end if
               end do
            end if
            !
         end do jones_loop2
         end do cs_loop2

         ! -----------------------------------------------------
         ! Now deal with orientationally-averaged cross-sections
         ! -----------------------------------------------------
         !

         if (scheme > 0) then
            groupname = '/Far-Field/oa_incidence      '
            if (scheme /= 3 .and. split_absOA) then !@
               do ndum = 1, 4!size(work2,1)
                  dsetname = filenames(ndum)
                  if (allocated(work5)) deallocate (work5)
                  allocate (work5(size(wavelen), 1))

                  do idum = 1, size(wavelen)
                     work5(idum, 1) = work2(ndum, 0, idum) !Total
                  end do
                  attribute = 'Total'
                  call h5_wrt2file(fname, groupname, dsetname, work5, attribute)
               end do
            else !@
               do ndum = 1, size(work2, 1)
                  dsetname = filenames(ndum)
                  if (allocated(work5)) deallocate (work5)
                  allocate (work5(size(wavelen), ncut(1) + 1))

                  do idum = 1, size(wavelen)
                     work5(idum, 1) = work2(ndum, 0, idum) !Total
                     work5(idum, 2:ncut(1) + 1) = work2(ndum, 1:ncut(1), idum)

                  end do
                  attribute = 'Total, n= 1, n=2, ...'
                  call h5_wrt2file(fname, groupname, dsetname, work5, attribute)
               end do
            end if !@

            if (split_absOA) then

               groupname = '/Far-Field/oa_incidence      '
               dsetname = 'csAbsOA_split'
               !write(11,'(A)') &
               ! ' Wavelength || Total || absOA(particle=1) | absOA(particle=2) | ...'
               attribute = 'Total, absOA(particle=1), absOA(particle=2),  ...'
               if (allocated(work5)) deallocate (work5)
               allocate (work5(size(wavelen), nscat + 1))
               do idum = 1, size(wavelen)
                  work5(idum, :) = jsig_abs_oa(:, 1, idum)
               end do
               call h5_wrt2file(fname, groupname, dsetname, work5, attribute)
            end if

         end if

      end if !if HDF5_out
!------------------------------------------------------------------
   elseif (mode == 1) then ! map near field at different wavelength and Incidence, Optical chirality(OC)
      ! & and orientation average of OC & E2 & B2
      !write(*,*)'nGridPoints', nGridPoints
      nf = .true.
      if (ndipoles > 0) then
         kdum = ndipoles
      else
         PWinc = .true.
         if (.not. allocated(incidences)) then
            allocate (incidences(4, 1))
            incidences = 0.0d0
            incidences(4, 1) = 1.0d0
         end if
         kdum = size(incidences, 2)
      end if

      ! if(nGridPoints /= 0) then

      if (nGridPoints > 0) then ! calculate regular grid
         allocate (p_label(nGridPoints, 2))

         p_label = 0
         allocate (work3(3, 5, nGridPoints, size(wavelen), kdum), &
                   work4(3, 5, nGridPoints, size(wavelen), kdum))
         work3 = 0
         work4 = 0
         if (dump_C) allocate (N_OC(1, 1, nGridPoints, size(wavelen), kdum))
         if (scheme /= 0 .AND. (PWinc)) then  !these quantities are just calculated for PWinc??
            if (dump_oaE2 .OR. dump_oaB2) allocate (orAvextEB_int(nGridPoints, 6, size(wavelen)))
            if (dump_oaLdoc) allocate (oa_ldoc(nGridPoints, 10, size(wavelen)))
         end if

         call calcGridPoints(work3(1:3, 1, 1:nGridPoints, 1, 1))
         ! write(*,*)'work3', work3(1:3,1,1:nGridPoints,1,1)
         work4 = work3
      else ! read from file
         inquire (file=trim(pfilename), exist=ldum)
         if (.not. ldum) then
            write (*, '(A,A)') 'ERROR: Missing file ', pfilename
            STOP
         end if
         open (unit=11, file=trim(pfilename), status='old')
         read (11, *) nGridPoints
         allocate (p_label(nGridPoints, 2))

         p_label = 0
         allocate (work3(3, 5, nGridPoints, size(wavelen), kdum), &
                   work4(3, 5, nGridPoints, size(wavelen), kdum))
         work3 = 0
         work4 = 0
         if (dump_C) allocate (N_OC(1, 1, nGridPoints, size(wavelen), kdum))

         if (scheme /= 0 .AND. (PWinc)) then   !these quantities are just calculated for PWinc??
            if (dump_oaE2 .OR. dump_oaB2) allocate (orAvextEB_int(nGridPoints, 6, size(wavelen)))
            if (dump_oaLdoc) allocate (oa_ldoc(nGridPoints, 10, size(wavelen)))
         end if

         do ndum = 1, nGridpoints
            read (11, *) work3(1:3, 1, ndum, 1, 1)

         end do
         close (11)
         work4(1:3, 1, :, 1, 1) = work3(1:3, 1, :, 1, 1)

      end if
      ! else
      !  allocate(work3(3,5,1,size(wavelen),kdum),&
      !            work4(3,5,1,size(wavelen),kdum), p_label(1,2))
      !  work3 = 0
      !  work4 = 0
      ! p_label=0
      ! if (dump_C) allocate(N_OC(1,1,1,size(wavelen),kdum))
      ! if (scheme/=0  .AND. (PWinc)) then
      !    if (dump_oaE2 .OR. dump_oaB2) allocate(orAvextEB_int(1,6,size(wavelen)))
      !    if (dump_oaLdoc) allocate(oa_ldoc(1,10,size(wavelen)))
      ! endif
      !endif
      if (dump_C) N_OC = 0
      if (scheme /= 0 .AND. (PWinc)) then
         if (dump_oaE2 .OR. dump_oaB2) orAvextEB_int = 0
         if (dump_oaLdoc) oa_ldoc = 0
         !   SpAveInE=0
      end if

      if (PWinc) then
         allocate (inc(4, size(incidences, 2)))
         inc = 0
         !  rdum(1) = real(pol_type,kind(incidences))
         do pdum = 1, size(incidences, 2)
            rdum(1) = real(pol_type, kind(incidences))  !(1:3,pdum)
            inc(1, pdum) = rdum(1)
            inc(2:4, pdum) = incidences(1:3, pdum)
         end do
      else
         allocate (inc(4, ndipoles))
         inc = 0
         do pdum = 1, ndipoles
            inc(2:4, pdum) = dipoles(4:6, pdum)
         end do
      end if
      call mapNF( &
         ncut=ncut, &
         wavelen=wavelen, &
         inc=inc, &
         ehost=ehost, &
         geometry=geometry, &
         scheme=scheme, &
         field=work3, &
         Bfield=work4, &
         tfiles_=tfiles, &
         escat_=escat, &
         nselect_=nselect, &
         verb_=verb, &
         noRTR_=noRTR, &
         N_OC=N_OC, &
         orAvextEB_int=orAvextEB_int, &
         oa_ldoc=oa_ldoc, &
         dump_oaE2=dump_oaE2, &
         dump_oaB2=dump_oaB2, &
         p_label=p_label, &
         dipoles=dipoles)

      if (.not. HDF5_out) then
         if (nGridpoints > 0) then
            if (scheme /= 0) then
               if (dump_oaE2 .OR. dump_oaB2 .OR. dump_oaLdoc) then
                  write (filename, '(A)') 'mapOaQuantity.dat'

                  if (verb > 1) write (*, '(A,A,A)') &
                     myname, '> Dumping orientation average quantities &
                               &to file ', filename

                  open (12, file=filename, status='replace')
                  write (12, '(A, 1x)', advance='NO') 'lambda, x1, y1, z1,'
                  if (dump_oaE2) write (12, '(A, 1x)', advance='NO') 'OaE2'
                  if (dump_oaB2) write (12, '(A, 1x)', advance='NO') 'OaB2'
                  if (dump_oaLdoc) then
                     if (pol_type == 1 .OR. pol_type == 2) then
                        write (12, '(A, 1x)') 'oa_n_ldoc(RCP), oa_n_ldoc(LCP)'
                     elseif (pol_type == 3) then

                        write (12, '(A, 1x)') 'oa_n_ldoc(RCP)'
                     elseif (pol_type == 4) then
                        write (12, '(A, 1x)') 'oa_n_ldoc(LCP)'

                     end if
                  end if
                  write (12, *)
                  do odum = 1, size(wavelen)
                     do ndum = 1, nGridpoints
                        write (12, '(f7.2,1x, 3(ES24.17E2, 1x))', advance='NO') &
                           wavelen(odum), work3(1:3, 1, ndum, 1, 1)

                        if (dump_oaE2) write (12, '((ES24.17E2, 1x))', advance='NO') &
                           orAvextEB_int(ndum, 1, odum)
                        if (dump_oaB2) write (12, '((ES24.17E2, 1x))', advance='NO') &
                           orAvextEB_int(ndum, 4, odum)
                        if (dump_oaLdoc) then
                           if (pol_type == 1 .OR. pol_type == 2) then !linear x or y polarization
                              write (12, '(2(ES24.17E2, 1x))', advance='NO') &
                                 oa_ldoc(ndum, 1, odum), oa_ldoc(ndum, 6, odum)
                           else
                              write (12, '((ES24.17E2, 1x))', advance='NO') &
                                 oa_ldoc(ndum, 1, odum)
                           end if
                        end if
                        write (12, *)
                     end do
                  end do
                  close (12)
                  !********
               end if
            end if
            if (dump_C) then

               write (filename, '(A)') 'normalised_ldoc.dat'
               if (verb > 1) write (*, '(A,A,A)') &
                  myname, '> Dumping normalized optical chirality to file ', filename
               open (13, file=filename, status='replace')
               if (size(incidences, 2) == 1) then  !++++++++

                  write (13, *) 'lambda, x1, y1, z1, Inc1: n_ldoc1, ...'
                  do odum = 1, size(wavelen)
                     do ndum = 1, nGridpoints
                        write (13, '(f7.2,1x, 10000(1x,ES24.17E2))') &
                           wavelen(odum), work3(1:3, 1, ndum, 1, 1), &
                           (N_OC(1, 1, ndum, odum, 1))
                     end do
                  end do
               else  !+++++

                  write (13, *) 'lambda, x1, y1, z1, Ave(n_ldoc), Inc1: n_ldoc1, ...'
                  do odum = 1, size(wavelen)
                     do ndum = 1, nGridpoints
                        Ave_N_OC = dot_product( &
                                   N_OC(1, 1, ndum, odum, :), &
                                   incidences(4, :))
                        write (13, '(f7.2,1x, 10000(1x,ES24.17E2))') &
                           wavelen(odum), work3(1:3, 1, ndum, 1, 1), &
                           Ave_N_OC, &
                           (N_OC(1, 1, ndum, odum, pdum), pdum=1, size(incidences, 2))
                     end do
                  end do

               end if  !+++++
               close (13)
            end if
            if ((nf)) then

               if (dump_E) then
                  filename = 'map_E.dat'
                  if (verb > 1) write (*, '(A,A,A)') &
                     myname, '> Dumping E field map to file ', filename
                  call dumpNFs2TXTFile( &
                     filename, &
                     incidences, &
                     Epower, &
                     wavelen, &
                     work3, &
                     Ef, &
                     p_label)
               end if
               Ef = .false.
               if (dump_B) then
                  filename = 'map_B.dat'
                  if (verb > 1) write (*, '(A,A,A)') &
                     myname, '> Dumping B field map to file ', filename
                  call dumpNFs2TXTFile( &
                     filename, &
                     incidences, &
                     Epower, &
                     wavelen, &
                     work4, &
                     Ef, &
                     p_label)
               end if

            end if
         else
            if (verb > 1) write (*, '(A,A,i2)') myname, &
               '> Not dumping field map or OC since nGridpoints=', nGridPoints
         end if

      else !HDF5 if
         fname = hdf5_filename
         groupname = 'Near-Field'
         call h5_crtgrp(fname, groupname)
         dsetname = 'Wavelengths'
         attribute = 'Wavelength(nm)'
         call h5_wrtvec2file(fname, groupname, dsetname, wavelen, attribute)
         dsetname = 'Incidences'
         attribute = 'Euler Angles:(alpha, beta, gamma) and Weights'
         call h5_wrt2file(fname, groupname, dsetname, TRANSPOSE(incidences), attribute)
         dsetname = 'Gridpoints'
         attribute = 'x,y,z'
         call h5_wrt2file(fname, groupname, dsetname, TRANSPOSE(work3(1:3, 1, :, 1, 1)), attribute)

         if(nGridpoints > 0) then  !&&
            !dump
            if (scheme /= 0) then  !@
               if (dump_oaE2 .OR. dump_oaB2 .OR. dump_oaLdoc) then
                  filename = 'mapOaQuantity'
                  attribute = 'lambda, x, y, z, '
                  if (verb > 1) write (*, '(A,A,A)') &
                    myname, '> Dumping orientation average quantities &
                               &to file ', filename
                  if (allocated(work5)) deallocate (work5)
                  pdum = 0
                  if (dump_oaE2) then
                     pdum = 1
                     attribute = trim(attribute)//'OaE2 ,'
                  end if
                  if (dump_oaB2) then
                     pdum = pdum + 1
                     attribute = trim(attribute)//'OaB2 ,'
                  end if
                  if (dump_oaLdoc) then
                     pdum = pdum + 1
                     if (pol_type == 1 .OR. pol_type == 2) then
                        pdum = pdum + 1
                        attribute = trim(attribute)//'Oa_norm_ldoc(RCP), Oa_norm_ldoc(LCP)'
                     elseif (pol_type == 3) then
                        attribute = trim(attribute)//'Oa_norm_ldoc(RCP)'
                     elseif (pol_type == 4) then
                        attribute = trim(attribute)//'Oa_norm_ldoc(LCP)'
                     end if
                  end if

                  allocate (work5(nGridpoints*size(wavelen), 4 + pdum))

                  do odum = 1, size(wavelen)
                     do ndum = 1, nGridpoints
                        work5(ndum + (odum - 1)*nGridpoints, 1) = wavelen(odum)
                        work5(ndum + (odum - 1)*nGridpoints, 2:4) = work3(1:3, 1, ndum, 1, 1)
                        idum = 5
                        if (dump_oaE2) then
                           work5(ndum + (odum - 1)*nGridpoints, idum) = orAvextEB_int(ndum, 1, odum)
                           idum = idum + 1
                        end if
                        if (dump_oaB2) then
                           work5(ndum + (odum - 1)*nGridpoints, idum) = orAvextEB_int(ndum, 4, odum)
                           idum = idum + 1
                        end if
                        if (dump_oaLdoc) then
                           if (pol_type == 1 .OR. pol_type == 2) then !linear pol x or y
                              work5(ndum + (odum - 1)*nGridpoints, idum) = oa_ldoc(ndum, 1, odum)
                              work5(ndum + (odum - 1)*nGridpoints, idum + 1) = oa_ldoc(ndum, 6, odum)
                           else
                              work5(ndum + (odum - 1)*nGridpoints, idum) = oa_ldoc(ndum, 1, odum)
                           end if
                        end if
                     end do
                  end do

                  call h5_wrt2file(fname, groupname, filename, work5, attribute)
               end if
            end if  !@
            if ((dump_C)) then  !****

               filename = 'normalised_ldoc'
               if (verb > 1) write (*, '(A,A,A)') &
                  myname, '> Dumping normalized optical chirality to file ', filename

               if (size(incidences, 2) == 1) then  !++++++++
                  if (allocated(work5)) deallocate (work5)
                  allocate (work5(nGridpoints*size(wavelen), 5))

                  ! write(12,*)'lambda, x1, y1, z1, Inc1: ldoc1, ...'
                  ! write(13,*)'lambda, x1, y1, z1, Inc1: n_ldoc1, ...'
                  do odum = 1, size(wavelen)
                     do ndum = 1, nGridpoints
                        work5(ndum + (odum - 1)*nGridpoints, 1) = wavelen(odum)
                        work5(ndum + (odum - 1)*nGridpoints, 2:4) = work3(1:3, 1, ndum, 1, 1)
                        work5(ndum + (odum - 1)*nGridpoints, 5) = N_OC(1, 1, ndum, odum, 1)
                     end do
                  end do
                  attribute = 'lambda, x, y, z, n_ldoc'
                  call h5_wrt2file(fname, groupname, filename, work5(:, 1:5), attribute)
               else  !+++++
                  ! write(12,*)'lambda, x1, y1, z1, Ave(OC), Inc1: ldoc11, ...'
                  ! write(13,*)'lambda, x1, y1, z1, Ave(N_OC), Inc1: n_ldoc1, ...'
                  pdum = size(incidences, 2)
                  if (allocated(work5)) deallocate (work5)
                  allocate (work5(nGridpoints*size(wavelen), 4 + (1 + pdum)))
                  do odum = 1, size(wavelen)
                     do ndum = 1, nGridpoints

                        Ave_N_OC = dot_product( &
                                   N_OC(1, 1, ndum, odum, :), &
                                   incidences(4, :))
                        !
                        work5(ndum + (odum - 1)*nGridpoints, 1) = wavelen(odum)
                        work5(ndum + (odum - 1)*nGridpoints, 2:4) = work3(1:3, 1, ndum, 1, 1)
                        work5(ndum + (odum - 1)*nGridpoints, 5) = Ave_N_OC
                        work5(ndum + (odum - 1)*nGridpoints, 6:5 + pdum) = N_OC(1, 1, ndum, odum, :)
                     end do
                  end do
                  attribute = 'lambda, x1, y1, z1, Ave(n_ldoc), Inc1: n_ldoc1, ...'
                  call h5_wrt2file(fname, groupname, filename, work5, attribute)

               end if  !+++++

            end if !****
            if ((nf)) then  !#
               if (verb > 1 .AND. dump_E) then
                  write (*, '(A,A,A)') &
                     myname, '> Dumping E  field map to files ', filename
                  filename = 'map_E '
               elseif (verb > 1 .AND. dump_B) then
                  write (*, '(A,A,A)') &
                     myname, '> Dumping B field map to files ', filename
                  filename = 'map_B '
               end if
               if (dump_E) then
                  filename = 'map_E'
                  call dumpNFs2HDF5File( &
                     fname, &
                     groupname, &
                     filename, &
                     incidences, &
                     Epower, &
                     wavelen, &
                     work=work3, &
                     p_label=p_label)
               end if
               if (dump_B) then
                  filename = 'map_B'
                  call dumpNFs2HDF5File( &
                     fname, &
                     groupname, &
                     filename, &
                     incidences, &
                     Epower, &
                     wavelen, &
                     work=work4, &
                     p_label=p_label)
               end if

            end if  !#

         else
            if (verb > 1) write (*, '(A,A,i2)') myname, &
               '> Not dumping field map or OC since nGridpoints=', nGridPoints
         endif  !&&&

      end if

      if (allocated(work3)) deallocate (work3)
      if (allocated(work4)) deallocate (work4)

      !------------------------------------------------------------------
!Stokes Scattering Vector, Stokes phase matrix, and diff_scat cross section
! at different wavelength or different scattering angle
   elseif (mode == 3) then
      if (.not. allocated(Sca_angles)) then
         allocate (sca_angles(3, size(incidences, 2)))
         sca_angles(1:3, :) = incidences(1:3, :)
      end if

      allocate (inc(4, size(incidences, 2)))
      inc = 0
      do pdum = 1, size(incidences, 2)
         rdum(1) = real(pol_type, kind(incidences))  !(1:3,pdum)
         inc(1, pdum) = rdum(1)
         inc(2:4, pdum) = incidences(1:3, pdum)
      end do

      if (size(incidences, 2) > 1) then
         allocate (StokesScaVec(4, 1, size(incidences, 2), size(sca_angles, 2)))
         StokesScaVec = 0
         allocate (StokesPhaseMat(4, 4, size(incidences, 2), size(sca_angles, 2)))
         StokesPhaseMat = 0
         allocate (diff_Sca_CS(4, 4, size(incidences, 2), size(sca_angles, 2)))
         diff_Sca_CS = 0
      else
         allocate (StokesScaVec(4, 1, size(wavelen), size(sca_angles, 2)))
         StokesScaVec = 0
         allocate (StokesPhaseMat(4, 4, size(wavelen), size(sca_angles, 2)))
         StokesPhaseMat = 0
         allocate (diff_Sca_CS(4, 4, size(wavelen), size(sca_angles, 2)))
         diff_Sca_CS = 0
      end if

      call calcStokesScaVec( &
         sca_angles=sca_angles, &
         inc2=inc, &
         ncut=ncut, &
         wavelen=wavelen, &
         ehost=ehost, &
         geometry=geometry, &
         scheme=scheme, &
         tfiles_=tfiles, &
         escat_=escat, &
         nselect_=nselect, &
         noRTR_=noRTR, &
         verb_=verb, &
         StokesPhaseMat=StokesPhaseMat, &
         StokesScaVec=StokesScaVec, &
         diff_sca=diff_Sca_CS)

      if (.not. HDF5_out) then  !****
         write (filename, '(A)') 'Stokes_Sca_Vec.dat'
         write (filename1, '(A)') 'Stokes_phase_Mat.dat'
         write (filename2, '(A)') 'diff_Sca_CS.dat'
         if (verb > 1) then
            write (*, '(A,A,A,/)') &
               myname, '> Dumping Stokes Scattering Vector to file ', filename
            write (*, '(A,A,A,/)') &
               myname, '> Dumping Phase matrix to file', filename1
            write (*, '(A,A,A,/)') &
               myname, '> Dumping differential Scattering cross section to file', filename2
         end if
         open (11, file=filename, status='replace')
         open (12, file=filename1, status='replace')
         open (13, file=filename2, status='replace')
         if (size(incidences, 2) > 1) then
            write (11, *) 'Inc_angles(alpha, beta, gamma),Sca_angles(alpha, beta, gamma) &
                       &: i, j   I, Q,U,V '
            write (12, *) 'Inc_angles(alpha, beta, gamma),Sca_angles(alpha, beta, gamma) &
                       &: i, j   z11, z21, z31, z41, ...,z44'
            write (13, *) 'Inc_angles(alpha, beta, gamma),Sca_angles(alpha, beta, gamma) &
                       &: i, j diff_Sca_CS'
            do odum = 1, size(incidences, 2)
               do pdum = 1, size(sca_angles, 2)
                  write (11, '(3(1x,F8.6),1x, 3(1x, F8.6),1x, 7000(1x,ES24.17E2))') &
                     incidences(1:3, odum), sca_angles(1:3, pdum), StokesScaVec(1:4, 1, odum, pdum)

                  write (12, '(3(1x,F8.6),1x, 3(1x, F8.6),1x, 7000(1x,ES24.17E2))') &
                     incidences(1:3, odum), sca_angles(1:3, pdum), StokesPhaseMat(1:4, 1:4, odum, pdum)

                  write (13, '(3(1x,F8.6),1x, 3(1x, F8.6),1x, 1(1x,ES24.17E2))') &
                     incidences(1:3, odum), sca_angles(1:3, pdum), diff_Sca_CS(1, 1, odum, pdum)
               end do
            end do

         else

            write (11, *) 'Sca_angle(alpha, beta, gamma)  lambda(nm)     I             Q              U              V'
            write (12, *) 'Sca_angle(alpha, beta, gamma)  lambda(nm)     z11 z21 z31 z41, ...,z44'
            write (13, *) 'Sca_angle(alpha, beta, gamma)  lambda(nm)     diff_Sca_CS'
            do pdum = 1, size(sca_angles, 2)
               do odum = 1, size(wavelen)
                  write (11, '(3(1x, F8.6),1x,f7.2,1x, 7000(1x,ES24.17E2))') &
                     sca_angles(1:3, pdum), wavelen(odum), StokesScaVec(1:4, 1, odum, pdum)
                  write (12, '(3(1x, F8.6),1x,f7.2,1x, 7000(1x,ES24.17E2))') &
                     sca_angles(1:3, pdum), wavelen(odum), StokesPhaseMat(1:4, 1:4, odum, pdum)
                  write (13, '(3(1x, F8.6),1x,f7.2,1x, 1(1x,ES24.17E2))') &
                     sca_angles(1:3, pdum), wavelen(odum), diff_Sca_CS(1, 1, odum, pdum)
               end do
            end do
         end if
      else !****
         fname = hdf5_filename
         groupname = 'Polarimetry'
         call h5_crtgrp(fname, groupname)
         dsetname = 'Wavelengths'
         attribute = 'Wavelength(nm)'
         call h5_wrtvec2file(fname, groupname, dsetname, wavelen, attribute)
         dsetname = 'Sca_angles'
         attribute = 'Euler Angles:(alpha, beta, gamma)'
         call h5_wrt2file(fname, groupname, dsetname, TRANSPOSE(sca_angles), attribute)
         dsetname = 'Incidences'
         attribute = 'Euler Angles:(alpha, beta, gamma) and Weights'
         call h5_wrt2file(fname, groupname, dsetname, TRANSPOSE(incidences), attribute)

         write (filename, '(A)') 'Stokes_Sca_Vec'
         write (filename1, '(A)') 'Stokes_phase_Mat'
         write (filename2, '(A)') 'diff_Sca_CS'
         if (verb > 1) then
            write (*, '(A,A,A,/)') &
               myname, '> Dumping Stokes Scattering Vector to file ', filename
            write (*, '(A,A,A,/)') &
               myname, '> Dumping Phase matrix to file', filename1
            write (*, '(A,A,A,/)') &
               myname, '> Dumping differential Scattering cross section to file', filename2
         end if

         !open(11, file=filename, status='replace')
         !open(12, file=filename1, status='replace')
         if (size(incidences, 2) > 1) then
            ! write(11,*)'Inc(i)&Sca(j)     I             Q              U              V'
            ! write(12,*)'Inc(i)&Sca(j)    z11 z12 z13 z14, ...,z44
            if (allocated(work5)) deallocate (work5)
            allocate (work5(size(incidences, 2)*size(sca_angles, 2), 10))
            if (allocated(work6)) deallocate (work6)
            allocate (work6(size(incidences, 2)*size(sca_angles, 2), 22))
            if (allocated(work7)) deallocate (work7)
            allocate (work7(size(incidences, 2)*size(sca_angles, 2), 7))

            do odum = 1, size(incidences, 2)
               do pdum = 1, size(sca_angles, 2)
                  work5(pdum + (odum - 1)*size(sca_angles, 2), 1:3) = incidences(1:3, odum)
                  work5(pdum + (odum - 1)*size(sca_angles, 2), 4:6) = sca_angles(1:3, pdum)
                  work5(pdum + (odum - 1)*size(sca_angles, 2), 7:10) = StokesScaVec(1:4, 1, odum, pdum)
                  work6(pdum + (odum - 1)*size(sca_angles, 2), 1:3) = incidences(1:3, odum)
                  work6(pdum + (odum - 1)*size(sca_angles, 2), 4:6) = sca_angles(1:3, pdum)
                  do mdum = 1, 4
                     work6(pdum + (odum - 1)*size(sca_angles, 2), 7 + (mdum - 1)*4:6 + mdum*4) = &
                     & StokesPhaseMat(1:4, mdum, odum, pdum)
                  end do
                  work7(pdum + (odum - 1)*size(sca_angles, 2), 1:3) = incidences(1:3, odum)
                  work7(pdum + (odum - 1)*size(sca_angles, 2), 4:6) = sca_angles(1:3, pdum)
                  work7(pdum + (odum - 1)*size(sca_angles, 2), 7) = diff_Sca_CS(1, 1, odum, pdum)
               end do
            end do
            dsetname = 'Stokes_Sca_Vec'
            attribute = 'Inc_angles(alpha, beta, gamma),Sca_angles(alpha, beta, gamma), &
                         &I, Q, U, V '
            call h5_wrt2file(fname, groupname, dsetname, work5, attribute)
            dsetname = 'Stokes_phase_Mat'
            attribute = 'Inc_angles(alpha, beta, gamma),Sca_angles(alpha, beta, gamma), &
                         &z11, z21, z31, z41, ...,z44'
            call h5_wrt2file(fname, groupname, dsetname, work6, attribute)
            dsetname = 'diff_Sca_CS'
            attribute = 'Inc_angles(alpha, beta, gamma),Sca_angles(alpha, beta, gamma), &
                        &diff_Sca_CS'
            call h5_wrt2file(fname, groupname, dsetname, work7, attribute)
         else

            ! write(11,*)'Sca_angle(j)  lambda(nm)     I             Q              U              V'
            !write(12,*)'Sca_angle(j)  lambda(nm)     z11 z12 z13 z14, ...,z44'
            if (allocated(work5)) deallocate (work5)
            allocate (work5(size(wavelen)*size(sca_angles, 2), 8))
            if (allocated(work6)) deallocate (work6)
            allocate (work6(size(wavelen)*size(sca_angles, 2), 20))
            if (allocated(work7)) deallocate (work7)
            allocate (work7(size(wavelen)*size(sca_angles, 2), 5))
            work5 = 0
            work6 = 0
            work7 = 0
            !write(*,*)'size(wavelen)=', size(wavelen)
            ! write(*,*)'size(sca_angles,2)=', size(sca_angles,2)
            do pdum = 1, size(sca_angles, 2)
               do odum = 1, size(wavelen)
                  work5(odum + (pdum - 1)*size(wavelen), 1:3) = sca_angles(1:3, pdum)
                  work5(odum + (pdum - 1)*size(wavelen), 4) = wavelen(odum)
                  work5(odum + (pdum - 1)*size(wavelen), 5:8) = StokesScaVec(1:4, 1, odum, pdum)
                  work6(odum + (pdum - 1)*size(wavelen), 1:3) = sca_angles(1:3, pdum)
                  work6(odum + (pdum - 1)*size(wavelen), 4) = wavelen(odum)
                  do mdum = 1, 4
                     work6(odum + (pdum - 1)*size(wavelen), 5 + (mdum - 1)*4:4 + mdum*4) = &
                     & StokesPhaseMat(1:4, mdum, odum, pdum)
                  end do
                  work7(odum + (pdum - 1)*size(wavelen), 1:3) = sca_angles(1:3, pdum)
                  work7(odum + (pdum - 1)*size(wavelen), 4) = wavelen(odum)
                  work7(odum + (pdum - 1)*size(wavelen), 5) = diff_Sca_CS(1, 1, odum, pdum)
               end do
            end do
            dsetname = 'Stokes_Sca_Vec'
            attribute = 'Sca_angle(alpha, beta, gamma), lambda(nm), I, Q, U, V'
            call h5_wrt2file(fname, groupname, dsetname, work5, attribute)
            dsetname = 'Stokes_phase_Mat'
            attribute = 'Sca_angle(alpha, beta, gamma), lambda(nm),&
                      &z11, z21, z31, z41, ...,z44'
            call h5_wrt2file(fname, groupname, dsetname, work6, attribute)
            dsetname = 'diff_Sca_CS'
            attribute = 'Sca_angle(alpha, beta, gamma), lambda(nm), diff_Sca_CS'
            call h5_wrt2file(fname, groupname, dsetname, work7, attribute)
         end if

      end if !****

   else
      !
      write (*, '(A,A)') myName, '> WARNING: No calculations triggered!'

   end if  !diff mode
   !
   call cpu_time(t2)
   call system_clock(toc, count_rate=tps)
   write (*, '(A,A,2(1x,es10.3E2))') &
      myName, '> Program run time (CPU & real in s):', &
      t2 - t1, real(toc - tic)/real(tps)
   !
   STOP
   !
contains
   !
   !*****************************************************************
   !
   subroutine readInputFile(inputfile)
      !
      ! ==============================================================
      ! Read the inputfile with keywords and the corresponding
      ! parameter values, and the scatterer geometry. NOTE: to be
      ! compatible with PyTERMS, global arrays should _NOT_ be
      ! allocated here, but they currently are.
      ! ==============================================================
      !
      use multiscat, only: dumpTmatCol_G, dumpStagedA_G, dumpPrestagedA_G, &
                           DumpScaCoeff_G, DumpIncCoeff_G, &
                           isolve_G, balScale_G, balStout_G, transInv_G, tfilename_G, cs_stout_G
      !use optim, only : dumpGeometry_G
      !
      implicit none
      !
      character(*), intent(in) :: inputfile
      character(*), parameter :: myname = 'readInputFile'
      character(len=1024) :: sentence
      character(len=64) :: words(18)
      character(len=64) :: keyword
      logical :: yes
      integer, parameter :: u4inpFile = 22982
      integer :: eof, i, j, k, nkeywords, nwords, nhi, nlo
      real(8) :: x(8)
      integer, allocatable :: selections(:, :, :, :)
      real(8), allocatable :: alphas(:), betas(:), gammas(:), alphas_bar(:), betas_bar(:), gammas_bar(:)
      character(len=256), allocatable :: tfilenames(:)
      !
      inquire (file=inputfile, exist=yes)
      if (.not. yes) then
         write (*, '(A,A,A)') myname, '> ERROR: Missing file ', &
            trim(inputfile)
         STOP
      else
         write (*, '(A,A,A,/)') myname, '> Parsing file ', trim(inputfile)
      end if
      !
      ! Set some defaults
      allocate (wavelen(1)); wavelen(1) = 666
      allocate (incidences(4, 1)); incidences(1:3, 1) = 0; incidences(4, 1) = 1
      !allocate(Sca_angles(3,1)); Sca_angles(1:3,1) = 0; !Sca_angles(4,1) = 1
      !
      open (unit=u4inpFile, file=inputfile, status='old')
      !
      nkeywords = 0
      read (u4inpFile, '(A)', IOSTAT=eof) sentence
      !
      sentences: do while (eof == 0)
         !
         ! Split sentence into space-separated words
         call sentence2words(sentence, words)
         keyword = words(1)
         if (keyword /= '') nkeywords = nkeywords + 1
         !
         if (keyword == 'ModeAndScheme') then
            !
            write (*, '(A,A,A)') myName, '> Detected keyword ', &
               trim(keyword)
            !
            ! If present, must be the first keyword.
            !
            if (nkeywords /= 1) then
               write (*, '(A,A)') &
                  myname, '> ERROR: ''Mode'' not the first keyword!'
               STOP
            end if
            !
            if (words(2) /= '' .and. words(3) /= '') then
               read (words(2:3), *, iostat=eof) mode, scheme
               if (eof /= 0) call errorParsingArguments(keyword)
            else
               write (*, '(A,A,A)') &
                  myname, '> ERROR: Missing argument(s) for ', keyword
               STOP
            end if
            !

            if (mode == 2) then ! spectrumFF
               write (*, '(15x,A)') 'mode=2 => spectrum_FF for far-field quantities'
            elseif (mode == 1) then ! map field on grid at diff. lambda and Inc. field
               write (*, '(15x,A)') 'mode=1 => mapNF at diff. lambda and diff. Inc. and Sca_angles'
            elseif (mode == 3) then ! Stokes Vector at diff. lambda and Inc. field
               write (*, '(15x,A)') 'mode=3 => Stokes Phase Matrixes & Stokes Scattering Vectors'
            else
               write (*, '(A,A,i2)') myname, '> ERROR: Unrecognised mode=', mode
               STOP
            end if
            !
            if (scheme == 0) then
               write (*, '(15x,A)') &
                  'scheme=0 => Do not seek T^(ji) and just solve Ax=b for x'
            elseif (scheme == 1) then
               write (*, '(15x,A)') &
                  'scheme=1 => Seek T^(ji) by direct inversion of A in Ax=b'
            elseif (scheme == 2) then
               write (*, '(15x,A)') &
                  'scheme=2 => Seek T^(ji) using Stout''s iterative scheme'
            elseif (scheme == 3) then
               write (*, '(15x,A)') &
                  'scheme=3 => Seek T^(j) using Mackowski''s approach'
            else
               write (*, '(A,A,2(1x,i2))') &
                  myName, '> ERROR: Unrecognised scheme ', scheme
               STOP
            end if
            !
         else if (keyword == 'DisableRTR') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            noRTR = .true.
            !
         else if (keyword == 'ScattererCentredCrossSections') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            cs_stout_G = .true.
            split_absOA = .true.
         else if (keyword == 'DielectricFunctions') then
            !
            write (*, '(A,A,A)') myName, '> Detected keyword ', &
               trim(keyword)
            read (words(2), *, IOSTAT=eof) idum
            if (eof /= 0) call errorParsingArguments(keyword)
            !
            allocate (dfuns(idum))
            do i = 1, idum
               read (u4inpfile, *, IOSTAT=eof) dfuns(i)
               if (eof /= 0) then
                  write (*, '(A,A,i2)') myname, &
                     '> ERROR reading dielectric file for i= ', i
                  STOP
               end if
            end do
            !
         else if (keyword == 'TmatrixFiles') then
            !
            write (*, '(A,A,A)') myName, '> Detected keyword ', &
               trim(keyword)
            read (words(2), *, IOSTAT=eof) idum
            if (eof /= 0) call errorParsingArguments(keyword)
            !
            allocate (tfilenames(idum))
            do i = 1, idum
               read (u4inpfile, *, IOSTAT=eof) tfilenames(i)
               if (eof /= 0) then
                  write (*, '(A,A,i2)') myname, &
                     '> ERROR reading tfilename for i= ', i
                  STOP
               end if
            end do
            !
         else if (keyword == "MultipoleSelections") then
            !
            write (*, '(A,A,A)') myName, '> Detected keyword ', &
               trim(keyword)
            read (words(2), *, IOSTAT=eof) idum
            if (eof /= 0) call errorParsingArguments(keyword)
            !
            allocate (selections(2, 2, 2, idum))
            selections = 0
            !
            do i = 1, idum
               !
               read (u4inpFile, '(A)', IOSTAT=eof) sentence
               if (eof /= 0) then
                  write (*, '(A,A,i4,A,i6)') &
                     myname, '> ERROR with IOSTAT= ', eof, &
                     ' when reading multipole selection ', i
                  STOP
               end if
               !
               call sentence2words(sentence, words, nwords_=nwords)
               !>>> TEMP >>>
               ! if(nwords==2) then
               !    read(words(1:2),*,IOSTAT=eof) j,k
               !    selections(1,:,:,i)=j ! nlo
               !    selections(2,:,:,i)=k ! nhi
               ! elseif(nwords==8) then
               !    read(words(1:8),*,IOSTAT=eof) selections(:,:,:,i)
               ! else
               !    write(*,'(A,A,i1)') myname, &
               !         '> ERROR bad word count in  multipole selection ', i
               ! endif
               !<<< TEMP <<<

               if (nwords /= 2) then
                  write (*, '(A,A,i1)') myname, &
                     '> ERROR: nwords != 2 for multipole selection', i
               else
                  write (*, '(15x,A,i1,A,1x,A,1x,A)') &
                     'selection ', i, ':', trim(words(1)), trim(words(2))
               end if
               ! Parse the selection's descriptor string string
               yes = .true.
               do while (yes)
                  keyword = words(1) (1:2) ! specifies the block MM/EM/ME/EE
                  j = index(words(1), ':')
                  read (words(1) (3:j - 1), *, IOSTAT=eof) nlo
                  if (eof /= 0) then
                     write (*, '(A,A,i1)') myname, &
                        '> ERROR: failed reading nlo for selection', i
                     STOP
                  end if
                  k = index(words(1), '_')
                  if (k > 0) then
                     read (words(1) (j + 1:k - 1), *, IOSTAT=eof) nhi
                     words(1) (1:k) = ' ' ! blank up to and including "_"
                     words(1) = adjustl(words(1))
                  else
                     yes = .false.
                     read (words(1) (j + 1:), *, IOSTAT=eof) nhi
                  end if
                  if (eof /= 0) then
                     write (*, '(A,A,i1)') myname, &
                        '> ERROR: failed reading nlo for selection', i
                     STOP
                  end if
                  if (trim(keyword) == 'MM') then
                     idum = 1; jdum = 1
                  elseif (trim(keyword) == 'EM') then
                     idum = 2; jdum = 1
                  elseif (trim(keyword) == 'ME') then
                     idum = 1; jdum = 2
                  elseif (trim(keyword) == 'EE') then
                     idum = 2; jdum = 2
                  else
                     write (*, '(A,A,A)') myname, &
                        '> ERROR: Unrecognised block label ', trim(keyword)
                     STOP
                  end if
                  if (nlo == 0 .and. nhi == 0) nlo = 1 ! make nlo > nhi
                  selections(1, idum, jdum, i) = nlo
                  selections(2, idum, jdum, i) = nhi
               end do
               if (words(2) (1:1) == 'r') then
                  ! Apply the selection to rows only, keeping all the columns
                  selections(2, :, :, i) = -selections(2, :, :, i)
               elseif (words(2) (1:1) == 'c') then
                  ! Apply the selection to columns only, keeping all the rows
                  selections(1, :, :, i) = -selections(1, :, :, i)
               elseif (words(2) (1:1) /= 'b') then
                  write (*, '(A,A,A)') myname, &
                     '> ERROR: Unrecognised selection type: ', trim(words(2))
                  STOP
               end if
               !
               write (*, '(28x,A)') "T-block | n_lo | n_hi "
               write (*, '(29x,A,2(4x,i3))') "MM/11", selections(1:2, 1, 1, i)
               write (*, '(29x,A,2(4x,i3))') "EM/21", selections(1:2, 2, 1, i)
               write (*, '(29x,A,2(4x,i3))') "ME/12", selections(1:2, 1, 2, i)
               write (*, '(29x,A,2(4x,i3))') "EE/22", selections(1:2, 2, 2, i)
               !
            end do
            !
         else if (keyword == 'Verbosity') then ! verb level
            !
            write (*, '(A,A,A)') myName, '> Detected keyword ', &
               trim(keyword)
            !
            read (words(2), *, IOSTAT=eof) verb
            if (eof /= 0) call errorParsingArguments(keyword)
            if (verb < 0 .or. verb > 3) then
               write (*, '(A,A,i3)') myName, &
                  '> ERROR: Unrecognised verbosity= ', verb
               STOP
            else
               write (*, '(15x,A,i5)', advance='NO') 'verbosity= ', verb
            end if
            !
            if (verb == 0) then
               write (*, '(A)') ' (Silent)'
            elseif (verb == 1) then
               write (*, '(A)') ' (Low)'
            elseif (verb == 2) then
               write (*, '(A)') ' (Medium)'
            elseif (verb == 3) then
               write (*, '(A)') ' (High)'
            end if
            !
         else if (keyword == 'DumpCollectiveTmatrix') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            !
            dumpTmatCol_G = .true.
            !
            if (words(2) /= '') tfilename_G = words(2)
            !
         else if (keyword == 'DumpStagedA') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            !
            dumpStagedA_G = .true.
            !
         else if (keyword == 'DumpPrestagedA') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            !
            dumpPrestagedA_G = .true.
            !
         else if (keyword == 'DumpScaCoeff') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            !
            DumpScaCoeff_G = .true.
         else if (keyword == 'DumpIncCoeff') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            !
            DumpIncCoeff_G = .true.
         else if (keyword == 'Solver') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            !
            read (words(2), *, IOSTAT=eof) isolve_G
            if (eof /= 0) call errorParsingArguments(keyword)
            write (*, '(15x,A,i2)') 'Requested isolve=', isolve_G
            !
         else if (keyword == 'TransInverse') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            transInv_G = .true.
            !
         else if (keyword == 'DisableStoutBalancing') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            balStout_G = .false.
            !
         else if (keyword == 'BalanceScale') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            !
            read (words(2), *, IOSTAT=eof) balScale_G
            if (eof /= 0) call errorParsingArguments(keyword)
            write (*, '(A,A,ES24.17E2)') myname, &
               '> Balancing weights will be scaled by', balScale_G
            !
         else if (keyword == 'MultipoleCutoff') then
            ! maximal multipole order(s)
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            !
            read (words(2), *, IOSTAT=eof) ncut(1)
            if (eof /= 0) call errorParsingArguments(keyword)
            write (*, '(15x,A,i2)') 'Supplied ncut(1)= ', ncut(1)
            if (words(3) /= '') then
               read (words(3), *, IOSTAT=eof) ncut(2)
               if (ncut(2) < ncut(1)) then
                  write (*, '(A,A)') myName, '> ERROR: ncut(2) < ncut(1)'
                  STOP
               else
                  write (*, '(15x,A,i2)') 'Supplied ncut(2)= ', ncut(2)
               end if
            else
               ncut(2) = ncut(1)
               write (*, '(15x,A)') 'Setting ncut(2)= ncut(1)'
            end if
            if (words(4) /= '') then
               read (words(4), *, IOSTAT=eof) ncut(3)
               if (ncut(3) >= 0) then
                  write (*, '(A,A)') myName, '> ERROR: ncut(3) >= 0'
                  STOP
               else
                  write (*, '(15x,A,i3)') 'Supplied ncut(3)= ', ncut(3)
               end if
            else
               ncut(3) = -8
               write (*, '(15x,A,i3)') 'Setting ncut(3)= ', ncut(3)
            end if
            !
         else if (keyword == 'MapQuantity') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            !
            !  read(words(2:3),*,IOSTAT=eof) Epower   !Atefeh
            read (words(2), *, IOSTAT=eof) Epower
            if (eof /= 0) then    !call errorParsingArguments(keyword)
               if (words(2) /= '') then
                  dump_E = .false.
                  do idum = 2, 4
                     if (words(idum) (1:1) == 'E') dump_E = .true.
                     if (words(idum) (1:1) == 'B') dump_B = .true.
                     if (words(idum) (1:1) == 'C') dump_C = .true.
                  end do
               end if
            else
               if (Epower == 0) then
                  write (*, '(15x,A)') &
                     'Map field components Ex_re Ex_im Ey_re Ey_im Ez_re Ez_im'
               else
                  write (*, '(15x,A,i3)') &
                     'Map enhancement |E|**p with p= ', Epower
               end if
               if (words(3) /= '') then
                  dump_B = .false.
                  dump_E = .false.
                  dump_C = .false.
                  do idum = 3, 5
                     if (words(idum) (1:1) == 'E') dump_E = .true.
                     if (words(idum) (1:1) == 'B') dump_B = .true.
                     if (words(idum) (1:1) == 'C') dump_C = .true.
                  end do
               end if
            end if

         else if (keyword == 'MapOaQuantity') then
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            if (words(2) /= '') then
               do idum = 2, 4
                  if (words(idum) (1:1) == 'E') dump_oaE2 = .true.
                  if (words(idum) (1:1) == 'B') dump_oaB2 = .true.
                  if (words(idum) (1:1) == 'C') dump_oaLdoc = .true.
               end do
            else
               call errorParsingArguments(keyword)
            end if
!
         else if (keyword == 'SpacePoints') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            !
            if (abs(mode) /= 1) then
               write (*, '(A,A,i2)') myName, &
                  '> ERROR: SpacePoints incompatible with mode ', mode
               STOP
            end if
            !
            if (words(3) == '') then ! only filename supplied
               !
               nGridPoints = -1
               read (words(2), *, IOSTAT=eof) pfilename
               if (eof /= 0) call errorParsingArguments(keyword)
               !
            else ! spatial grid supplied
               !
               read (words(2:10), *, IOSTAT=eof) &
                  (gridLB(i), gridUB(i), gridNbins(i), i=1, 3)
               if (eof /= 0) call errorParsingArguments(keyword)
               do i = 1, 3
                  write (*, '(15x,A,i1,2(1x,ES24.17E2),i5)') &
                     'i,lb_i,ub_i,npts_i= ', &
                     i, gridLB(i), gridUB(i), gridNbins(i)
               end do
               !
               nGridPoints = &
                  (gridNbins(1) + 1)*(gridNbins(2) + 1)*(gridNbins(3) + 1)
               write (*, '(15x,A,i9)') 'nGridPoints= ', nGridPoints
            end if
            !
         else if (keyword == 'Wavelength') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)

            ! baptiste 20/04/2022: added option to pass a filename
            if (words(2) (1:1) == 'f' .or. words(2) (1:1) == 'F') then ! only filename supplied
               
               if (trim(words(3)) == '') then
                  write (*, '(A,A)') myname, &
                     '> ERROR: Missing wavelength filename'
                  STOP
               end if
               !
               write (*, '(15x,A,A)') 'Wavelength filename= ', trim(words(3))
               !
               inquire (file=trim(words(3)), exist=yes)
               if (.not. yes) then
                  write (*, '(A,A,A)') myName, &
                     '> ERROR: Missing wavelength file ', &
                     trim(words(3))
                  STOP
               end if
               !
               open (unit=11, file=trim(words(3)), status='old')
               read (11, *, IOSTAT=eof) k
               if (eof /= 0) then
                  write (*, '(A,A,A)') myName, '> ERROR: In header of ', trim(words(3))
                  STOP
               elseif (k < 1) then
                  write (*, '(A,A)') myName, '> ERROR: Wavelength count < 1'
                  STOP
               else
                  write (*, '(15x,A,i6)') 'Expected wavelength count= ', k
               end if
               if (allocated(wavelen)) deallocate (wavelen)
               allocate (wavelen(k))
               if (allocated(ehost)) deallocate (ehost)
               allocate (ehost(k))
               do j = 1, k
                  read (11, *, IOSTAT=eof) wavelen(j)
                  if (eof /= 0) then
                     write (*, '(A,A,i6)') myName, &
                        '> ERROR reading wavelength ', j
                     STOP
                  end if
               end do
               close (11)

            else if (words(2) /= '' .and. words(3) /= '' .and. words(4) /= '') then
               !
               !
               read (words(2:4), *, IOSTAT=eof) x(1), x(2), i
               if (eof /= 0) call errorParsingArguments(keyword)
               write (*, '(15x,A,f12.6)') 'Wavelength LB (nm): ', x(1)
               write (*, '(15x,A,f12.6)') 'Wavelength UB (nm): ', x(2)
               ! When reading from inputfile, x(2) is initially
               ! interpreted as the end-point of the wavelength
               ! range, so it needs to be converted into a step
               ! for discretising the range.
               if (i > 0) then
                  x(2) = (x(2) - x(1))/i
               else
                  i = 0
                  x(2) = 0.0d0
               end if
               write (*, '(15x,A,i9,f10.4)') 'nsteps, step: ', i, x(2)
               !
               if (allocated(wavelen)) deallocate (wavelen)
               allocate (wavelen(i + 1))
               if (allocated(ehost)) deallocate (ehost)
               allocate (ehost(i + 1))
               !
               do j = 0, i
                  wavelen(j + 1) = x(1) + j*x(2)
               end do
               !
            else if (words(2) /= '') then
               !
               if (allocated(wavelen)) deallocate (wavelen)
               allocate (wavelen(1))
               if (allocated(ehost)) deallocate (ehost)
               allocate (ehost(1))
               !
               read (words(2), *, IOSTAT=eof) wavelen(1)
               if (eof /= 0) call errorParsingArguments(keyword)
               write (*, '(15x,A,f10.4)') 'Wavelength (nm): ', wavelen(1)
               !
            else
               !
               write (*, '(15x,A,f10.4)') 'not happy', wavelen(1)
               call errorParsingArguments(keyword)
               !
            end if
            !
         else if (keyword == 'Medium') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            !
            read (words(2), *, IOSTAT=eof) ehost_const
            if (eof /= 0) call errorParsingArguments(keyword)
            if (ehost_const < 0) ehost = ehost**2
            write (*, '(15x,A,es10.4e1)') 'Constant host epsilon= ', ehost_const
            !
         else if (keyword == 'Incidence') then
            PWinc = .true.
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            !
            if (words(2) (1:1) == 'f' .or. words(2) (1:1) == 'F') then
               !

               if (trim(words(3)) == '') then
                  write (*, '(A,A)') myname, &
                     '> ERROR: Missing incidence filename'
                  STOP
               end if
               !
               write (*, '(15x,A,A)') 'Incidence filename= ', trim(words(3))
               incfile = .true.
               !
               inquire (file=trim(words(3)), exist=yes)
               if (.not. yes) then
                  write (*, '(A,A,A)') myName, &
                     '> ERROR: Missing incidence file ', &
                     trim(words(3))
                  STOP
               end if
               !
               open (unit=11, file=trim(words(3)), status='old')
               read (11, *, IOSTAT=eof) k
               if (eof /= 0) then
                  write (*, '(A,A,A)') myName, '> ERROR: In header of ', trim(words(3))
                  STOP
               elseif (k < 1) then
                  write (*, '(A,A)') myName, '> ERROR: Incidence count < 1'
                  STOP
               else
                  write (*, '(15x,A,i6)') 'Expected incidence count= ', k
               end if
               if (allocated(incidences)) deallocate (incidences)
               allocate (incidences(4, k))
               do j = 1, k
                  read (11, *, IOSTAT=eof) incidences(1:4, j)
                  if (eof /= 0) then
                     write (*, '(A,A,i6)') myName, &
                        '> ERROR reading incidence ', j
                     STOP
                  end if
               end do
               close (11)
               if (mode == 1 .OR. mode == 3) then
                  read (words(4), *, IOSTAT=eof) pol_type
                  if (eof /= 0) then
                     write (*, '(A,A)') &
                        myName, '> Polarization is X-linear by default '
                  elseif (pol_type == 1) then
                     write (*, '(A,A)') myName, '> X-linear polarization'
                  elseif (pol_type == -1) then
                     write (*, '(A,A)') myName, '> Y-linear polarization'
                     pol_type = 2 ! remap -1 -> 2
                  elseif (pol_type == 2) then
                     write (*, '(A,A)') myName, '> R-circular polarization'
                     pol_type = 3 ! remap 2 -> 3
                  elseif (pol_type == -2) then
                     write (*, '(A,A)') myName, '> L-circular polarization'
                     pol_type = 4 ! remap -2 -> 4
                  else
                     write (*, '(A,A)') myName, &
                        '> ERROR: Unrecogised polarization'
                     STOP
                  end if
               end if

            elseif (words(2) /= '' .and. words(3) /= '' .and. words(4) /= '') then
               !
               yes = .false.
               !
               ! Euler angles alpha (phi), beta (theta), and gamma
               x(1:3) = 0
               read (words(2:4), *, IOSTAT=eof) x(1:3)
               if (eof /= 0) then
                  write (*, '(A,A)') myName, &
                     '> ERROR reading Euler angles'
                  STOP
               end if
               !
               ! check if single incidence or grid...
               allocate (alphas(1)); alphas = 0
               if (x(1) < -1.5d0) then
                  k = -nint(x(1))
                  if (k > 1) yes = .true.
                  write (*, '(A,A,i4)') myName, &
                     '> Euler angle alpha [0,2Pi) grid-points: ', k
                  deallocate (alphas); allocate (alphas(k))
                  x(1) = tpi ! full range
                  if (mode == 1 .OR. mode == 3) then
                     if (words(6) /= '') then
                        read (words(6), *, IOSTAT=eof) j
                        !write(*,*) 'TEST!!!!! j=', j
                        if (eof == 0) then
                           if (j > 0) then
                              x(1) = x(1)/j
                              write (*, '(A,A,f12.8)') myName, &
                                 '> Modified alpha maximum: ', x(1)
                           end if
                        end if
                     end if
                  else
                     if (words(5) /= '') then
                        read (words(5), *, IOSTAT=eof) j
                        !write(*,*) 'TEST!!!!! j=', j
                        if (eof == 0) then
                           if (j > 0) then
                              x(1) = x(1)/j
                              write (*, '(A,A,f12.8)') myName, &
                                 '> Modified alpha maximum: ', x(1)
                           end if
                        end if
                     end if
                  end if
                  x(1) = x(1)/k ! step size in alpha (phi)
                  rdum(1) = 0.5d0*x(1)
                  do i = 0, k - 1
                     alphas(i + 1) = i*x(1) + rdum(1)
                  end do
               elseif (x(1) < 0 .or. x(1) >= tpi) then
                  write (*, '(A,A)') myName, &
                     '> ERROR: Badly valued Euler angle alpha'
                  STOP
               else
                  write (*, '(A,A,f12.8)') myName, &
                     '> Euler angle alpha [0,2Pi) value= ', x(1)
                  alphas(1) = x(1)
               end if
               !
               allocate (betas(1)); betas = 0
               if (x(2) < -1.5d0) then
                  k = -nint(x(2))
                  if (k > 1) yes = .true.
                  write (*, '(A,A,i4)') myName, &
                     '> Euler angle beta  [0,Pi]  grid-points: ', k
                  deallocate (betas); allocate (betas(k))
                  x(2) = pi ! full range
                  if (mode == 1 .OR. mode == 3) then
                     if (words(7) /= '') then
                        read (words(7), *, IOSTAT=eof) j
                        if (eof == 0) then
                           if (j > 0) then
                              x(2) = x(2)/j
                              write (*, '(A,A,f12.8)') myName, &
                                 '> Modified beta maximum: ', x(2)
                           end if
                        end if
                     end if
                  else
                     if (words(6) /= '') then
                        read (words(6), *, IOSTAT=eof) j
                        if (eof == 0) then
                           if (j > 0) then
                              x(2) = x(2)/j
                              write (*, '(A,A,f12.8)') myName, &
                                 '> Modified beta maximum: ', x(2)
                           end if
                        end if
                     end if
                  end if
                  rdum(2) = cos(x(2)) ! cos(beta_max)
                  ! Compute regular step size in cos(beta) [cos(theta)]
                  x(2) = (1 - rdum(2))/k
                  ! Find mid-point of lowest-lying increment
                  rdum(1) = rdum(2) + 0.5d0*x(2)
                  do i = 0, k - 1
                     betas(i + 1) = acos(rdum(1) + i*x(2))
                  end do
               elseif (x(2) < 0 .or. x(2) > pi) then
                  write (*, '(A,A)') myName, &
                     '> ERROR: Badly valued Euler angle beta'
                  STOP
               else
                  write (*, '(A,A,f12.8)') myName, &
                     '> Euler angle beta  [0,Pi]  value= ', x(2)
                  betas(1) = x(2)
               end if
               !
               allocate (gammas(1)); gammas = 0
               if (x(3) < -1.5d0) then
                  k = -nint(x(3))
                  if (k > 1) yes = .true.
                  write (*, '(A,A,i4)') myName, &
                     '> Euler angle gamma [0,2Pi) grid-points: ', k
                  deallocate (gammas); allocate (gammas(k))
                  x(3) = tpi ! full range
                  if (mode == 1 .OR. mode == 3) then
                     if (words(8) /= '') then
                        read (words(8), *, IOSTAT=eof) j
                        if (eof == 0) then
                           if (j > 0) then
                              x(3) = x(3)/j
                              write (*, '(A,A,f12.8)') myName, &
                                 '> Modified gamma maximum: ', x(3)
                           end if
                        end if
                     end if
                  else
                     if (words(7) /= '') then
                        read (words(7), *, IOSTAT=eof) j
                        if (eof == 0) then
                           if (j > 0) then
                              x(3) = x(3)/j
                              write (*, '(A,A,f12.8)') myName, &
                                 '> Modified gamma maximum: ', x(3)
                           end if
                        end if
                     end if
                  end if
                  x(3) = x(3)/k ! step size in gamma
                  rdum(1) = 0.5d0*x(3)
                  do i = 0, k - 1
                     gammas(i + 1) = i*x(3) + rdum(1)
                  end do
               elseif (x(3) < 0 .or. x(3) >= tpi) then
                  write (*, '(A,A)') myName, &
                     '> ERROR: Badly valued Euler angle gamma'
                  STOP
               else
                  write (*, '(A,A,f12.8)') myName, &
                     '> Euler angle gamma [0,2Pi) value= ', x(3)
                  gammas(1) = x(3)
               end if
               !
               if (abs(mode) == 1 .OR. abs(mode) == 3) then
                  !
                  read (words(5), *, IOSTAT=eof) pol_type
                  !
                  if (eof /= 0) then
                     write (*, '(A,A)') &
                        myName, '> ERROR: Failed to parse polarization'
                     STOP
                  elseif (pol_type == 1) then
                     write (*, '(A,A)') myName, '> X-linear polarization'
                  elseif (pol_type == -1) then
                     write (*, '(A,A)') myName, '> Y-linear polarization'
                     pol_type = 2 ! remap -1 -> 2
                  elseif (pol_type == 2) then
                     write (*, '(A,A)') myName, '> R-circular polarization'
                     pol_type = 3 ! remap 2 -> 3
                  elseif (pol_type == -2) then
                     write (*, '(A,A)') myName, '> L-circular polarization'
                     pol_type = 4 ! remap -2 -> 4
                  else
                     write (*, '(A,A)') myName, &
                        '> ERROR: Unrecogised polarization'
                     STOP
                  end if
               end if
               !
               ndum = size(alphas)*size(betas)*size(gammas)
               if (allocated(incidences)) deallocate (incidences)
               allocate (incidences(4, ndum))
               k = 0
               do idum = 1, size(alphas)
                  do jdum = 1, size(betas)
                     do kdum = 1, size(gammas)
                        k = k + 1
                        incidences(1, k) = alphas(idum)
                        incidences(2, k) = betas(jdum)
                        incidences(3, k) = gammas(kdum)
                        incidences(4, k) = real(1, kind(incidences))/ndum
                     end do
                  end do
               end do
               !
               deallocate (alphas, betas, gammas)
               !
            elseif (words(3) /= '' .or. words(4) /= '' .or. words(5) /= '') then
               write (*, '(A,A)') myName, &
                  '> ERROR: Some Euler angles missing'
               STOP
            end if
            !
            if (.not. allocated(incidences)) then
               write (*, '(A,A)') myName, &
                  '> ERROR: incidence array not allocated!'
               STOP
            end if
            !
            write (*, '(15x,A)') 'Incident Euler angles and weights:'
            write (*, '(10x,4(8x,A))') 'alpha', 'beta', 'gamma', 'weight'
            do k = 1, size(incidences, 2)
               write (*, '(12x,4(1x,f12.8))') incidences(1:4, k)
            end do

            !----------------------Atefeh add this section----------------
         else if (keyword == 'ScatteringAngles') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            !
            if (words(2) (1:1) == 'f' .or. words(2) (1:1) == 'F') then
               !
               if (mode /= 3) then
                  write (*, '(A,A)') myname, &
                     '> ERROR: scattering angles file incompatible with mode= ', mode
                  STOP
               elseif (trim(words(3)) == '') then
                  write (*, '(A,A)') myname, &
                     '> ERROR: Missing Scattering angles filename'
                  STOP
               end if
               !
               write (*, '(15x,A,A)') 'Scattering angles filename= ', trim(words(3))
               Scafile = .true.
               !
               inquire (file=trim(words(3)), exist=yes)
               if (.not. yes) then
                  write (*, '(A,A,A)') myName, &
                     '> ERROR: Missing Scattering angles file ', &
                     trim(words(3))
                  STOP
               end if
               !
               open (unit=11, file=trim(words(3)), status='old')
               read (11, *, IOSTAT=eof) k
               if (eof /= 0) then
                  write (*, '(A,A,A)') myName, '> ERROR: In header of ', trim(words(3))
                  STOP
               elseif (k < 1) then
                  write (*, '(A,A)') myName, '> ERROR: Scattering angles count < 1'
                  STOP
               else
                  write (*, '(15x,A,i6)') 'Expected Scattering angles count= ', k
               end if
               if (allocated(Sca_angles)) deallocate (Sca_angles)
               allocate (Sca_angles(3, k))
               do j = 1, k
                  read (11, *, IOSTAT=eof) Sca_angles(1:3, j)
                  if (eof /= 0) then
                     write (*, '(A,A,i6)') myName, &
                        '> ERROR reading Scattering angles ', j
                     STOP
                  end if
               end do
               close (11)
               !
            elseif (words(2) /= '' .and. words(3) /= '' .and. words(4) /= '') then
               !
               yes = .false.
               !
               ! Euler angles alpha (phi), beta (theta), and gamma
               x(1:3) = 0
               read (words(2:4), *, IOSTAT=eof) x(1:3)
               if (eof /= 0) then
                  write (*, '(A,A)') myName, &
                     '> ERROR reading Euler angles'
                  STOP
               end if
               !
               ! check if single Sca_angles or grid...
               allocate (alphas_bar(1)); alphas_bar = 0
               if (x(1) < -1.5d0) then
                  k = -nint(x(1))
                  if (k > 1) yes = .true.
                  write (*, '(A,A,i4)') myName, &
                     '> Euler angle alpha [0,2Pi) grid-points: ', k
                  deallocate (alphas_bar); allocate (alphas_bar(k))
                  x(1) = tpi ! full range
                  if (words(5) /= '') then
                     read (words(5), *, IOSTAT=eof) j
                     !write(*,*) 'TEST!!!!! j=', j
                     if (eof == 0) then
                        if (j > 0) then
                           x(1) = x(1)/j
                           write (*, '(A,A,f12.8)') myName, &
                              '> Modified alpha maximum: ', x(1)
                        end if
                     end if
                  end if
                  x(1) = x(1)/k ! step size in alpha (phi)
                  rdum(1) = 0.5d0*x(1)
                  do i = 0, k - 1
                     alphas_bar(i + 1) = i*x(1) + rdum(1)
                  end do
               elseif (x(1) < 0 .or. x(1) >= tpi) then
                  write (*, '(A,A)') myName, &
                     '> ERROR: Badly valued Euler angle alpha'
                  STOP
               else
                  write (*, '(A,A,f12.8)') myName, &
                     '> Euler angle alpha [0,2Pi) value= ', x(1)
                  alphas_bar(1) = x(1)
               end if
               !
               allocate (betas_bar(1)); betas_bar = 0
               if (x(2) < -1.5d0) then
                  k = -nint(x(2))
                  if (k > 1) yes = .true.
                  write (*, '(A,A,i4)') myName, &
                     '> Euler angle beta  [0,Pi]  grid-points: ', k
                  deallocate (betas_bar); allocate (betas_bar(k))
                  x(2) = pi ! full range
                  if (words(6) /= '') then
                     read (words(6), *, IOSTAT=eof) j
                     if (eof == 0) then
                        if (j > 0) then
                           x(2) = x(2)/j
                           write (*, '(A,A,f12.8)') myName, &
                              '> Modified beta maximum: ', x(2)
                        end if
                     end if
                  end if
                  !rdum(2) = x(2)
                  x(2) = x(2)/k ! step size in beta (phi)
                  ! Find mid-point of lowest-lying increment
                  rdum(1) = 0.5d0*x(2)
                  do i = 0, k - 1
                     betas_bar(i + 1) = rdum(1) + i*x(2)
                  end do
               elseif (x(2) < 0 .or. x(2) > pi) then
                  write (*, '(A,A)') myName, &
                     '> ERROR: Badly valued Euler angle beta'
                  STOP
               else
                  write (*, '(A,A,f12.8)') myName, &
                     '> Euler angle beta  [0,Pi]  value= ', x(2)
                  betas_bar(1) = x(2)
               end if
               !
               allocate (gammas_bar(1)); gammas_bar = 0
               if (x(3) < -1.5d0) then
                  k = -nint(x(3))
                  if (k > 1) yes = .true.
                  write (*, '(A,A,i4)') myName, &
                     '> Euler angle gamma [0,2Pi) grid-points: ', k
                  deallocate (gammas_bar); allocate (gammas_bar(k))
                  x(3) = tpi ! full range
                  if (words(7) /= '') then
                     read (words(7), *, IOSTAT=eof) j
                     if (eof == 0) then
                        if (j > 0) then
                           x(3) = x(3)/j
                           write (*, '(A,A,f12.8)') myName, &
                              '> Modified gamma maximum: ', x(3)
                        end if
                     end if
                  end if
                  x(3) = x(3)/k ! step size in gamma
                  rdum(1) = 0.5d0*x(3)
                  do i = 0, k - 1
                     gammas_bar(i + 1) = i*x(3) + rdum(1)
                  end do
               elseif (x(3) < 0 .or. x(3) >= tpi) then
                  write (*, '(A,A)') myName, &
                     '> ERROR: Badly valued Euler angle gamma'
                  STOP
               else
                  write (*, '(A,A,f12.8)') myName, &
                     '> Euler angle gamma [0,2Pi) value= ', x(3)
                  gammas_bar(1) = x(3)
               end if
               !

               !
               ndum = size(alphas_bar)*size(betas_bar)*size(gammas_bar)
               if (allocated(Sca_angles)) deallocate (Sca_angles)
               allocate (Sca_angles(3, ndum))
               k = 0
               do idum = 1, size(alphas_bar)
                  do jdum = 1, size(betas_bar)
                     do kdum = 1, size(gammas_bar)
                        k = k + 1
                        Sca_angles(1, k) = alphas_bar(idum)
                        Sca_angles(2, k) = betas_bar(jdum)
                        Sca_angles(3, k) = gammas_bar(kdum)
                        !Sca_angles(4,k) = real(1,kind(Sca_angles))/ndum
                     end do
                  end do
               end do
               !
               deallocate (alphas_bar, betas_bar, gammas_bar)
               !
            elseif (words(3) /= '' .or. words(4) /= '' .or. words(5) /= '') then
               write (*, '(A,A)') myName, &
                  '> ERROR: Some Euler angles missing'
               STOP
            end if
            if (mode == 3) then
               if (.not. allocated(Sca_angles)) then
                  write (*, '(A,A)') myName, &
                     '> ERROR: Scattering angles array not allocated!'
                  STOP
               end if
            end if
            !
            write (*, '(15x,A)') 'Scattering Euler angles:'
            write (*, '(10x,4(8x,A))') 'alpha', 'beta', 'gamma'!, 'weight'
            do k = 1, size(Sca_angles, 2)
               write (*, '(12x,4(1x,f12.8))') Sca_angles(1:3, k)
            end do
         else if (keyword == 'OutputFormat') then
            !
            write (*, '(A,A,A)') &
               myName, '> Detected keyword ', trim(keyword)
            if (words(2) == 'HDF5' .or. words(2) == 'TXT') then
               if (words(2) == 'HDF5') then
                  HDF5_out = .true.
                  write (*, '(15x,A,A)') &
                     'OutputFormat=', trim(words(2))
                  read (words(3), *, IOSTAT=eof) hdf5_filename
                  if (eof /= 0) then
                     write (*, '(15x,A)') &
                        'All output files are stored in file "results.h5"'
                     hdf5_filename = 'results.h5'

                  else
                     hdf5_filename = trim(words(3))//'.h5'
                     write (*, '(15x,3A)') &
                        'All output files are stored in file "', hdf5_filename, '"'

                  end if
               else
                  write (*, '(15x,A,A)') &
                     'OutputFormat=', trim(words(2))
                  write (*, '(15x,A)') &
                     'All output files are stored in Txt files'
               end if
            else
               write (*, '(/,A,A)') &
                  myname, '> ERROR: Wrong output Format!'
               STOP
            end if
            ! else if(keyword == 'differentialScatteringCS' ) then

            !  write(*,'(A,A,A)') &
            !   myName,'> Detected keyword ',trim(keyword)
            ! diffscattcs=.true.
            !--------------------------end section----------
            ! Dipole source
         else if (keyword == 'DipoleSource') then
            !
            write (*, '(/,A,A,A)', advance='NO') &
               myName, '> Detected keyword ', trim(keyword)

            read (words(2), *, IOSTAT=eof) ndipoles
            if (eof /= 0 .or. ndipoles < 1) then
               write (*, '(/,A,A)') &
                  myname, '> ERROR: Failed to read ndipoles!'
               STOP
            end if
            !
            write (*, '(1x,A,i4,/)') 'with ndipoles= ', ndipoles
            !
            allocate (dipoles(6, ndipoles))
            dipoleread: do i = 1, ndipoles
               read (u4inpFile, '(A)', IOSTAT=eof) sentence
               if (eof /= 0) then
                  write (*, '(A,A,i4,A,i6)') &
                     myname, '> ERROR with IOSTAT= ', eof, &
                     ' when reading line for dipole ', i
                  STOP
               end if
               !
               call sentence2words(sentence, words)
               read (words(1:6), *, IOSTAT=eof) dipoles(1:6, i)
               if (eof /= 0) then ! triggered by word='' ?!?
                  write (*, '(A,A,i4,A,i6)') &
                     myname, '> ERROR with IOSTAT= ', eof, &
                     ' when parsing (x,y,z,alpha, beta, gamma) for dipole ', i
                  write (*, *) words(1:6)
                  STOP
               end if
            end do dipoleread

            !-------------------------------------------------------
         else if (keyword == 'Scatterers') then
            !
            write (*, '(/,A,A,A)', advance='NO') &
               myName, '> Detected keyword ', trim(keyword)
            !
            ! Expected to be the last keyword!
            write (*, *)
            !
            read (words(2), *, IOSTAT=eof) nscat
            if (eof /= 0 .or. nscat < 1) then
               write (*, '(/,A,A)') &
                  myname, '> ERROR: Failed to read nscat!'
               STOP
            end if
            !
            write (*, '(1x,A,i4,/)') 'with nscat= ', nscat
            !
            allocate (geometry(8, nscat))
            allocate (labels(nscat, 4), string(nscat))
            allocate (escat(nscat, 4, size(wavelen)))
            allocate (tfiles(nscat))
            allocate (nselect(2, 2, 2, nscat))
            escat = 0
            geometry = 0
            labels = ''
            tfiles = '' ! empty string
            nselect = 0
            !
            scatread: do i = 1, nscat
               !
               read (u4inpFile, '(A)', IOSTAT=eof) sentence
               if (eof /= 0) then
                  write (*, '(A,A,i4,A,i6)') &
                     myname, '> ERROR with IOSTAT= ', eof, &
                     ' when reading line for scatterer ', i
                  STOP
               end if
               !
               call sentence2words(sentence, words)
               read (words(1), *, IOSTAT=eof) string(i)
               if (eof /= 0) then
                  write (*, '(A,A,i4,A,i6)') &
                     myname, '> ERROR with IOSTAT= ', eof, &
                     ' when parsing descriptor for scatterer ', i
                  STOP
               end if
               !
               !--------------------------------------------------------
               !read(words(2:9),*,IOSTAT=eof) geometry(1:8,i)
               read (words(2:5), *, IOSTAT=eof) geometry(1:4, i) !location x0, y0, z0, R0
               if (eof /= 0) then ! triggered by word='' ?!?
                  write (*, '(A,A,i4,A,i6)') &
                     myname, '> ERROR with IOSTAT= ', eof, &
                     ' when parsing (x,y,z,R) for scatterer ', i
                  write (*, *) words(2:5)
                  STOP
               end if
               !
               do j = 1, i - 1 ! Check for overlap
                  x(1) = norm2(geometry(1:3, i) - geometry(1:3, j))
                  if (x(1) < geometry(4, i) + geometry(4, j)) then
                     write (*, '(A,A,i4,A,i4,A)') &
                        myname, '> WARNING: circumscribing spheres ', &
                        i, ' and ', j, ' overlap!'
                  end if
               end do
               !
               ! Check for multipole selection
               j = index(string(i), '_S', back=.true.)
               if (j > 0) then
                  keyword = trim(string(i) (j + 2:))
                  string(i) (j:) = ' '
                  read (keyword, *) j
                  if (allocated(selections)) then
                     if (j <= size(selections, 4)) then
                        nselect(:, :, :, i) = selections(:, :, :, j)
                     else
                        write (*, '(A,A)') myname, &
                           '> ERROR: selection out of bounds!'
                        STOP
                     end if
                  else
                     write (*, '(A,A)') myname, '> ERROR: selections unallocated!'
                     STOP
                  end if
               end if
               !
               scatcheck: if (string(i) (1:2) == 'TF') then
                  !
                  ! T-matrix file
                  !
                  labels(i, 1) = trim(string(i))
                  read (labels(i, 1) (3:3), *) j
                  if (.not. allocated(tfilenames)) then
                     write (*, '(A,A)') myname, &
                        '> ERROR: Unallocated tfilenames'
                  elseif (j > size(tfilenames)) then
                     write (*, '(A,A)') myname, &
                        '> ERROR: T-matrix ID out of bounds'
                  else
                     tfiles(i) = trim(tfilenames(j))
                  end if
                  !
                  read (words(6:8), *, IOSTAT=eof) geometry(5:7, i)
                  if (eof /= 0) then
                     write (*, '(/,A,A,i4,A,i6)') &
                        myname, '> ERROR with IOSTAT= ', eof, &
                        ' when reading Euler angles for scatterer ', i
                     STOP
                  end if
                  !
                  if (words(9) /= '') then
                     read (words(9), *, IOSTAT=eof) geometry(8, i)
                     if (eof /= 0) then
                        write (*, '(/,A,A,i4,A,i6)') &
                           myname, '> ERROR with IOSTAT= ', eof, &
                           ' when reading asp.rat. of scatterer ', i
                        STOP
                     end if
                  else
                     geometry(8, i) = 1
                  end if
                  ! elseif(string(i)(1:2) == 'DP') then  !dipole source---------------------
                  !ndipoles=ndipoles+1
                  ! Dipole source polarizability file :Now we do not need this file
                  !
                  !labels(i,1) = trim(string(i))
                  !read(labels(i,1)(3:3), *) j
                  !if(.not. allocated(tfilenames)) then
                  !   write(*,'(A,A)') myname,&
                  !        '> ERROR: Unallocated tfilenames'
                  ! elseif(j > size(tfilenames)) then
                  !   write(*,'(A,A)') myname,&
                  !    '> ERROR: T-matrix ID out of bounds'
                  !else
                  !   tfiles(i) = trim(tfilenames(j))
                  !endif
                  !
                  !read(words(6:8),*,IOSTAT=eof) geometry(5:7,i)
                  !if(eof /= 0) then
                  !  write(*, '(/,A,A,i4,A,i6)') &
                  !        myname,'> ERROR with IOSTAT= ',eof, &
                  !        ' when reading Euler angles for scatterer ',i
                  !   STOP
                  !endif
                  !
                  ! if(words(9) /= '') then
                  !   read(words(9),*,IOSTAT=eof) geometry(8,i)
                  !   if(eof /= 0) then
                  !     write(*, '(/,A,A,i4,A,i6)') &
                  !         myname,'> ERROR with IOSTAT= ',eof, &
                  !       ' when reading asp.rat. of scatterer ',i
                  !   STOP
                  ! endif
                  ! else
                  !    geometry(8,i) = 1
                  ! endif

                  !-----------------------------------------------------------------------
               else ! Mie scatterer
                  !
                  MiePresent = .true.
                  k = 0 ! shell count
                  keyword = string(i)
                  ! read string from back!
                  j = index(keyword, '@', back=.true.)
                  labels(i, k + 1) = trim(keyword(j + 1:))
                  do while (j /= 0)
                     keyword(j:) = ''
                     k = k + 1
                     j = index(keyword, '@', back=.true.)
                     labels(i, k + 1) = trim(keyword(j + 1:))
                  end do
                  if (k > 3) then
                     write (*, '(/,A,A,i6)') myName, &
                        '> ERROR: Too many shells in scatterer ', i
                     STOP
                  else
                     !
                     read (words(5 + 1:5 + k), *, IOSTAT=eof) &
                        geometry(4 + 1:4 + k, i)
                     !
                     if (k > 0) geometry(8, i) = -k
                     !
                  end if
                  !
               end if scatcheck
            end do scatread
            !
            ! Now just print a bunch of stuff
            !
            write (*, '(A,A,/)') myname, &
               '> Descriptor(s) and circumscribing sphere(s):'
            write (*, '(1x,A,1x,A,18x,A,10x,A,10x,A,9x,A)') &
               'scatID', 'String', 'x', 'y', 'z', 'R_0'
            do i = 1, nscat
               write (*, '(i5,3x,A,4(1x,es10.4E1))') &
                  i, adjustl(string(i)), geometry(1:4, i)
            end do
            !
            write (*, '(/,A,A,/)') myname, &
               '> Individual geometry characteristic(s):'
            write (*, '(1x,A,1x,A)') &
               'scatID', 'Details'
            do i = 1, nscat
               if (geometry(8, i) > 1.0d-8) then
                  write (*, '(i5,3x,A,3(1x,f10.7),1x,A,1x,f10.7)') &
                     i, 'Euler angles:', geometry(5:7, i), &
                     'asp. ratio:', geometry(8, i)
               else
                  write (*, '(i5,3x,A,1x,i1)', advance='NO') &
                     i, 'Mie with ncoats=', nint(-geometry(8, i))
                  if (nint(-geometry(8, i)) > 0) then
                     write (*, '(A)', advance='NO') ' R_{-k}: '
                     do j = 1, nint(-geometry(8, i))
                        write (*, '(1x,es9.4E1)', advance='NO') &
                           geometry(4 + j, i)
                     end do
                  end if
                  write (*, *)
               end if
            end do
            !
            ! Print material parameters:
            if (MiePresent) then
               write (*, '(/,A,A,/)') myname, &
                  '> Dielectric functions for (coated) Mie scatterer(s):'
               write (*, '(1x,A,1x,A,1x,A)') &
                  'scatID', 'volID', 'Label'
               do i = 1, nscat
                  do k = 1, nint(-geometry(8, i)) + 1
                     write (*, '(i5,1x,i5,5x,A)') &
                        i, 1 - k, trim(labels(i, k))
                  end do
               end do
            end if
            !
            exit sentences
            !
         elseif (keyword /= '') then
            !
            write (*, '(A,A,A)') &
               myname, '> ERROR: Bad keyword ', keyword
            STOP
            !
         end if
         !
         read (u4inpFile, '(A)', IOSTAT=eof) sentence
         !
      end do sentences

      close (u4inpFile)
      if (mode == 3) then
         if (.not. allocated(Sca_angles)) then
            write (*, '(A,A)') myName, &
               '> WARNING: Scattering angles array not allocated!', &
               '>Scattering angles are the same as incidence angles'
         end if

         if (size(incidences, 2) > 1 .AND. size(wavelen) > 1) then
            write (*, '(A,A)') myName, &
               '> ERROR: Both incidence and wavelength cannot vary at mode 3 '
            STOP
         end if

      end if

      write (*, '(/,A,A,i2,A,/)') myname, '> Finished parsing ', &
         nkeywords, ' keywords'
      !
      call calcEpsilon()

      ehost(:) = ehost_const
      !
      if (allocated(tfilenames)) deallocate (tfilenames)
      if (allocated(selections)) deallocate (selections)
      !if (ndipoles > 0 .AND. PWinc ) then
      !   write(*,'(A,A)') myName, &
      !    '> ERROR: Both plane wave incidence and dipole source are not allowed'
      !  STOP
      ! elseif (ndipoles == 0 .AND. (.not. PWinc) ) then
      !   write(*,'(A,A)') myName, &
      !    '> ERROR: Please determine the source (Plane Wave or dipole)'
      !  STOP
      !endif
   end subroutine readInputFile
   !
   !*****************************************************************g
   !
   subroutine errorParsingArguments(keyword)
      !
      character(*), intent(in) :: keyword
      character(*), parameter :: myname = 'errorParsingArguments'
      !
      write (*, '(A,A,A)') myname, '> ', trim(keyword)
      STOP
      !
   end subroutine errorParsingArguments
   !
   !*****************************************************************
   !
   subroutine calcEpsilon()
      !
      ! ==============================================================
      ! Update the escat to a given wavelength.
      ! ==============================================================
      !
      use eps, only: epsAu, epsAg, epsAuRaschke, epsAgRaschke, epsPd, epsPt, epsSi, epsAl, epsCr, epsWater, interp1
      !
      implicit none
      !
      character(*), parameter :: myname = 'calcEpsilon'
      integer :: i, j, k, l, io
      real(8), allocatable :: lambdas(:)
      complex(8), allocatable :: epsilons(:)
      !
      shells: do k = 1, 4
         scatterers: do i = 1, nscat
            if (trim(labels(i, k)) == '' .or. &
                labels(i, k) (1:2) == 'TF') then ! no update
               cycle scatterers
            elseif (trim(labels(i, k)) == 'Au') then
               do j = 1, size(wavelen)
                  escat(i, k, j) = epsAu(wavelen(j))
               end do
            elseif (trim(labels(i, k)) == 'AuRaschke') then
               do j = 1, size(wavelen)
                  escat(i, k, j) = epsAuRaschke(wavelen(j))
            end do
            elseif (trim(labels(i, k)) == 'AgRaschke') then
               do j = 1, size(wavelen)
                  escat(i, k, j) = epsAgRaschke(wavelen(j))
            end do
            elseif (trim(labels(i, k)) == 'Ag') then
               do j = 1, size(wavelen)
                  escat(i, k, j) = epsAg(wavelen(j))
               end do
            elseif (trim(labels(i, k)) == 'Pd') then
               do j = 1, size(wavelen)
                  escat(i, k, j) = epsPd(wavelen(j))
               end do
            elseif (trim(labels(i, k)) == 'Pt') then
               do j = 1, size(wavelen)
                  escat(i, k, j) = epsPt(wavelen(j))
               end do
            elseif (trim(labels(i, k)) == 'Si') then

               escat(i, k, 1:size(wavelen)) = epsSi(wavelen)

            elseif (trim(labels(i, k)) == 'Al') then

               escat(i, k, 1:size(wavelen)) = epsAl(wavelen)

            elseif (trim(labels(i, k)) == 'Cr') then

               escat(i, k, 1:size(wavelen)) = epsCr(wavelen)
            elseif (trim(labels(i, k)) == 'Water') then

               escat(i, k, 1:size(wavelen)) = epsWater(wavelen)

            elseif (labels(i, k) (1:2) == 'DF') then
               !
               read (labels(i, k) (3:), *) j
               if (.not. allocated(dfuns)) then
                  write (*, '(A,A)') myname, &
                     '> ERROR: dfuns array no allocated!'
                  STOP
               end if
               ! check if dfuns(j) contains a value or filename
               read (dfuns(j), *, IOSTAT=io) rdum(1), rdum(2)

               if (io == 0) then
                  escat(i, k, :) = cmplx(rdum(1), rdum(2), kind(escat))

               else
                  l = countLines(trim(dfuns(j)))
                  if (l < 1) then
                     write (*, '(A,A,A)') myname, &
                        '> ERROR: Missing/empty dielectric file ', &
                        trim(dfuns(j))
                     STOP
                  end if
                  allocate (lambdas(l), epsilons(l))
                  open (unit=11, file=trim(dfuns(j)), iostat=io, &
                        status='old')
                  do j = 1, size(lambdas)
                     read (11, *, iostat=io) rdum(1:3)
                     lambdas(j) = rdum(1)
                     epsilons(j) = cmplx(rdum(2), rdum(3), kind(epsilons))
                  end do
                  close (11)
                  call interp1(lambdas, epsilons, &
                               wavelen, escat(i, k, :))
                  deallocate (lambdas, epsilons)
               end if
            else
               write (*, '(A,A,A,A,i3,A,i3,A,i3)') myname, &
                  '> ERROR: Unrecognised label ''', &
                  trim(labels(i, k)), &
                  ''' of len=', len(trim(labels(i, k))), &
                  ' for scatID=', i, &
                  ' volID=', k - 1
               STOP
            end if
         end do scatterers
      end do shells
      !
   end subroutine calcEpsilon
   !
   !*****************************************************************
   !
   subroutine calcGridPoints(points)
      !
      real(8), intent(inout) :: points(3, nGridPoints)
      character(*), parameter :: myName = 'calcGridPoints'
      integer :: ix, iy, iz, npts
      real(8) :: r(3), dr(3)
      !
      do ix = 1, 3
         if (gridNbins(ix) > 0) then
            dr(ix) = (gridUB(ix) - gridLB(ix))/dble(gridNbins(ix))   !step along each axis
         else
            dr(ix) = 0.0d0
         end if
      end do
      !
      npts = 0
      do ix = 0, gridNbins(1)
         r(1) = gridLB(1) + dble(ix)*dr(1)
         do iy = 0, gridNbins(2)
            r(2) = gridLB(2) + dble(iy)*dr(2)
            do iz = 0, gridNbins(3)
               r(3) = gridLB(3) + dble(iz)*dr(3)
               npts = npts + 1
               points(:, npts) = r
            end do
         end do
      end do
      !
   end subroutine calcGridPoints
   !
   !*****************************************************************
   !
   subroutine sentence2words(sentence, words, nwords_)
      !
      character(*), intent(inout) :: sentence
      character(len=64), intent(inout) :: words(:)
      integer, optional, intent(out) :: nwords_
      integer :: nwords, i
      !
      nwords = 0
      words = ''
      do while (len_trim(sentence) > 0 .and. nwords < size(words))
         nwords = nwords + 1
         sentence = adjustl(sentence) ! remove leading spaces if needed
         i = index(sentence, ' ')  ! look for 1st non-leading space
         words(nwords) = trim(sentence(:i))
         if (words(nwords) (1:1) == '#') then ! discard comment
            words(nwords) = ''
            exit ! ignore all subsequent words
         else
            sentence(:i) = ' ' ! blank first sequence of characters
         end if
      end do
      !
      if (present(nwords_)) nwords_ = nwords
      !
   end subroutine sentence2words
!********************************************************
   subroutine dumpNFs2TXTFile(filename, incidences, Epower, &
                              wavelen, work, Ef, p_label)
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      CHARACTER(LEN=32), intent(in) :: filename
      real(8), intent(in) :: wavelen(:), incidences(:, :), work(:, :, :, :, :)
      integer, intent(in) :: Epower, p_label(:, :)
      logical, intent(in) :: Ef
      ! Local variables
      character(*), parameter :: myName = 'dumpNFs2TXTFile'
      real(8) :: rdum(3), abs_Ave
      real(8), dimension(3, 4) :: Ave
      integer ::  idum, ndum, odum, pdum
      integer ::  nGridPoints

      nGridPoints = size(work, 3)

      open (11, file=filename, status='replace')
      if (size(incidences, 2) == 1) then  !*******
         if (Epower == 0) then
            if (Ef) then
               write (11, *) 'lambda, x, y, z, ScatID, volID, Re(Ex_total), Im(Ex_total), Re(Ey_total), Im(Ey_total),&
                     &Re(Ez_total), Im(Ez_total), Re(Ex_sca), Im(Ex_sca), Re(Ey_sca ), &
                     &Im(Ey_sca), Re(Ez_sca) , Im(Ez_sca)'
            else
               write (11, *) 'lambda, x, y, z, ScatID, volID, Re(Bx_total), Im(Bx_total), Re(By_total), Im(By_total),&
                     &Re(Bz_total), Im(Bz_total), Re(Bx_sca), Im(Bx_sca), Re(By_sca ), &
                     &Im(By_sca), Re(Bz_sca) , Im(Bz_sca)'
            end if
         else
            if (Ef) then
               write (11, *) 'lambda, x,y,z, ScatID, volID, Inc1: |E1|^n'
            else
               write (11, *) 'lambda, x,y,z, ScatID, volID, Inc1: |B1|^n'
            end if
         end if
         do odum = 1, size(wavelen)
            do ndum = 1, nGridpoints
               rdum(1) = sqrt(sum(work(1:3, 2:3, ndum, odum, 1)**2))
               if (rdum(1) == rdum(1)) then ! catches out NaN  ***
                  if (Epower == 0) then
                     write (11, '(f7.2,1x,3(1x,es15.4E2),(i4,1x,i2), 10000(1x,ES24.17E2,1x,ES24.17E2))') &
                        wavelen(odum), work(1:3, 1, ndum, 1, 1), p_label(ndum, 1:2), &
                        ((work(idum, 2, ndum, odum, pdum), work(idum, 3, ndum, odum, pdum), idum=1, 3), &
                         (work(idum, 4, ndum, odum, pdum), work(idum, 5, ndum, odum, pdum), idum=1, 3), &
                         pdum=1, size(incidences, 2))
                  else
                     write (11, '(f7.2,1x,3(1x,es15.4E2),(i4,1x,i2),1000(1x,ES24.17E2))') &
                        wavelen(odum), work(1:3, 1, ndum, 1, 1), p_label(ndum, 1:2), &
                        ((sqrt(sum(work(1:3, 2:3, ndum, odum, pdum)**2)))**Epower, pdum=1, size(incidences, 2))

                  end if

               elseif (verb > 1) then
                  write (*, '(A,A,3(1x,es10.4E1))') myname, &
                     '> enhancement is NaN for r=', work(1:3, 1, ndum, odum, 1)
               end if  !***
            end do
         end do
         close (11)
      else !*******
         if (Epower == 0) then
            if (Ef) then
               write (11, *) 'lambda, x,y,z, ScatID, volID, Ave[Re(Ex_total), Im(Ex_total), Re(Ey_total), Im(Ey_total),&
                       &Re(Ez_total), Im(Ez_total), Re(Ex_sca), Im(Ex_sca), Re(Ey_sca ), &
                       &Im(Ey_sca), Re(Ez_sca) , Im(Ez_sca)], Inc1: ..., Inc2: ..., ...'
            else
               write (12, *) 'lambda, x,y,z, ScatID, volID, Ave[Re(Bx_total), Im(Bx_total), Re(By_total), Im(By_total),&
                       &Re(Bz_total), Im(Bz_total), Re(Bx_sca), Im(Bx_sca), Re(By_sca ), &
                       &Im(By_sca), Re(Bz_sca) , Im(Bz_sca)], Inc1: ..., Inc2: ..., ...'
            end if
         else
            if (Ef) then
               write (11, *) 'lambda,x,y,z, ScatID, volID, Ave(|E|^n), Inc1:|E1|^n, Inc2:|E1|^n ,..'
            else
               write (12, *) 'lambda,x,y,z, ScatID, volID, Ave(|B|^n), Inc1:|B1|^n, Inc2:|B1|^n ,..'
            end if
         end if

         do odum = 1, size(wavelen)
            do ndum = 1, nGridpoints
               rdum(1) = sqrt(sum(work(1:3, 2:3, ndum, odum, 1)**2))
               !
               if (rdum(1) == rdum(1)) then ! catches out NaN  ***
                  if (Epower == 0) then
                     Ave = 0

                     do pdum = 1, size(incidences, 2)
                        Ave(1:3, 1:2) = Ave(1:3, 1:2) + work(1:3, 2:3, ndum, odum, pdum)*incidences(4, pdum)
                        Ave(1:3, 3:4) = Ave(1:3, 3:4) + work(1:3, 4:5, ndum, odum, pdum)*incidences(4, pdum)
                     end do
                     write (11, '(f7.2,1x,3(1x,es15.4E2),(i4,1x,i2),1500(1x,ES24.17E2,1x,ES24.17E2))') &
                        wavelen(odum), work(1:3, 1, ndum, 1, 1), p_label(ndum, 1:2), &
                        (Ave(idum, 1), Ave(idum, 2), idum=1, 3), &
                        (Ave(idum, 3), Ave(idum, 4), idum=1, 3), &
                        ((work(idum, 2, ndum, odum, pdum), work(idum, 3, ndum, odum, pdum), idum=1, 3), &
                         (work(idum, 4, ndum, odum, pdum), work(idum, 5, ndum, odum, pdum), &
                          idum=1, 3), pdum=1, size(incidences, 2))
                  else

                     abs_Ave = 0
                     do pdum = 1, size(incidences, 2)
                        abs_Ave = abs_Ave + ((sqrt(sum(work(1:3, 2:3, ndum, odum, pdum)**2)))**Epower)*incidences(4, pdum)
                     end do
                     write (11, '(f7.2,1x,3(es15.4E2,1x),(i4,1x,i2),700(1x,ES24.17E2))') &
                        wavelen(odum), work(1:3, 1, ndum, 1, 1), p_label(ndum, 1:2), &
                        abs_Ave, &
                        ((sqrt(sum(work(1:3, 2:3, ndum, odum, pdum)**2)))**Epower, pdum=1, size(incidences, 2))

                  end if
               elseif (verb > 1) then
                  write (*, '(A,A,3(1x,es10.4E1))') myname, &
                     '> enhancement is NaN for r=', work(1:3, 1, ndum, odum, 1)
               end if  !***
            end do
         end do
         close (11)

      end if  !*******

   end subroutine dumpNFs2TXTFile

!*******************************************************
   subroutine dumpNFs2HDF5File(fname, groupname, filename, incidences, Epower, &
                               wavelen, work, p_label)

      use HDFfive, only: h5_wrt2file
      !---------------------------------------------------
      ! Start of variable declarations.
      !---------------------------------------------------
      CHARACTER(LEN=64), intent(in) :: fname
      CHARACTER(LEN=32), intent(in) :: groupname, filename
      real(8), intent(in) :: wavelen(:), incidences(:, :), work(:, :, :, :, :)
      integer, intent(in) :: Epower, p_label(:, :)
      !
      ! Local variables
      character(*), parameter :: myName = 'dumpNFs2HDF5File'
      real(8) :: rdum(3), abs_Ave
      real(8), dimension(3, 4) :: Ave
      integer ::  idum, jdum, kdum, mdum, ndum, odum, pdum
      real(8), allocatable ::  work6(:, :)
      integer:: nGridPoints

      nGridPoints = size(work, 3)

      if (allocated(work6)) deallocate (work6)

      if (size(incidences, 2) == 1) then  !*******

         if (Epower == 0) then
            allocate (work6(nGridpoints*size(wavelen), 18))
            if (filename == 'map_E') then
               attribute = 'lambda, x, y, z, ScatID, volID, Re(Ex_total), Im(Ex_total), Re(Ey_total), Im(Ey_total),&
                      &Re(Ez_total), Im(Ez_total), Re(Ex_sca), Im(Ex_sca), Re(Ey_sca ), &
                      &Im(Ey_sca), Re(Ez_sca) , Im(Ez_sca)'
            else
               attribute = 'lambda, x, y, z, ScatID, volID, Re(Bx_total), Im(Bx_total), Re(By_total), Im(By_total),&
                      &Re(Bz_total), Im(Bz_total), Re(Bx_sca), Im(Bx_sca), Re(By_sca ), &
                      &Im(By_sca), Re(Bz_sca) , Im(Bz_sca)'
            end if
         else
            allocate (work6(nGridpoints*size(wavelen), 7))
            if (filename == 'map_E') then
               attribute = 'lambda, x, y, z, ScatID, volID, Inc1: |E1|^n'
            else
               attribute = 'lambda, x, y, z, ScatID, volID, Inc1: |B1|^n'
            end if
         end if

         work6 = 0
         jdum = 0
         do odum = 1, size(wavelen)
            do ndum = 1, nGridpoints
               rdum(1) = sqrt(sum(work(1:3, 2:3, ndum, odum, 1)**2))
               if (rdum(1) == rdum(1)) then ! catches out NaN  ***
                  jdum = jdum + 1
                  if (Epower == 0) then
                     work6(jdum, 1) = wavelen(odum)
                     work6(jdum, 2:4) = work(1:3, 1, ndum, 1, 1)
                     work6(jdum, 5:6) = p_label(ndum, 1:2)
                     work6(jdum, 7:8) = work(1, 2:3, ndum, odum, 1)
                     work6(jdum, 9:10) = work(2, 2:3, ndum, odum, 1)
                     work6(jdum, 11:12) = work(3, 2:3, ndum, odum, 1)

                     work6(jdum, 13:14) = work(1, 4:5, ndum, odum, 1)
                     work6(jdum, 15:16) = work(2, 4:5, ndum, odum, 1)
                     work6(jdum, 17:18) = work(3, 4:5, ndum, odum, 1)

                  else
                     work6(jdum, 1) = wavelen(odum)
                     work6(jdum, 2:4) = work(1:3, 1, ndum, 1, 1)
                     work6(jdum, 5:6) = p_label(ndum, 1:2)
                     work6(jdum, 7) = ((sqrt(sum(work(1:3, 2:3, ndum, odum, 1)**2)))**Epower)
                  end if

               elseif (verb > 1) then
                  write (*, '(A,A,3(1x,es10.4E1))') myname, &
                     '> enhancement is NaN for r=', work(1:3, 1, ndum, odum, 1)
               end if  !***
            end do
         end do
         call h5_wrt2file(fname, groupname, filename, work6, attribute)
      else !*******

         if (Epower == 0) then
            allocate (work6(nGridpoints*size(wavelen), 1 + 3 + 2 + 12 + 12*size(incidences, 2)))
            if (filename == 'map_E') then
               attribute = 'lambda, x,y,z, ScatID, volID, Ave[Re(Ex_total), Im(Ex_total), Re(Ey_total), Im(Ey_total),&
                         &Re(Ez_total), Im(Ez_total), Re(Ex_sca), Im(Ex_sca), Re(Ey_sca ), &
                         &Im(Ey_sca), Re(Ez_sca) , Im(Ez_sca)], Inc1: ..., Inc2: ..., ...'
            else
               attribute = 'lambda, x,y,z, ScatID, volID, Ave[Re(Bx_total), Im(Bx_total), Re(By_total), Im(By_total),&
                         &Re(Bz_total), Im(Bz_total), Re(Bx_sca), Im(Bx_sca), Re(By_sca ), &
                         &Im(By_sca), Re(Bz_sca) , Im(Bz_sca)], Inc1: ..., Inc2: ..., ...'
            end if
         else
            allocate (work6(nGridpoints*size(wavelen), 1 + 3 + 2 + 1 + size(incidences, 2)))
            if (filename == 'map_E') then
               attribute = 'lambda,x,y,z,ScatID, volID,Ave(|E|^n), Inc1:|E1|^n, Inc2:|E1|^n ,..'
            else
               attribute = 'lambda,x,y,z,ScatID, volID,Ave(|B|^n), Inc1:|B1|^n, Inc2:|B1|^n ,..'
            end if
         end if

         work6 = 0
         jdum = 0
         do odum = 1, size(wavelen)
            do ndum = 1, nGridpoints
               rdum(1) = sqrt(sum(work(1:3, 2:3, ndum, odum, 1)**2))
               !
               if (rdum(1) == rdum(1)) then ! catches out NaN  ***
                  jdum = jdum + 1
                  if (Epower == 0) then
                     Ave = 0
                     do pdum = 1, size(incidences, 2)
                        Ave(1:3, 1:2) = Ave(1:3, 1:2) + work(1:3, 2:3, ndum, odum, pdum)*incidences(4, pdum) !Ave(E_total)
                        Ave(1:3, 3:4) = Ave(1:3, 3:4) + work(1:3, 4:5, ndum, odum, pdum)*incidences(4, pdum) !Ave(E_sca)
                     end do

                     work6(jdum, 1) = wavelen(odum)
                     work6(jdum, 2:4) = work(1:3, 1, ndum, 1, 1)
                     work6(jdum, 5:6) = p_label(ndum, 1:2)
                     mdum = 7
                     do idum = 1, 3, 2
                        do kdum = 1, 3
                           work6(jdum, mdum:mdum + 1) = Ave(kdum, idum:idum + 1)
                           mdum = mdum + 2
                        end do

                     end do

                     do pdum = 1, size(incidences, 2)
                        do idum = 2, 4, 2
                           do kdum = 1, 3
                              work6(jdum, mdum:mdum + 1) = work(kdum, idum:idum + 1, ndum, odum, pdum)
                              mdum = mdum + 2
                           end do

                        end do
                     end do

                  else

                     abs_Ave = 0
                     do pdum = 1, size(incidences, 2)
                        abs_Ave = abs_Ave + ((sqrt(sum(work(1:3, 2:3, ndum, odum, pdum)**2)))**Epower)*incidences(4, pdum)
                     end do
                     work6(jdum, 1) = wavelen(odum)
                     work6(jdum, 2:4) = work(1:3, 1, ndum, 1, 1)
                     work6(jdum, 5:6) = p_label(ndum, 1:2)
                     work6(jdum, 7) = abs_Ave

                     do pdum = 1, size(incidences, 2)
                        work6(jdum, pdum + 7) = (sqrt(sum(work(1:3, 2:3, ndum, odum, pdum)**2)))**Epower
                     end do

                  end if

               elseif (verb > 1) then
                  write (*, '(A,A,3(1x,es10.4E1))') myname, &
                     '> enhancement is NaN for r=', work(1:3, 1, ndum, odum, 1)
               end if  !***
            end do
         end do
         call h5_wrt2file(fname, groupname, filename, work6, attribute)
      end if !****************

   end subroutine dumpNFs2HDF5File
!********************************************************
   !
   function countLines(filename) result(nlines)
      !
      ! Function for counting lines in file
      !
      character(*), intent(in) :: filename
      integer :: nlines, io
      !
      open (unit=11, file=filename, iostat=io, status='old')
      if (io /= 0) then
         nlines = -1
      else
         nlines = 0
         do
            read (11, *, iostat=io)
            if (io /= 0) exit
            nlines = nlines + 1
         end do
         close (11)
      end if
      !
   end function countLines

   !
end program termsProgram
