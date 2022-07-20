!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  The is a clump-finding algorithm that is based upon the MRB's disc
!  finding algorithm.  We take the densest particle, and find all the bound
!  particles.  We then take the next densest particle that is not part of
!  the clump, and repeat.
!
!  REFERENCES: None
!
!  OWNER: James Wurster
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, part, prompting
!+
!--------------------------------------------------------------------------

module analysis
 use dim,          only: maxptmass, maxvxyzu, mhd, maxp_hard ! maximum number of ptmasses, values of vxyzu, mhd => using mhd, maximum number of particles
 use part,         only: Bxyz, xyzmh_ptmass, vxyz_ptmass, nptmass, poten, ihacc, periodic, iorig ! Vectors of B-field, ptmass details, number of ptmasses, particle potentials, accretion radius, period boundaries accounted for, positions of original particles (indices)
 use sortutils,    only: indexx ! tool for sorting particles according to ascending parameter

 implicit none

 character(len=20), parameter, public :: analysistype = 'clumpfind_kat'

 public  :: do_analysis

 integer, parameter  :: nclumpmax = 1000 ! maximum number of individual clumps which can be identified
 integer :: idclumpold(maxp_hard), idclump(maxp_hard), iorigold(maxp_hard) ! id of clumps in previous dump, id of clumps in current dump, indices of particles originally
 integer :: idclumpsinkold(maxptmass), idclumpsink(maxptmass) ! id of clumps containing sinks in previous and current dumps respectively
 integer :: idclumpmax ! clump running total
 integer :: nmin = 57 ! minimum number of particles required for a pure-gas clump
 integer :: idtmp = -10 ! temporary id
 integer :: dims = 3 ! Number of dimensions

 logical :: firstcall = .true., init_conds = .true. ! tells us whether we are on first dump analysed or not

 real :: rhominpeak_cgs = 1.0d-17 ! particles with rho < rhominpeak cannot be clump particle lead (cgs units)
 real :: rhominbkg_cgs = 1.0d-18!d-30!-18 ! particles with rho < rhominbkg cannot be considered for clump membership (cgs units)
 real :: cloudmax_pc = 200. ! maximum extent of cloud (in pc)
 real :: cloudmax2, cloudmax ! square of cloudmax, and cloudmax
 real :: ekin_coef = 2. ! clump must satisfy the condition ekin_coef*ekin + epot < 0 => can change ekin_coef
! real :: yexc = 0. ! remnant from analysis of cloud collision simulations by james - might remove this

 ! define a new parameter type: sphclump
 type sphclump
    integer :: nmember, pointmass(nclumpmax), ID ! number of members in the clump, number of sinks, ID of the clump
    integer :: conv_ID_v(nclumpmax + 1), conv_ID_a(nclumpmax + 1) ! ID of clumps a given clump with converge with on a timescale < tdyn (two different arrays depending on the prescription of tdyn used) - index 1 is for the final assigned value (i.e., clump which is collapsing onto this one first)
    integer, allocatable :: partIDs(:)
    real, dimension(3) :: r, v, dr, a ! clump coordinates, clump velocity vector, separation, acceleration?
    real :: mass, size ! mass of clump, physical extent of clump
    real :: kinetic, potential, thermal, magnetic, bound ! energies and bound state of the clump
    real :: time_ff, time_dyn, time_dynv, time_dyna
    real :: tform, tdead
    real :: ptmass_pos(3, nclumpmax)
    ! real :: 
    real, allocatable :: particle_positions(:, :) ! number should be the same as nmembermax defined below
    logical :: col(2) = .false., sink_only = .false., reassigned = .false. ! Logical which tells you whether clump will collide with another one on a timescale < tff, whether clump is sink only, whether sink_only clump was reassigned to another clump
    logical :: SURV = .FALSE., NEW = .TRUE.
 end type sphclump

 ! assign parameter to the type sphclump, with a length equal to the number of possible clumps
 type(sphclump) :: clump(nclumpmax)

!----- New additions for working out t_dyn -----!

 integer :: nptmass_prev

 character(len=50) :: dump_cur, dump_prev ! Current and previous dump strings

 real :: vxyzu_prev(3, maxp_hard), vxyz_ptmass_prev(3, maxptmass), time_prev ! Previous particle and pointmass velocities, previous dumpfile time
 real :: xyzh_init(4, maxp_hard), vxyzu_init(4, maxp_hard) ! Initial conditions from *_00000 dumpfile
 real :: pi

 logical :: time_comp = .false., prev_file = .false. ! Logical to check whether two dumpfiles are supplied => required for t_dyn vs t_ff comparison
 integer :: clumps_cur(nclumpmax), clumps_prev(nclumpmax) !!!!!!!!!!
 !----------------------------------------------------------------------------!

 private ! not shared beyond this module

contains

!----------------------------------------------------------------------------------------------------------------------
! Module contains 3 subroutines which can be used: do_analysis (which draws on get_membership and reset_clump), obtain_membership (works out whether a particle belongs to a clump) and reset_clump (resets clump properties)

subroutine do_analysis(dumpfile, num, xyzh, vxyzu, particlemass, npart, time, iunit)
 use part,          only: rhoh ! works out the particle density based on h
 use physcon,       only: pc ! pc units
 use units,         only: udist, umass, utime, unit_energ ! unit distance, unit mass
 use prompting,     only: prompt

 ! parameters which are fed in/out
 character(len=*), intent(in) :: dumpfile ! file to analyse

 integer, intent(in) :: num, npart, iunit ! , number of particles, unit to open

 real, intent(in) :: xyzh(:,:) ! particle positions
 real, intent(inout) :: vxyzu(:,:) ! velocity vectors
 real, intent(in) :: particlemass, time ! particle masses, time in simulation

! character(len=*), optional :: dumpfile2

 ! parameters contained in SR:
 character(len=128) :: filename, prefix ! strings for saving data
 character(len=8) :: fmt, str_count, sink_count, nparts_str ! format descriptor - definitions for when we are lookping over density thresholds for the same dump
 character(len=5) :: dump_index ! string store for current dump
 character(len=8) :: id_fmt ! format descriptor for clump ID
 character(len=6) :: sink_fmt ! format descriptor for sink ID
 !character(len=128) :: init_dump ! Name of initial dumpfile

 integer :: ii(3), i, j, k, io, iclump, iclumptmp, idx, izero, ndense, tmpctr ! ii = particle being analysed, ijk = counting integers, , clump index, temporary "", dumpfile index, , number of particles located in dense regions => potential clump members
 integer :: lst_rho(npart), lst_x(npart), lst_y(npart), lst_z(npart), idense(npart), iis(3, maxptmass) ! particle indices in order of density, "" xyz-position, dense particle indices, indices of point masses, added 3D indexing to iis
 integer :: lst_rad(npart), drad, dradmin(maxptmass) ! radial distance as opposed to only y
 integer :: comp_count ! counter for keeping track of the number of comparison analyses done for a single dump
 integer :: dump_no_cur, dump_no_prev ! Number of current and previously analysed dumpfile
 integer :: most_bound, Nsync, sinkIDs(nclumpmax), sink_locs(3, maxptmass), resctr ! ID of clump a sink is most bound to, number of sinks per clump, IDs of sinks, locations of sinks, number of sinks reassigned to another clump
 integer :: colv, cola ! For writing collision status of clump to file

 logical :: keep_searching = .false., iexist ! logical value for checking whether to keep looking through particles, check whether a parameter exists
 logical :: init_dump = .false. ! Is this the initial dumpfile?

 real :: rhomin_peak, rhomin_bkg, xi, yi, zi, hi, dx, dy, dz, rtmp ! peak and background densities in code units, specific particle coordinates, xyz-separation, temporary value holder
 real :: dxmin(maxptmass), dymin(maxptmass), dzmin(maxptmass) ! Minimum separation of particle from a sink in xyz
 real :: rho(npart), rho_dense(npart) ! densities of all particles, densities of all particles considered to be in dense regions, particle r positions, minimum separation between a particle and a point mass
 real :: rhomax, rhomin, rhomean, rhomean_dense, rhobkg_range(5), rho_stddev ! maximum, minimum, and mean densities achieved (all particles and dense particles), array for storing a range of background density thtresholds to use, and density standard deviation
 real :: tff_cm(npart), m_cm(npart), clump_T(npart), clump_U(npart), rho_cont(2), rpos(3,npart) ! Freefall timescale, mass, kinetic, potential of a clump per particle, array for storing density contour max and min values, and particle xyz positions
 real :: tdyn_a(nclumpmax, nclumpmax), tdyn_v(nclumpmax, nclumpmax), rij(3), vij(3), aij(3), rijabs, vijabs, adotr ! dynamical timescale of a clump wrt the other clumps based on accelerations and velocities, vector clump separations and relative velocities and accelerations
! absolute value of rij, dot product of a and r
 real :: pot_tmp, etot_tmp(2) ! temporary store of clump potentials and total energies
 real :: part_vx, part_vx_init, part_vy, part_vy_init, part_vz, part_vz_init, part_v_tot, part_v_init_tot ! Arrays for initial and current particle velocities
 real :: mwms(3), ptmasssum, mass_sum, vcom(1:3) ! mass weighted momentum sum, sum of point masses, mass sum of system, centre of velocity
 real :: colvrs(3), colvvs(3), colvas(3), colars(3), colavs(3), colaas(3)
! TEST ZONE -------------------------------------------------------------------
 ! TEST ZONE -------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------------------------------------

 if (firstcall) then
     clumps_prev = 0
 endif

 clumps_cur = 0

! SHIFT TO THE CENTRE-OF-MOMENTUM FRAME:
 mwms = 0.! Mass-weighted momentum sum
 ! Contribution from gas particles
 do i = 1, npart
     do j = 1, 3
         mwms(j) = mwms(j) + vxyzu(j, i)*particlemass*particlemass
     enddo
 enddo
 ! Contribution from sink particles
 ptmasssum = 0.
 if (nptmass > 0) then
     do i = 1, nptmass
         ptmasssum = ptmasssum + xyzmh_ptmass(4, i)
         do j = 1, 3
             mwms(j) = mwms(j) + vxyz_ptmass(j, i)*xyzmh_ptmass(4, i)*xyzmh_ptmass(4, i)
         enddo
     enddo
 endif
 ! Total system mass
 mass_sum = npart*particlemass + ptmasssum*nptmass
 mwms = mwms/mass_sum
 ! System COV
 vcom = mwms/mass_sum
 ! New velocities
 do i = 1, npart
     vxyzu(1:3, i) = vxyzu(1:3, i) - vcom(1:3)
 enddo
 do i = 1, nptmass
     vxyz_ptmass(1:3, i) = vxyz_ptmass(1:3, i) - vcom(1:3)
 enddo
! SHIFT TO THE CENTRE-OF-MOMENTUM FRAME

 if (init_conds) then
     !print*, init_conds
     call prompt('Initial dumpfile name for velocity comparisons (must provide initial dumpfile followed by dumpfile &
                 & you want analysed)', init_dump)
     if (init_dump) then
         xyzh_init(:, 1:npart) = xyzh(:, 1:npart)
         vxyzu_init(1:3, 1:npart) = vxyzu(1:3, 1:npart)
         init_conds = .false.
         return
     else
         init_conds = .false.
     endif
 endif

 ! Establish whether a previous comparison dumpfile has been supplied to allow particle accelerations to be determined - assumes first file in command line is the reference file and second file is the file to be analysed
 if ((.not. prev_file) .and. (.not. init_conds)) then
     call prompt('Do you wish to compare a clumps t_dyn with t_ff as an additional clump measure? (will not work unless &
               & a previous comparison dumpfile is supplied)', time_comp)

     if (time_comp) then
         prev_file = .true. ! Previous file has been supplied
         time_prev = time
         vxyzu_prev(:, 1:npart) = vxyzu(:, 1:npart)
         nptmass_prev = nptmass
         vxyz_ptmass_prev(:, 1:nptmass) = vxyz_ptmass(:, 1:nptmass)
         return ! Move onto next dumpfile for proper analysis
     endif
 endif         
 
 !fmt = '(I2.2)' ! Integer of width 2 with zeros on the left
 fmt = '(I6)'
 comp_count = 0
 rhomax = tiny(rhomax)
 rhomin = huge(rhomin)
 rhomean = 0
 rhomean_dense = 0
 rho_stddev = 0
 tdyn_a = huge(tdyn_a)
 tdyn_v = huge(tdyn_v)
 izero = 0
 iclump = 0
 idense = 0
 rhomin_peak = rhominpeak_cgs*udist**3/umass
 rhomin_bkg = rhominbkg_cgs*udist**3/umass
 cloudmax = cloudmax_pc*pc/udist
 cloudmax2 = cloudmax**2
 tmpctr = 0
 ndense = 0

 if (firstcall) then ! using parameters from the first dumpfile
     idclumpold = 0
     idclumpsinkold = 0
     idclumpmax = 0
     iorigold(1:npart) = iorig(1:npart)
     inquire(file='clumpfind_KEcoef.in', exist=iexist) ! checks whether there is a file to read in which provides a value for the  clumpfind_KEcoef - could make this a prompt?
     if (iexist) then
         open(unit=122, file = 'clumpfind_KEcoef.in')
         read(122,*, iostat=io) rtmp
         if (io == 0) ekin_coef = rtmp ! read in a value for the coefficient
         close(122)
     endif
 endif

 idclump = 0
 idclumpsink = 0
 idx = index(dumpfile, '_')
 prefix = dumpfile(1:idx-1)
 print*, 'Groups of particles are considered a clump if', ekin_coef, '*ekin + epot < 0'

 ! First exclude background particles from the calculation of mean particle density (they will bring it down too low)
 do i = 1, npart
     yi = xyzh(2, i)
     hi = xyzh(4, i)
     if (hi > 0.) then
         rho(i) = rhoh(hi, particlemass) ! found density of particle
         rhomax = max(rhomax, rho(i))
         rhomin = min(rhomin, rho(i))
         if (rho(i) > rhomin_bkg) then
             ndense = ndense + 1 ! have identified a particle in a region above minimum density threshold
             rhomean = rhomean + log10(rho(i))
         endif
     endif
 enddo

 rho_cont(1) = rhomin
 rho_cont(2) = rhomax

 ! Save the maximum and minimum density values for setting contour plots in splash
 open (unit=26,file='splash.contours_lims', status='unknown')
 do i = 1, 2
        write(26,'(Es18.6, a)') rho_cont(i)
 enddo
 close(26)

 print*, 'Initial background threshold; ', rhominbkg_cgs , 'g cm^-3'
 print*, 'Initial peak threshold: ', rhominpeak_cgs, 'g cm^-3'

 if (ndense == 0.) return
 rhomean = rhomean/real(ndense)
 rhomean = 10**(rhomean)
 rhomin_bkg = rhomean ! Set the background threshold to rhomean

 rhomean = 0 ! reset for new definitions
 ndense = 0

 ! Now calculate the mean using particles which were likely part of the original cloud

 do i = 1, npart
    hi = xyzh(4, i)
    if (hi > 0) then
        rho(i) = rhoh(hi, particlemass)
        if (rho(i) > rhomin_bkg) then
            ndense = ndense + 1 ! have identified a particle in a region above minimum density threshold
            idense(ndense) = i ! save the index of this particle
            rho_dense(ndense) = rho(i) ! save density of this particle
            do k = 1, 3
                rpos(k, ndense) = xyzh(k, i) ! save location of this particle (in this case y-coordinate but it can be changed - also requires modifying lst_y though)
            enddo
            rhomean_dense = rhomean_dense + log10(rho(i))
        endif
     endif
 enddo
 print*, 'Dense gas located; ', ndense, ' particles, out of ', npart, ' total particles'
 
 rhomean_dense = rhomean_dense/real(ndense)
 rhomean_dense = 10**(rhomean_dense)
 rhomin_peak = rhomean_dense ! Set the minimum density required for a particle to be considered as clump lead via log mean of dense particle densities

 print*, 'Background density threshold reset to ', rhomin_bkg*(umass/udist**3), 'g cm^-3'
 print*, 'Particles lead if densities are higher than ', rhomin_peak*(umass/udist**3), 'g cm^-3'

 ! Sort dense particles by density and position
 call indexx(ndense, rho_dense, lst_rho)
 call indexx(ndense, rpos(1, :), lst_x)
 call indexx(ndense, rpos(2, :), lst_y)
 call indexx(ndense, rpos(3, :), lst_z)

 ! Moving along the y-axis, we find the closest gas particle to the sink
 iis(:, :) = 0 ! sink index
 dxmin = huge(dxmin)
 dymin = huge(dymin) ! set dymin to a very large value of the same type as dymin for comparison
 dzmin = huge(dzmin)
 do j = 1, ndense ! look at all these dense particles
    do i = 1, nptmass ! for a given particle, work out the closest sink
         !dy = abs(xyzmh_ptmass(2, i) - xyzh(2, idense(lst_y(j))))
         dx = abs(xyzmh_ptmass(1, i) - xyzh(1, idense(lst_rho(j))))
         dy = abs(xyzmh_ptmass(2, i) - xyzh(2, idense(lst_rho(j))))
         dz = abs(xyzmh_ptmass(3, i) - xyzh(3, idense(lst_rho(j))))
         drad = (dx*dx + dy*dy + dz*dz)**0.5
        ! Position of particle j in lst_y => find it's position along y-axis in the cloud and use this to index the particle identity from idense
        ! Moving up lst_y incrementally means you get the indices of original particles based on their location e.g., the particle closest to the centre is particle x, with an index in lst_y of 1. Use this to find original particle index from the original list of particles.
         if (dx < dxmin(i)) then ! Ordered list in 3D
             iis(1, i) = j ! save the index of the particle closest to the sink
             dxmin(i) = dx
         endif
         if (dy < dymin(i)) then
             iis(2, i) = j ! save the index of the particle closest to the sink
             dymin(i) = dy
         endif
         if (dz < dzmin(i)) then
             iis(3, i) = j ! save the index of the particle closest to the sink
             dzmin(i) = dz
         endif
     enddo
 enddo
 print*, 'Found gas nearest to each sink'

 ! Time to look for clumps around sinks
 if (nptmass > 0) then
    do i = 1, nptmass
        ! Save nptmass
        iclump = iclump + 1 ! update the clump index
        xi = xyzmh_ptmass(1, i)
        yi = xyzmh_ptmass(2, i)
        zi = xyzmh_ptmass(3, i) ! identify sink position

        ! initialise the clump
        call reset_clump(iclump)
        ! localise clump properties onto that of the point mass
        clump(iclump)%r(1:3) = xyzmh_ptmass(1:3, i)
        clump(iclump)%v(1:3) = vxyz_ptmass(1:3, i)
        clump(iclump)%mass = xyzmh_ptmass(4, i)
        clump(iclump)%size = xyzmh_ptmass(ihacc, i) ! size set by accretion radius
        clump(iclump)%pointmass(1) = i
        print*, 'Sink mass defining clump = ', clump(iclump)%mass

        ! work out which particles belong to the sink-containing clump:
        call obtain_membership(dumpfile, iis(:, i), xi, yi, zi, particlemass, iclump, npart, lst_x, lst_y, lst_z, idense, ndense,&
                i, xyzh, vxyzu, rho, iis, time) ! i is the sink index, iis is the next particle on
        print*, 'CLUMP SINK', iclump, clump(iclump)%ID, xyzmh_ptmass(4, i), clump(iclump)%nmember, &
                clump(iclump)%kinetic, clump(iclump)%thermal, clump(iclump)%potential, &
                clump(iclump)%mass, clump(iclump)%size, time*(utime/(3600*24*365))
    enddo
    print*, 'Found clumps around sinks'
 else
    print*, 'No sinks found'
 endif

 ! Now we want to find clumps which don't contain sinks
 k = ndense
 if (k > 2) keep_searching = .true.
 do while (keep_searching)
    do while (idclump(idense(lst_rho(k)))/=0 .and. k > 2) ! Initially this is not satisfied as idclump is initialised to zero
        ! This means that the lead particle of a clump has already been identified(?)/is already part of a clump...
        k = k - 1
        if (k <= 2) keep_searching = .false.
    enddo
    i = idense(lst_rho(k)) ! We obtain the index of the most dense particle in the dump
    if (rho(i) > rhomin_peak) then
        iclump = iclump + 1
        ! New clump because we have moved onto next densest particle
        iclumptmp = iclump
        xi = xyzh(1, i)
        yi = xyzh(2, i)
        zi = xyzh(3, i) ! Clump centred on this particle

        ! initialise the clump
        call reset_clump(iclump)
        ! localise the clump properties on the lead particle
        clump(iclump)%r(1:3) = xyzh(1:3, i)
        clump(iclump)%v(1:3) = vxyzu(1:3, i)
        clump(iclump)%nmember = 1 ! this particle in question
        clump(iclump)%mass = particlemass
        if (maxvxyzu >= 4) clump(iclump)%thermal = vxyzu(4, i)*particlemass ! eos conditions
        if (mhd) clump(iclump)%magnetic = 0.5*particlemass*dot_product(Bxyz(1:3, i), Bxyz(1:3, i))/rho(i)
        idclump(i) = idtmp ! I think this means we have located our clump lead

        ! want to identify the clump members
        ii(1) = minloc(abs(idense(lst_x(1:ndense)) - i), 1)
        ii(2) = minloc(abs(idense(lst_y(1:ndense)) - i), 1) ! Minloc finds the index of the minimum value of an array ! might need to change this one back to lst_y?
        ii(3) = minloc(abs(idense(lst_z(1:ndense)) - i), 1)
        ! Take idense and subtract the index of the most dense particle in the dump
        ! Take the absolute value of this, and the smallest of these numbers sets ii (just gives you the index of i in idense)
        ! find which particles belong to the clump set by ii
        call obtain_membership(dumpfile, ii(:), xi, yi, zi, particlemass, iclump, npart, lst_x, lst_y, lst_z, idense, ndense,&
                izero, xyzh, vxyzu, rho, iis, time)
        ! remove the peak clump canidate from the list, or print the clump details
        if (iclump < iclumptmp) then
            idclump(i) = -1
        else
            idclump(i) = clump(iclump)%ID
            print*, 'CLUMP!', iclump, clump(iclump)%ID, rho(i)*umass/udist**3, clump(iclump)%nmember, &
                clump(iclump)%kinetic, clump(iclump)%thermal, clump(iclump)%potential, clump(iclump)%mass, &
                clump(iclump)%size, time*(utime/(3600*24*365))
        endif
    else
        keep_searching = .false. ! no point in continuing the search because the particle does not rise above background density criterion
    endif
 enddo
 resctr = 0
! want to check here if sink is contained with xyz of another clump
 do i = 1, iclump
     etot_tmp = 0
     if (clump(i)%sink_only) then ! for clumps which only have a sink in them - no gas
         do j = 1, iclump
             if (i == j) cycle
             if (.not. clump(j)%sink_only) then ! only do this for clumps which have additional gas members
                 ! Check if sink could be bound to another clump
                 rij = clump(i)%r - clump(j)%r
                 rijabs = dot_product(rij(1:3), rij(1:3))
                 pot_tmp = clump(j)%potential - (clump(i)%mass*clump(j)%mass)/sqrt(rijabs)
                 etot_tmp(2) = clump(i)%kinetic + clump(j)%kinetic + pot_tmp
                 if (etot_tmp(2) < 0.) then
                     resctr = resctr + 1
                     clump(i)%reassigned = .true.
                     clump(i)%SURV = .false.
                     if ((.not. clump(i)%SURV) .and. (.not. clump(i)%NEW)) then
                         clump(i)%tdead = time
                     endif
                     if (etot_tmp(2) < etot_tmp(1)) then
                         most_bound = clump(j)%ID
                         etot_tmp(1) = etot_tmp(2) ! Save this as the minimum comparison value for next time
                     endif
                 endif
             endif
         enddo
         Nsync = count(clump(most_bound)%pointmass /= 0) ! number of sinks reassigned to the clump - where 'most_bound' is the ID of the clump the sink is most bound to
         ! Need to update clump properties to include those of the sinks it has taken on
         clump(most_bound)%pointmass(Nsync + 1) = clump(i)%pointmass(1) ! add sink to the clump it is most bound to (+ 1 as Nsync is automatically 1 for each clump already containing a sink)
         if (clump(most_bound)%pointmass(1) /= 0) then
             clump(most_bound)%ptmass_pos(:, 1) = clump(most_bound)%r
         endif
         ! Save sink position wrt clump:
         clump(most_bound)%ptmass_pos(:, Nsync + 1) = clump(i)%r
         ! Update clump mass:
         clump(most_bound)%mass = clump(most_bound)%mass + clump(i)%mass
         ! Update clump potential:
         rij = clump(i)%r - clump(most_bound)%r
         rijabs = sqrt(dot_product(rij(1:3), rij(1:3)))
         clump(most_bound)%potential = clump(most_bound)%potential - (clump(i)%mass*clump(most_bound)%mass)/sqrt(rijabs)
         ! Update clump kinetic energy:
         clump(most_bound)%thermal = clump(i)%thermal + clump(most_bound)%thermal
         ! Update clump thermal energy:
         clump(most_bound)%magnetic = clump(i)%magnetic + clump(most_bound)%magnetic
         ! Update clump magnetic energy:
         clump(most_bound)%kinetic = clump(i)%kinetic + clump(most_bound)%kinetic
         ! Update clump size:
         clump(most_bound)%size = max(rijabs, clump(most_bound)%size)
         ! Update clump member number:
         clump(most_bound)%nmember = clump(most_bound)%nmember + 1
     endif
 enddo
 ! We have now identified all the clumps initially, but we want to check whether they will collide before they collapse due to free-fall
 do i = 1, iclump
     do j = 1, iclump
         if (i == j) cycle ! Do not compare a clump with itself
         if (clump(i)%reassigned) cycle
         if (clump(j)%reassigned) cycle ! these clumps should have been reassigned...
         rij = clump(i)%r - clump(j)%r
         vij = clump(i)%v - clump(j)%v
         vijabs = sqrt(dot_product(vij(1:3), vij(1:3)))
         rijabs = sqrt(dot_product(rij(1:3), rij(1:3)))
         tdyn_v(i, j) = rijabs/vijabs
 !----- comment back in if we want multiple collision comparisons -----!
         if ((tdyn_v(i, j) < clump(i)%time_ff)) then
             clump(i)%conv_ID_v(j + 1) = clump(j)%ID ! SO as not to redefine the first value which corresponds to the clump we are assessing
         endif
 !---------------------------------------------------------------------!
     enddo
     if ((tdyn_v(i, minloc(tdyn_v(i, 1:idclumpmax), 1)) < clump(i)%time_ff)) then ! clump i will collide with clump j before collapsing
         clump(i)%col(1) = .true.
         clump(i)%conv_ID_v(1) = clump(i)%conv_ID_v(minloc(tdyn_v(i, 1:idclumpmax), 1) + 1)
     endif
     clump(i)%time_dynv = tdyn_v(i, minloc(tdyn_v(i, 1:idclumpmax), 1))
 enddo

 ! We have now identified all the clumps initially, but we want to check whether they will collide before they collapse due to free-fall
 do i = 1, iclump
     do j = 1, iclump
         if (i == j) cycle ! Do not compare a clump with itself
         if (clump(i)%reassigned) cycle
         if (clump(j)%reassigned) cycle ! these clumps should have been reassigned...
         rij = clump(i)%r - clump(j)%r
         aij = clump(i)%a - clump(j)%a
         adotr = dot_product(aij(1:3), rij(1:3))
         rijabs = sqrt(dot_product(rij(1:3), rij(1:3)))
         tdyn_a(i, j) = ((rijabs**2)/(abs(adotr)))**0.5
!----- comment back in if we want multiple collision comparisons -----!
         if ((tdyn_a(i, j) < clump(i)%time_ff)) then
             clump(i)%conv_ID_a(j + 1) = clump(j)%ID ! SO as not to redefine the first value which corresponds to the clump we are assessing
         endif
!---------------------------------------------------------------------!
     enddo
     if ((tdyn_a(i, minloc(tdyn_a(i, 1:idclumpmax), 1)) < clump(i)%time_ff)) then ! clump i will collide with clump j before collapsing
         clump(i)%col(2) = .true.
         clump(i)%conv_ID_a(1) = clump(i)%conv_ID_a(minloc(tdyn_a(i, 1:idclumpmax), 1) + 1)
     endif
     clump(i)%time_dyna = tdyn_a(i, minloc(tdyn_a(i, 1:idclumpmax), 1))
 enddo

 if (iclump == 1) then
     clump(iclump)%time_dynv = 0.
     clump(iclump)%time_dyna = 0.
 endif

 ! Save the clupm IDs for comparison
 do i = 1, iclump
     if (.not. clump(i)%reassigned) then
         clumps_cur(i) = clump(i)%ID
     endif
 enddo

 do i = 1, nclumpmax
     ! Work out whether the identified clump formed in this dump or a not
     if (clumps_cur(i) > 0) then
         do j = 1, nclumpmax
             if (clumps_cur(i) == clumps_prev(j)) then
                 do k = 1, iclump
                     if (clump(k)%ID == clumps_cur(i)) clump(k)%NEW = .false.
                 enddo
             endif
         enddo
     endif

     ! Work out whether any clumps disappeared since the last dump
     if (clumps_prev(i) > 0) then
         print*, clumps_cur(i)
         do j = 1, nclumpmax
             if (clumps_prev(i) == clumps_cur(j)) then
                 do k = 1, iclump
                     if (clump(k)%ID == clumps_prev(i)) clump(k)%SURV = .true.
                 enddo
             endif
         enddo
     endif
 enddo
 do i = 1, iclump
     if ((.not. clump(i)%SURV) .and. (.not. clump(i)%NEW)) clump(i)%tdead = time
     if (clump(i)%NEW) clump(i)%tform = time
 enddo


 ! Write the results to file - taken straight from original script
 !write(str_count, fmt) comp_count
 
!---------------------------------------------------------- USEFUL WRITING TO FILE ----------------------------------------------------------!

! fmt = '(I8)'
! write(nparts_str, fmt) npart
! print*, nparts_str
! read(*,*)

 do i = 1, iclump
     if (clump(i)%SURV) then
         write(filename, '(1a,I5.5,1a)') 'clumpID_',clump(i)%ID,'_param_evol'
         inquire(file=trim(filename),exist=iexist)
         if (iexist) then
             open(26,file=trim(filename), position='append')
         else
             open(26,file=trim(filename))
         endif
         if ((.not. clump(i)%col(1)) .or. (clump(i)%conv_ID_v(1) < 1)) then
            colv = 0
            colvrs = 0.
            colvvs = 0.
            colvas = 0.
         else
            colv = 1
            colvrs = clump(clump(i)%conv_ID_v(1))%r
            colvvs = clump(clump(i)%conv_ID_v(1))%v
            colvas = clump(clump(i)%conv_ID_v(1))%a
         endif
         if ((.not. clump(i)%col(2)) .or. (clump(i)%conv_ID_a(1) < 1)) then
            cola = 0
            colars = 0.
            colavs = 0.
            colaas = 0.
         else
            cola = 1
            colars = clump(clump(i)%conv_ID_a(1))%r
            colavs = clump(clump(i)%conv_ID_a(1))%v
            colaas = clump(clump(i)%conv_ID_a(1))%a
         endif
         Nsync = count(clump(i)%pointmass /= 0)! - 1
         if (Nsync /= 0) then
             sinkIDs(1:Nsync) = clump(i)%pointmass(1:Nsync)
!             do j = 1, Nsync
!                 write(26,'(7Es18.6,1I8)') clump(sinkIDs(j))%r, 0., 0., 0., 0., sinkIDs(j)
!             enddo
         endif
         write(26, '(1I8,15Es18.6,1I8,9Es18.6,1I8,9Es18.6,2I8,3Es18.6,1I8,4Es18.6)') &
            clump(i)%ID, clump(i)%r, clump(i)%v, clump(i)%a, clump(i)%mass, clump(i)%size, & ! 1 - 12
            clump(i)%time_ff, clump(i)%time_dyn, clump(i)%time_dynv, clump(i)%time_dyna, & ! 13 - 16
            clump(i)%conv_ID_v(1), colvrs, colvvs, & ! 17 - 23
            colvas, clump(i)%conv_ID_a(1), colars, & ! 24 - 30
            colavs, colaas, colv, cola, & ! 31 - 38
            clump(i)%tform*utime, clump(i)%tdead*utime, time*utime, clump(i)%pointmass(1), & ! 39 - 42
            clump(i)%kinetic, clump(i)%potential, clump(i)%thermal, clump(i)%magnetic ! 43 - 46


         write(filename, '(1a,I5.5,1a)') 'clump_part_store/clumpID_',clump(i)%ID,'_part_v_t'
         inquire(file=trim(filename),exist=iexist)
         if (iexist) then
             open(27,file=trim(filename), position='append')
         else
             open(27,file=trim(filename))
         endif

         fmt = '(I8)'
         write(nparts_str, fmt) count(clump(i)%partIDs /= 0)
         print*, trim(nparts_str)
         write(27, '(1Es18.6,'//trim(nparts_str)//'I8)') &
            time*utime, clump(i)%partIDs(1:count(clump(i)%partIDs /= 0))

     endif
     colv = 1
     cola = 1
 enddo
! read(*,*)
 close(26)
 close(27)

          




!  write(filename,'(2a)') trim(dumpfile),'_timescale_comparisons'
!  open (unit=26,file=trim(filename))
!  do i = 1, iclump
!      if (.not. clump(i)%col(1)) colv = 0
!      if (.not. clump(i)%col(2)) cola = 0
!      if (.not. clump(i)%sink_only) then
!          write(26, '(1I8,14Es18.6,1I8,9Es18.6,1I8,9Es18.6,2I8)') clump(i)%ID, clump(i)%r, clump(i)%v, clump(i)%a, clump(i)%mass,&
!                             clump(i)%time_ff, clump(i)%time_dyn, clump(i)%time_dynv, clump(i)%time_dyna, &
!                             clump(i)%conv_ID_v(1), clump(clump(i)%conv_ID_v(1))%r, clump(clump(i)%conv_ID_v(1))%v, &
!                             clump(clump(i)%conv_ID_v(1))%a, clump(i)%conv_ID_a(1), clump(clump(i)%conv_ID_a(1))%r, &
!                             clump(clump(i)%conv_ID_a(1))%v, clump(clump(i)%conv_ID_a(1))%a, colv, cola
!      endif
!      colv = 1
!      cola = 1
!  enddo

!  do i = 1,iclump
!      if (.not. clump(i)%sink_only) then
!          write(filename,'(2a,I5.5,1a)') trim(prefix),'_clumpID_',clump(i)%ID,'_particle_positions'
!          open(unit=26,file=trim(filename))
!          Nsync = count(clump(i)%pointmass /= 0)! - 1
!          if (Nsync /= 0) then
!              sinkIDs(1:Nsync) = clump(i)%pointmass(1:Nsync)
!              do j = 1, Nsync
!                  write(26,'(7Es18.6,1I8)') clump(sinkIDs(j))%r, 0., 0., 0., 0., sinkIDs(j)
!              enddo
!          endif
!          do j = 2, clump(i)%nmember - Nsync ! Avoids zero positions being saved (as sinks count as members)
!              part_vx = vxyzu(1, clump(i)%partIDs(j))
!              part_vy = vxyzu(2, clump(i)%partIDs(j))
!              part_vz = vxyzu(3, clump(i)%partIDs(j))
!              part_v_tot = sqrt(part_vx**2 + part_vy**2 + part_vz**2)
    
!              part_vx_init = vxyzu_init(1, clump(i)%partIDs(j))
!              part_vy_init = vxyzu_init(2, clump(i)%partIDs(j))
!              part_vz_init = vxyzu_init(3, clump(i)%partIDs(j))
!              part_v_init_tot = sqrt(part_vx_init**2 + part_vy_init**2 + part_vz_init**2)
!                  write(26,'(7Es18.6,1I8)') clump(i)%particle_positions(:, j), poten(clump(i)%partIDs(j))*unit_energ, &
!                  part_vx/part_vx_init, part_vy/part_vy_init, part_vz/part_vz_init, 0
!          enddo
!      endif
!  enddo
!  close(26)


! if (time_comp) then
!     write(filename,'(2a)') trim(dumpfile),'_clumps_timecomp'
! else
!     write(filename,'(2a)') trim(dumpfile),'_clumps'
! endif
! open (unit=26,file=trim(filename))
! write(26,'(2(a,I6),a,Es18.6,a,Es18.6,a,Es18.6)') '#Nclumps = ',iclump,'; Nclumpunique = ',&
!                            idclumpmax,'; Time = ',time,'; Background density = ',rhomin_bkg*(umass/udist**3),&
!                                    '; Mean density = ',rhomean*(umass/udist**3)! Time was added 21 Aug 2020
! do i = 1,iclump
!    write(26,'(3I8,16Es18.6)') clump(i)%ID,clump(i)%nmember,clump(i)%pointmass(1),                        & ! 1-3
!                            clump(i)%r(:),clump(i)%v(:),clump(i)%mass,clump(i)%size,cloudmax, & ! 4-6,7-9,10,11,12
!                            clump(i)%dr(:), & ! 13-15
!                            clump(i)%kinetic,clump(i)%potential,clump(i)%thermal,clump(i)%magnetic!, &   ! 16-18
!                            !rhomin_bkg*(umass/udist**3) ! 19
! enddo
! close(26)
!!
! do i = 1,iclump
!    if (.not. clump(i)%sink_only) then
!        if (time_comp) then
!             write(filename,'(2a,I5.5)') trim(prefix),'_clumpID_timecomp_bs_',clump(i)%ID
!        else
!             write(filename,'(2a,I5.5)') trim(prefix),'_clumpID_',clump(i)%ID
!        endif
!        inquire(file=trim(filename),exist=iexist)
!        if (iexist) then
!           open(26,file=trim(filename), position='append')
!        else
!           open(26,file=trim(filename))
!        endif
!        Nsync = count(clump(i)%pointmass /= 0)
!        if (Nsync > 0) then
!            sink_fmt = '(I4)'
!            write(sink_count, sink_fmt) Nsync - 1
!            sinkIDs(1:Nsync) = clump(i)%pointmass(1:Nsync)
!            print*, Nsync, sink_count, clump(i)%pointmass(2:Nsync) !, sinkIDs
!            write(26,'(Es18.6,5I8,18Es18.6,'//trim(sink_count)//'I4)') time,num,clump(i)%ID,clump(i)%nmember,clump(i)%pointmass(1),& ! 1-5
!                                             count(clump(i)%pointmass /= 0), clump(i)%r(:),clump(i)%v(:),clump(i)%mass, & ! 6,7-9,10-12,13
!                                             clump(i)%size,cloudmax,clump(i)%dr(:), & ! 14,15,16-18
!                                             clump(i)%kinetic,clump(i)%potential,clump(i)%thermal,clump(i)%magnetic, &    ! 19-22
!                                             rhomin_bkg*(umass/udist**3), rhomean*(umass/udist**3), &
!                                             sinkIDs(2:Nsync) ! sinkIDs - why doens't this one work??
!            do j = 1, Nsync
!                write(filename,'(1a,I5.5,1a,I4.4)') 'spare_sink_xyz_',clump(i)%ID,'_sink_',sinkIDs(j)
!                print*, trim(filename), 'bok', Nsync, j, sinkIDs(1), sinkIDs(2), sinkIDs(3), sinkIDs(j)
!                if (iexist) then
!                   open(27,file=trim(filename), position='append')
!                else
!                   open(27,file=trim(filename))
!                endif
!                write(27,'(4Es18.6)') time, clump(sinkIDs(j))%r
!            enddo
!            close(27)
!        else
!            write(26,'(Es18.6,5I8,18Es18.6)') time,num,clump(i)%ID,clump(i)%nmember,clump(i)%pointmass(1),& ! 1-5
!                                             count(clump(i)%pointmass /= 0), clump(i)%r(:),clump(i)%v(:),clump(i)%mass, & ! 6,7-9,10-12,13
!                                             clump(i)%size,cloudmax,clump(i)%dr(:), & ! 14,15,16-18
!                                             clump(i)%kinetic,clump(i)%potential,clump(i)%thermal,clump(i)%magnetic, &    ! 19-21
!                                             rhomin_bkg*(umass/udist**3), rhomean*(umass/udist**3)
!        endif
!        close(26)
!    endif
! enddo

!-----------------------------------------------------------------------------------------------------------------------------!
 
!----- Uncomment if you want to save comparison values between dumps -----!
! write(filename, '(2a)') trim(dumpfile),'rho_vs_Nclumps4'
! inquire(file = trim(filename), exist = iexist)
! if (iexist) then
!       open(26, file = trim(filename), position = 'append')
! else
!       open(26, file = trim(filename))
! endif
! write(26, '(I6, 3Es18.6)') idclumpmax, rhomin_bkg*(umass/udist**3), rhomax*(umass/udist**3), rhomean*(umass/udist**3)
! close(26)
!-------------------------------------------------------------------------!

  if (.false.) then
     do i = 1,nptmass
        write(368,*) 0,idclumpsink(i),xyzmh_ptmass(1:4,i)
        write(370,*) 0,idclumpsink(i),xyzmh_ptmass(1:4,i)
     enddo
     do i = 1,npart
        if (rho(i) > rhomin_bkg) write(368,*) i,idclump(i),xyzh(1:3,i),rho(i)*umass/udist**3
        if (idclump(i) > 0) write(369,*) i,idclump(i),xyzh(1:3,i),rho(i)*umass/udist**3
     enddo
  endif

 !--Save values for next dumpfile
 idclumpold = idclump
 idclumpsinkold = idclumpsink
 iorigold(1:npart) = iorig(1:npart)
 clumps_prev = clumps_cur
 firstcall = .false.

end subroutine do_analysis

!----------------------------------------------------------------------------------------------------------------------
! Now we want a SR for determining clump properties
subroutine obtain_membership(dumpfile, ii, xi, yi, zi, pmassi, iclump, npart, lst_x, lst_y, lst_z, idense, ndense, isink,&
            xyzh, vxyzu, rho, iis, time)

#ifdef PERIODIC
 use boundary, only: dxbound, dybound, dzbound
#endif
 use units, only: udist, umass, utime, unit_energ ! unit distance, unit mass

 ! parameters passed in/out of the SR
 integer, intent(inout) :: iclump ! clump index
 integer, intent(in) :: npart, isink, ndense, idense(:), lst_x(:), lst_y(:), lst_z(:), ii(3) ! isink is the sink index
 integer, intent(in) :: iis(3, maxptmass) ! Index list of ptmasses

 real, intent(in) :: xi, yi, zi, pmassi, xyzh(:,:), vxyzu(:,:), rho(:)
 real, intent(in) :: time

 character(len=*), intent(in) :: dumpfile
 character(len=128) :: filename_t
 
 ! parameters local to SR
 integer, parameter :: nmembermax = 500000 ! maximum number of particles which can be assigned to a clump
 integer :: i, j, jj, k, p, q, iiorig, imember, inotmember, nctr, ienc, imax, ival(2), ctr ! generally 'i' prefix => index value, jj is the index of the next closest particle, nctr and ctr are counter integers, and ival is a 2-value array, the first item is for storing a comparison value
 integer :: icandidate(npart), lst(npart), iprev(nclumpmax) ! icandidate stores array of candidate particle indices,
 integer :: part_enc ! Number of particles enclosed by the clump
 integer :: partID_array(npart), icandidate_xyz(3, npart), kx, ky, kz ! array storing particle clump IDs, candidate particle IDs in the xyz direction, and number of candidates in xyz

 logical :: keep_searching, full_list
 logical :: killed = .false., iny = .false., inz = .false. ! Logicals for checking clump life status, and whether a particle is in a given list
 
 real :: dist2tmp, sep, dv2, epot, etherm, ekin, dx, dy, dz, dr, menc, dxyz(3, nmembermax) ! temporary squared distance holder, separation, difference in velocities squared, energy values, separations, enclosed mass, storage of separations in 3D for member particles
 real :: dist2(4, npart), xyz(4, nmembermax) ! dist2 index 4 is for the final k list
 real :: dvxyz(3), axyz(3), a_tot, dtime , r_cen!, r_cen_abs ! Change in clump velocity between dumps, clump acceleration (vector and absolute), time difference between dumps
 real :: rho_clump, vol_clump, t_ff, clump_bulk_v ! pi, clump mean densities and volume, clump free-fall time and bulk velocity, dynamical timescale
 real :: venc, rho_cm(npart), t_dyn ! total volume enclosed by clump, current clump density, clump dynamical timescale
 real :: tff_cm(npart), m_cm(npart), clump_T(npart), clump_U(npart)! Running arrays for plotting how parameters evolve with clump mass e.g., volume, density, free-fall time, mass, kinetic and thermal energies
 real :: t_rat, e_rat ! Values storing time ratios, energy ratios, and the difference in these ratios
 real :: r_cm(3, npart), v_cm(3, npart), dr_cm(3, npart), size_cm(npart), clump_Th(npart), clump_Mag(npart) ! running totals of clump properties per particle added

 real :: mass_bins(5), part_loc(3, npart), ddist_tmp(3)!, U_part(npart) ! Array for sorting mass bins, 3D position array for each particle, 3D distance temporary holder, array for storing individual particle potentials

 ! initialise the values
 dist2 = 0.
 ddist_tmp = 0.
 k = 0
 icandidate = 0
 icandidate_xyz = 0
 etherm = 0.
 xyz = 0.
 xyz(1, 1) = xi
 xyz(2, 1) = yi
 xyz(3, 1) = zi
 xyz(4, 1) = pmassi
 iprev = 0
 imember = 0
 tff_cm = 0
 m_cm = 0
 rho_cm = 0
 clump_T = 0
 clump_U = 0
 r_cm = 0
 v_cm = 0
 dr_cm = 0
 size_cm = 0
 clump_Th = 0
 clump_Mag = 0
 !U_part = 0
 pi = 4*atan(1.0) ! Define pi
 partID_array = 0

 if (isink > 0) then
    xyz(4, 1) = xyzmh_ptmass(4, isink)
    if (idclumpsinkold(isink) > 0) iprev(idclumpsinkold(isink)) = 1 ! sink identified in a previous dump?
 endif

 ! Use the array of particles sorted in y to form a list of particles within r = cloudmax of the lead particles
 dx = 0.
 dy = 0.
 dz = 0.
 ddist_tmp = 0.
 jj = 1
 full_list = .false. ! we have not exxhausted the entire particle list
 if (isink > 0) jj = 0 ! sink is lead
 do q = 1, 3
     ctr = 0
     do p = 1, 2
        !print*, 'here too??'
        do while(abs(ddist_tmp(q)) < cloudmax .and. ctr < ndense)
            ctr = ctr + 1 ! failsafe for loop exit
            if (full_list) then
                exit ! Have checked all available particles
            endif
            ! Move out in one direction
            if (p == 1) then
                ival(2) = ii(q) + jj ! ii is the position in the lst_rho
                if (ival(2) > ndense) then
                    ival(2) = ival(2) - ndense
                    if (.not. periodic) ctr = ndense + 1 ! exit loop
                endif
            ! Move out in opposite direction
            else
                ival(2) = ii(q) - jj ! We are moving up and down lst_rho
                if (ival(2) < 1) then
                    ival(2) = ival(2) + ndense
                    if (.not. periodic) ctr = ndense + 1 ! exit loop
                endif
            endif
            if (ctr == 1) then
                ival(1) = ival(2)
            else
                if (ival(2) == ival(1) - 1) full_list = .true.
            endif
            if (q == 1) then
                j = idense(lst_x(ival(2))) ! find original index of particle
            else if (q == 2) then
                j = idense(lst_y(ival(2)))
            else if (q == 3) then
                j = idense(lst_z(ival(2)))
            endif
            if (idclump(j) == 0) then ! j not assigned to a clump
                ! separation of lead particle and current particle
                dx = xi - xyzh(1, j)
                dy = yi - xyzh(2, j)
                dz = zi - xyzh(3, j)

                if (q == 1) then
                    ddist_tmp(q) = dx
                else if (q == 2) then
                    ddist_tmp(q) = dy
                else if (q == 3) then
                    ddist_tmp(q) = dz
                endif

#ifdef PERIODIC
                ! account for particles on the edge of the box
                if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0, dx)
                if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0, dy)
                if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0, dz)
#endif
                dist2tmp = dx*dx + dy*dy + dz*dz
                if (dist2tmp < cloudmax2) then
                    k = k + 1 ! counter for number of clump candidates
                    dist2(q, k) = dist2tmp ! save position of candidate particle ! current issue
                    icandidate_xyz(q, k) = j ! save index of candidate particle
                endif
            endif
            jj = jj + 1 ! Move along
        enddo
        ! reset properties for search in the opposite direction
        ddist_tmp = 0.
        jj = 1
     enddo
 if (q == 1) then
     kx = k ! Save the maximum number of particles checked by x
 else if (q == 2) then
     ky = k
 else if (q == 3) then
     kz = k
 endif
 ddist_tmp = 0.
 jj = 1
 k = 0
 full_list = .false.
 enddo
 k = 0
 do i = 1, kx ! Don't need to do for all 3 as if they must appear in all 3
     ! check if the x index candidate appears in the other two arrays
     do j = 1, ky
         if (icandidate_xyz(2, j) == icandidate_xyz(1, i)) then
             iny = .true.
         endif
     enddo
     if (iny .eqv. .true.) then
         do j = 1, kz
             if (icandidate_xyz(3, j) == icandidate_xyz(1, i)) then
                 inz = .true.
                 k = k + 1
                 dist2(4, i) = dist2(1, i)
                 icandidate(i) = icandidate_xyz(1, i)
             endif
         enddo
     endif
     iny = .false. ! Reset for next kx
     inz = .false.

 enddo

 if (k > 0) then ! we have clump candidates
    call indexx(k, dist2(4, :), lst) ! sort them in order of distance from the clump leader
    keep_searching = .true.
    inotmember = 0 ! for storing the number of rejected particles
    j = 1
    nctr = 1

    do while (keep_searching)
        i = icandidate(lst(j)) ! original particle index
        dv2 = (clump(iclump)%v(1) - vxyzu(1, i))**2 + (clump(iclump)%v(2) - vxyzu(2, i))**2 + (clump(iclump)%v(3) - vxyzu(3, i))**2 ! particle velocity relative to clump
        nctr = nctr + 1
        xyz(1, nctr) = xyzh(1, i)
        xyz(2, nctr) = xyzh(2, i)
        xyz(3, nctr) = xyzh(3, i)
        xyz(4, nctr) = pmassi

        epot = 0.
        sep = 0.
        iiorig = minloc(abs(iorigold(1:npart) - iorig(i)), 1) ! value of this should just be zero

        ! particle contibution to the potential energy
        ! n**2 operation => very slow
        do p = 1, nctr - 1
            dx = xyz(1, p) - xyz(1, nctr)
            dy = xyz(2, p) - xyz(2, nctr)
            dz = xyz(3, p) - xyz(3, nctr)
            dr = sqrt(dx*dx + dy*dy + dz*dz)
            sep = max(sep, dr)
            if (dr > 0.) epot = epot + xyz(4, p)/dr ! internal potential
        enddo

        ! find other energy contributions
        ekin = clump(iclump)%kinetic + 0.5*pmassi*dv2
        epot = clump(iclump)%potential - epot*pmassi
        !U_part(nctr) = epot ! Save individual particle potential
        if (maxvxyzu >= 4) etherm = clump(iclump)%thermal + vxyzu(4, i)*pmassi

        ! check particle boundness - if bound it is added to the clump
        if (ekin_coef*ekin + epot < 0. .or. clump(iclump)%nmember < nmin) then

            ! second condition is where clump is not yet big enough
            idclump(i) = idtmp
            if (idclumpold(iiorig) > 0) iprev(idclumpold(iiorig)) = iprev(idclumpold(iiorig)) + 1 ! Move up a clump index

            ! centre of mass/velocity
            do p = 1, 3
                clump(iclump)%r(p) = (clump(iclump)%mass*clump(iclump)%r(p) + pmassi*xyzh(p, i))/(clump(iclump)%mass + pmassi)
                clump(iclump)%v(p) = (clump(iclump)%mass*clump(iclump)%v(p) + pmassi*vxyzu(p, i))/(clump(iclump)%mass + pmassi)
                r_cm(p, nctr) = (clump(iclump)%mass*clump(iclump)%r(p) + pmassi*xyzh(p, i))/(clump(iclump)%mass + pmassi)
                v_cm(p, nctr) = (clump(iclump)%mass*clump(iclump)%v(p) + pmassi*vxyzu(p, i))/(clump(iclump)%mass + pmassi) ! Save running totals for clump COM and COV
            enddo
            
            ! update clump properties
            clump(iclump)%nmember = nctr
            clump(iclump)%kinetic = ekin
            clump(iclump)%potential = epot
            clump(iclump)%thermal = etherm

            ! Save running total values
            clump_T(nctr) = ekin
            clump_U(nctr) = epot
            clump_Th(nctr) = etherm
            

            if (mhd) then
                clump(iclump)%magnetic = clump(iclump)%magnetic + 0.5*pmassi*dot_product(Bxyz(1:3, i), Bxyz(1:3, i))/rho(i)
                clump_Mag(nctr) = clump(iclump)%magnetic + 0.5*pmassi*dot_product(Bxyz(1:3, i), Bxyz(1:3, i))/rho(i)
            endif
            clump(iclump)%mass = clump(iclump)%mass + pmassi
            clump(iclump)%size = max(clump(iclump)%size, 0.5*sep)
            size_cm(nctr) = max(clump(iclump)%size, 0.5*sep)
            venc = (4/3)*pi*((clump(iclump)%size))**3 ! Assume nctr is just the index of clump particles (as it is reset below if a particle is not added)
            m_cm(nctr) = (0.9*nctr*pmassi)
            rho_cm(nctr) = m_cm(nctr)/venc
            tff_cm(nctr) = (3*pi/(32*rho_cm(nctr)))**0.5 ! take out make loop
            partID_array(nctr) = i ! Particle is assigned to the clump
        else
            ! candidate is not a clump member
            inotmember = inotmember + 1
            nctr = nctr - 1
        endif

        ! stop clump identification if non-members outweigh members or we have exhausted our list
        if (inotmember >= clump(iclump)%nmember .or. j == k) keep_searching = .false.
        j = j + 1
    enddo
    ekin = clump(iclump)%kinetic
    epot = clump(iclump)%potential
    etherm = clump(iclump)%thermal


    if (ekin_coef*ekin + epot < 0. .and. ekin > 0.) then
    ! Second requirement prevents clumps from being identified in initially unmmoving background
    ! Check if sink was identified in a previous dump, use its imember value for clump id
        if (isink > 0) then
            if (idclumpsinkold(isink) > 0) then
                imember = idclumpsinkold(isink)
                idclumpsink(isink) = imember
            endif
        endif
    
    ! Same but for sinkless clumps: check which id from previous dumps was most used (if it was not already used) => most number of particles
        if (imember == 0) then
            imax = maxval(iprev, 1)
            if (imax > 0) then
                imember = maxloc(iprev, 1)
                do jj = 1, iclump - 1
                    if (clump(jj)%ID == imember) imember = 0 ! If the cloud is split we don't want to duplicate ID numbers
                enddo
            endif
            if (imember == 0) then
                idclumpmax = idclumpmax + 1
                imember = idclumpmax
            endif
        endif
        do jj = 1, k
            i = icandidate(lst(jj)) ! index of candidate particle from original list
            if (idclump(i) == idtmp) idclump(i) = imember
        enddo
        clump(iclump)%ID = imember

        ! find radius of clump which contains 90% of its mass, and update the number of members it contains
        dist2 = 0.
        dxyz = 0.
        clump(iclump)%dr(:) = 0.
        do j = 1, nctr
            dx = clump(iclump)%r(1) - xyz(1, j)
            dy = clump(iclump)%r(2) - xyz(2, j)
            dz = clump(iclump)%r(3) - xyz(3, j)
#ifdef PERIODIC
            if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0, dx)
            if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0, dy)
            if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0, dz)
#endif
            dxyz(1, j) = dx*dx
            dxyz(2, j) = dy*dy
            dxyz(3, j) = dz*dz
            dist2(4, j)  = dxyz(1, j) + dxyz(2, j) + dxyz(3, j) ! check this doesn't come back to bite me in the ass
        
            do p = 1, 3
                clump(iclump)%dr(p) = max(clump(iclump)%dr(p), dxyz(p, j))
                dr_cm(p, j) = max(clump(iclump)%dr(p), dxyz(p, j))
            enddo
        enddo
        
        if (isink > 0) then
            ! need to account for sink mass
            menc = 0.9*(clump(iclump)%mass - xyzmh_ptmass(4, isink))
        else
            ! just gas contribution to mass
            menc = 0.9*(clump(iclump)%mass)
        endif

        ienc = int(menc/pmassi + 1.0) ! number of particles enclosed
        call indexx(nctr, dist2(4, :), lst) ! order particles enclosed by distance from centre
        clump(iclump)%size = sqrt(dist2(4, lst(ienc)))

        ! Want to get the clump free-fall time
        ! Going to assume a constant density throughout, and that the clump is spherical => get density by dividing mass by volume
        print*, 'Clump member number before timecomp = ', clump(iclump)%nmember
        if (time_comp) then
            vol_clump = (4/3)*pi*(clump(iclump)%size)**3 ! Determine clump volume
            rho_clump = (menc)/venc
            t_ff = (3*pi/(32*rho_clump))**0.5 ! is G = 1 in code units - yes
            tff_cm(nctr) = t_ff ! Reset this final value to extend to full clump size
            ! Clump should be killed if t_ff > t_dyn
            dtime = abs(time - time_prev)
            if (isink > 0) then
                if (isink > nptmass_prev) then
                    dvxyz(1) = vxyz_ptmass(1, isink) - vxyzu_prev(1, idense(lst_y(iis(1, isink))))
                    dvxyz(2) = vxyz_ptmass(2, isink) - vxyzu_prev(2, idense(lst_y(iis(1, isink))))
                    dvxyz(3) = vxyz_ptmass(3, isink) - vxyzu_prev(2, idense(lst_y(iis(1, isink)))) ! Taking velocity of particle closest to current sink in previous dump if sink did not exist in the previous dump
                    axyz(:) = dvxyz(:)/dtime
                else
                    dvxyz(1) = vxyz_ptmass(1, isink) - vxyz_ptmass_prev(1, isink)
                    dvxyz(2) = vxyz_ptmass(2, isink) - vxyz_ptmass_prev(2, isink)
                    dvxyz(3) = vxyz_ptmass(3, isink) - vxyz_ptmass_prev(3, isink) ! Sink acceleration defines clump acceleration
                    axyz(:) = dvxyz(:)/dtime
                endif
            else
                dvxyz(1) = vxyzu(1, ii(1)) - vxyzu_prev(1, ii(1))
                dvxyz(2) = vxyzu(2, ii(2)) - vxyzu_prev(2, ii(2))
                dvxyz(3) = vxyzu(3, ii(3)) - vxyzu_prev(3, ii(3)) ! If no sink we just use the acceleration of the lead particle in the clump
                axyz(:) = dvxyz(:)/dtime
            endif
            clump(iclump)%a = axyz ! Store clump acceleration
            a_tot = sqrt((axyz(1)*axyz(1)) + (axyz(2)*axyz(2)) + (axyz(3)*axyz(3))) ! absolute value of clump acceleration
            clump_bulk_v = sqrt((clump(iclump)%v(1))**2 + (clump(iclump)%v(2))**2 + (clump(iclump)%v(3))**2)

            r_cen = (clump_bulk_v**2)/a_tot

            t_dyn = clump_bulk_v/a_tot

            clump(iclump)%time_ff = t_ff
            clump(iclump)%time_dyn = t_dyn

        endif
                
        do p = 1, 3
            clump(iclump)%dr(p) = sqrt(clump(iclump)%dr(p))
        enddo
    endif
    print*, 'eligible particles near candidate particle found, ', 'number of candidates =', k, 'host =', ii
 else
    print*, 'no eligible particles near the candidate particle'
 endif

! If we don't have a bound clump we want to kill it
 if (imember == 0) then
    killed = .true.
    do jj = 1, k
        i = icandidate(lst(jj))
        if (idclump(i) == idtmp) then
            if (jj < nmin) then
                idclump(i) = -1 ! not enough members
            else
                idclump(i) = 0
            endif
        endif
    enddo

    ! resetting all properties

    call reset_clump(iclump)
   
    if (isink > 0) then
        print*, 'ruhoh there was a sinky winky in the clumpy wumpy'
        clump(iclump)%r(1:3) = xyzmh_ptmass(1:3, isink)
        clump(iclump)%v(1:3) = vxyz_ptmass(1:3, isink)
        clump(iclump)%mass = xyzmh_ptmass(4, isink)
        clump(iclump)%size = xyzmh_ptmass(ihacc, isink)
        clump(iclump)%pointmass(1) = isink
        clump(iclump)%sink_only = .true.
       
        ! same bookkeeping as before (for multiple dumpfile analysis)
        if (idclumpsinkold(isink) > 0) then
            clump(iclump)%ID = idclumpsinkold(isink)
        else
            idclumpmax = idclumpmax + 1
            clump(iclump)%ID = idclumpmax
        endif
        idclumpsink(isink) = clump(iclump)%ID
    else
        iclump = iclump - 1
    endif
 endif
 if (.not. killed) then
    allocate(clump(iclump)%partIDs(nctr))
    clump(iclump)%partIDs = partID_array

    do i = 1, nctr
            e_rat = abs(clump_U(i)/clump_T(i))
            t_rat = tff_cm(i)/t_dyn
            if (e_rat < t_rat) then ! Move up i until this condition is met
                print*, 'Resetting clump based on free-fall timescale', clump(iclump)%nmember, i, e_rat, t_rat
                clump(iclump)%nmember = i
                nctr = i
                clump(iclump)%r(1:3) = r_cm(1:3, i)
                clump(iclump)%v(1:3) = v_cm(1:3, i)
                clump(iclump)%dr(1:3) = dr_cm(1:3, i)
                clump(iclump)%mass = m_cm(i)
                clump(iclump)%size = size_cm(i)
                clump(iclump)%kinetic = clump_T(i)
                clump(iclump)%potential = clump_U(i)
                clump(iclump)%thermal = clump_Th(i)
                clump(iclump)%magnetic = clump_Mag(i)
                clump(iclump)%time_ff = tff_cm(i)
                clump(iclump)%partIDs = partID_array(1:i)
                exit
            endif
        enddo
        allocate(clump(iclump)%particle_positions(dims, clump(iclump)%nmember))
        clump(iclump)%particle_positions = xyz(1:3, 1:clump(iclump)%nmember)

!---------------------------------------------------------- USEFUL WRITING TO FILE ----------------------------------------------------------!
! IF YOU WANT MASS BINS COMMENT THIS BIT IN -------------------------------------------------------------------------------------------------!
!     do i = 1, 5
!         mass_bins(i) = ((clump(iclump)%mass)/5)*i
!         part_enc = ((clump(iclump)%nmember)/5)*i
!         do j = 1, part_enc
!             part_loc(1:3, j) = xyz(1:3, j)
!             print*, part_loc(1:3, j), xyz(1:3, j)
!             read(*,*)
!         enddo
!         write(filename_t, '(2a,I1,a,I5.5)') trim(dumpfile),'_part_pos_mass_bin_',i,'_',clump(iclump)%ID
!         open(unit = 26, file = trim(filename_t))
!         write(26, '(a,I8,a,Es18.6)') '# Mass enclosed by ', part_enc, ' particles = ', mass_bins(i)
!         do j = 1, part_enc
!             write(26, '(5Es18.6)') part_loc(1:3, j), rho(partID_array(j))*(umass/udist**3), U_part(j)
!         enddo
!     enddo
!     close(26)
!---------------------------------------------------------------------------------------------------------------------------------------------!
! IF YOU WANT SINGLE CLUMP DATA COMMENT THIS BIT IN ------------------------------------------------------------------------------------------!
!     write(filename_t, '(2a,I1,a,I5.5)') trim(dumpfile),'_part_pot_minima_',clump(iclump)%ID
!     open(unit = 26, file = trim(filename_t))
!     do j = 1, (clump(iclump)%nmember)
!         print*, xyz(1:3, j)
!         read(*,*)
         !write(26, '(5Es18.6)') xyz(1:3, j), rho(partID_array(j))*(umass/udist**3), U_part(j)*unit_energ
!     enddo
!     close(26)
!---------------------------------------------------------------------------------------------------------------------------------------------!
!---------------------------------------------------------------------------------------------------------------------------------------------!
 endif
 killed = .false.

end subroutine obtain_membership

!----------------------------------------------------------------------------------------------------------------------

! Subroutine for resetting clump parameters
subroutine reset_clump(i)
 integer, intent(in) :: i

 clump(i)%r(:) = 0.
 clump(i)%v(:) = 0.
 clump(i)%dr(:) = 0.
 clump(i)%a(:) = 0.
 clump(i)%ID = 0
 clump(i)%mass = 0.
 clump(i)%size = 0.
 clump(i)%nmember = 0
 clump(i)%kinetic = 0.
 clump(i)%potential = 0.
 clump(i)%thermal = 0.
 clump(i)%magnetic = 0.
 clump(i)%pointmass = 0
 clump(i)%time_dynv = 0.
 clump(i)%time_dyna = 0.
 clump(i)%conv_ID_v = 0.
 clump(i)%conv_ID_a = 0.
 clump(i)%SURV = .false.
 clump(i)%NEW = .true.
 if (allocated(clump(i)%partIDs)) then
     deallocate(clump(i)%partIDs)
 endif
 if (allocated(clump(i)%particle_positions)) then
     deallocate(clump(i)%particle_positions)
 endif

end subroutine reset_clump

!----------------------------------------------------------------------------------------------------------------------

end module


