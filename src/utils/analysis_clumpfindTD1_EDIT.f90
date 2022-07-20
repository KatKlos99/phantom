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
 use dim,        only:maxptmass,maxvxyzu,mhd,maxp_hard
 use part,       only:Bxyz,xyzmh_ptmass,vxyz_ptmass,nptmass,ihacc,periodic,iorig
 use sortutils,  only:indexx
 implicit none
 character(len=20), parameter, public :: analysistype = 'clumpfindTD'
 public  :: do_analysis
 integer, parameter :: nclumpmax = 1000 ! the maximum number of unique clumps
 integer :: idclumpold(maxp_hard),idclump(maxp_hard),iorigold(maxp_hard)
 integer :: idclumpsinkold(maxptmass),idclumpsink(maxptmass)
 integer :: idclumpmax                  ! the running total of clumps
 integer :: nmin           = 57         ! minimum number of particles to define a pure-gas clump
 integer :: idtmp          = -10        ! a temporary id 
 logical :: firstcall      = .true.     ! required logical; do not change
 real    :: rhominpeak_cgs = 1.0d-20    ! particles with rho < rhominpeak will not be considered for the lead particle in a clumps
 real    :: rhominbkg_cgs  = 1.0d-20    ! particles with rho < rhominbkg  will not be considered for membership in a clumps
 real    :: cloudmax_pc    = 20.        ! maximum possible radius of the cloud
 real    :: cloudmax2,cloudmax
 real    :: ekin_coef      = 2.         ! Will be a clump if ekin_coef*ekin + epot < 0; ekin_coef can be read in
 real    :: yexc           = 0.         ! Will exclude gas particles with -yexc < y < yexc

 type sphclump
    integer           :: nmember,pointmass,ID
    real,dimension(3) :: r,v,dr
    real              :: mass,size
    real              :: kinetic,potential,thermal,magnetic,bound
 end type sphclump
 type(sphclump)       :: clump(nclumpmax)

 private

contains
!--------------------------------------------------------------------------
subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,       only:rhoh
 use physcon,    only:pc
 use units,      only:udist,umass
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer                      :: ii,i,j,k,io,iclump,iclumptmp,idx,izero,ndense, COUNTER
 integer                      :: lst_rho(npart),lst_y(npart),idense(npart),iis(maxptmass)
 real                         :: rhomin_peak,rhomin_bkg,xi,yi,zi,hi,dy,rtmp
 real                         :: rho(npart),rho_dense(npart),rpos(npart),dymin(maxptmass)
 logical                      :: keep_searching,iexist
 character(len=128)           :: filename,prefix

 !--Initialise values
 izero       = 0
 iclump      = 0
 rhomin_peak = rhominpeak_cgs*udist**3/umass
 rhomin_bkg  = rhominbkg_cgs*udist**3/umass
 cloudmax    = cloudmax_pc*pc/udist
 cloudmax2   = cloudmax**2
 if (firstcall) then
    idclumpold     = 0
    idclumpsinkold = 0
    idclumpmax     = 0
    iorigold(1:npart)       = iorig(1:npart)
    inquire(file='clumpfind_KEcoef.in',exist=iexist)
    if (iexist) then
       open(unit=122,file='clumpfind_KEcoef.in')
       read(122,*,iostat=io) rtmp
       if (io==0) ekin_coef = rtmp
       close(122)
    endif
    inquire(file='clumpfind_yexclude.in',exist=iexist)
    if (iexist) then
       open(unit=122,file='clumpfind_yexclude.in')
       read(122,*,iostat=io) rtmp
       if (io==0) yexc = abs(rtmp)
       close(122)
    endif
 endif
 idclump     = 0
 idclumpsink = 0
 idx         = index(dumpfile,'_')
 prefix      = dumpfile(1:idx-1)
 print*, 'Groups of particles are considered a clump if ',ekin_coef,'*ekin + epot < 0'
 if (yexc > 0.) then
    print*, 'Excluding particles with -',yexc,' < y/code < ',yexc
 endif

 !--calculate all densities & y-position
 rho    = 0.0
 rpos   = 0.0
 idense = 0
 ndense = 0
 do i = 1,npart
    yi = xyzh(2,i)
    hi = xyzh(4,i)
    if (hi > 0.0) then
       rho(i) = rhoh(hi,particlemass)
       if (rho(i) > rhomin_bkg .and. (yi >= yexc .or. yi <= -yexc)) then
          ndense            = ndense + 1
          idense(ndense)    = i
          rho_dense(ndense) = rho(i)
          rpos(ndense)      = xyzh(2,i)    ! This could be any of the three position; choose the best one for your structure (but need to modify other location that use lst_y)
       endif
    endif
 enddo
 print*, 'Dense gas located', ndense, npart ! have definitely found dense gas

! do i = 1,npart
!    if (idense(i) > npart) print*, idense(i)
! end do
! read(*,*)

 !--Sort the dense particles by position and density
 call indexx(ndense,rho_dense,lst_rho)
 call indexx(ndense,rpos,lst_y)
 print*, 'array sorting complete'

 !--Find the gas particle nearest to the sink in the y-direction
 iis(:)   = 0
 dymin = huge(dymin)
 do j = 1,ndense
    do i = 1,nptmass
       dy = abs(xyzmh_ptmass(2,i) - xyzh(2,idense(lst_y(j))))
       if (dy < dymin(i)) then
          iis(i)   = j
          dymin(i) = dy
       endif
    enddo
 enddo
 print*, 'found gas nearest to each sink'

 print*, SIZE(iis), nptmass
 !read(*,*)

 !--Find the clumps around sinks
 if (nptmass > 0) then
     do i = 1,nptmass
        iclump = iclump + 1
        xi = xyzmh_ptmass(1,i)
        yi = xyzmh_ptmass(2,i)
        zi = xyzmh_ptmass(3,i)
        ! initialise clump
        call reset_clump(iclump)
        clump(iclump)%r(1:3)    = xyzmh_ptmass(1:3,i)
        clump(iclump)%v(1:3)    = vxyz_ptmass(1:3,i)
        clump(iclump)%mass      = xyzmh_ptmass(4,i)
        clump(iclump)%size      = xyzmh_ptmass(ihacc,i)
        clump(iclump)%pointmass = i
        ! find members
        call obtain_membership(iis(i),xi,yi,zi,particlemass,iclump,npart,lst_y,idense,ndense,i,xyzh,vxyzu,rho)
        print*, 'CLUMP_SINK!',iclump,clump(iclump)%ID,xyzmh_ptmass(4,i),clump(iclump)%nmember, &
                 clump(iclump)%kinetic,clump(iclump)%thermal,clump(iclump)%potential,clump(iclump)%mass,clump(iclump)%size
     enddo
     print*, 'found clumps around sinks'
 else
     print*, 'no sinks in dump'
 end if

 !--Find the clumps around dense gas particles
 k              = ndense
 keep_searching = .true.
 do while(keep_searching)
    do while (idclump(idense(lst_rho(k)))/=0 .and. k > 2)
       k = k - 1
    enddo
    i=idense(lst_rho(k))  ! the densest particle
    if (rho(i) > rhomin_peak) then
       iclump    = iclump + 1
       iclumptmp = iclump
       xi = xyzh(1,i)
       yi = xyzh(2,i)
       zi = xyzh(3,i)
       ! initialise clump
       call reset_clump(iclump)
       clump(iclump)%r(1:3)  = xyzh(1:3,i)
       clump(iclump)%v(1:3)  = vxyzu(1:3,i)
       clump(iclump)%nmember = 1
       clump(iclump)%mass    = particlemass
       if (maxvxyzu>=4) clump(iclump)%thermal = vxyzu(4,i)*particlemass
       if (mhd)         clump(iclump)%magnetic= 0.5*particlemass*dot_product(Bxyz(1:3,i),Bxyz(1:3,i))/rho(i)
       idclump(i) = idtmp
       ! find members
!      ii = findloc(idense(lst_y(1:ndense)),i,1)          ! this inline function is not available on all compilers
       ii = minloc(abs(idense(lst_y(1:ndense)) - i), 1)
       print*, ii, i, minval(abs(idense(lst_y(1:ndense)) - i), 1)
       !read(*,*)
       !print*, ii
       call obtain_membership(ii,xi,yi,zi,particlemass,iclump,npart,lst_y,idense,ndense,izero,xyzh,vxyzu,rho)
       !print*, 'bonk!'
       ! print clump details, or remove peak clump candidate from list
       if (iclump < iclumptmp) then
          idclump(i) = -1
       else
          idclump(i) = clump(iclump)%ID
          print*, 'CLUMP!     ',iclump,clump(iclump)%ID,rho(i)*umass/udist**3,clump(iclump)%nmember, &
                  clump(iclump)%kinetic,clump(iclump)%thermal,clump(iclump)%potential,clump(iclump)%mass,clump(iclump)%size
       endif
    else
       keep_searching = .false.
    endif
 enddo

 !--Write results to file (both all the clumps at the current time, and to the file for each clump)


write(filename, '(1a)') 'iclump_per_dump' !'reassigned_per_dump'
if (iclump > 0) then
do i = 1, iclump
!     print*, clump(i)%ID
    !if (.not. clump(i)%reassigned) then
        !if (clump(i)%reassigned) then
                
            inquire(file=trim(filename),exist=iexist)
            if (iexist) then
                open(26,file=trim(filename), position='append')
            else
                open(26,file=trim(filename))
            endif
            
            write(26, '(1I8)') iclump!clump(i)%reassigned
        !endif
enddo
endif
close(26)
inquire(file=trim(filename),exist=iexist)
if (iexist) open(26,file=trim(filename), position='append')
write(26, '(1a)') dumpfile
close(26)


! write(filename,'(2a)') trim(dumpfile),'clumps'
! open (unit=26,file=trim(filename))
! write(26,'(2(a,I6),a,Es18.6)') '#Nclumps = ',iclump,'; Nclumpunique = ',idclumpmax,'; Time = ',time     ! Time was added 21 Aug 2020
! do i = 1,iclump
!    write(26,'(3I8,15Es18.6)') clump(i)%ID,clump(i)%nmember,clump(i)%pointmass,                        & ! 1-3
!                               clump(i)%r(:),clump(i)%v(:),clump(i)%mass,clump(i)%size,clump(i)%dr(:), & ! 4-6,7-9,10,11,12-14
!                               clump(i)%kinetic,clump(i)%potential,clump(i)%thermal,clump(i)%magnetic    ! 15-17
! enddo
! close(26)
! do i = 1,iclump
!    write(filename,'(2a,I5.5)') trim(prefix),'clumpID',clump(i)%ID
!    inquire(file=trim(filename),exist=iexist)
!    if (iexist) then
!       open(26,file=trim(filename), position='append')
!    else
!       open(26,file=trim(filename))
!    endif
!    write(26,'(Es18.6,4I8,15Es18.6)') time,num,clump(i)%ID,clump(i)%nmember,clump(i)%pointmass,               & ! 1-5
!                                      clump(i)%r(:),clump(i)%v(:),clump(i)%mass,clump(i)%size,clump(i)%dr(:), & ! 6-8,9-11,12,13,14-16
!                                      clump(i)%kinetic,clump(i)%potential,clump(i)%thermal,clump(i)%magnetic    ! 17-19
!    close(26)
! enddo

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
 idclumpold     = idclump
 idclumpsinkold = idclumpsink
 iorigold       = iorig
 firstcall      = .false.

end subroutine do_analysis
!--------------------------------------------------------------------------
! This is the subroutine that will determine the clump properties
subroutine obtain_membership(ii,xi,yi,zi,pmassi,iclump,npart,lst_y,idense,ndense,isink,xyzh,vxyzu,rho)
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
 integer, intent(inout) :: iclump
 integer, intent(in)    :: npart,isink,ndense,idense(:),lst_y(:), ii
 real,    intent(in)    :: xi,yi,zi,pmassi,xyzh(:,:),vxyzu(:,:),rho(:)
 integer, parameter     :: nmembermax = 500000
 integer                :: i,j,jj,k,p,iiorig,imember,inotmember,nctr,ienc,imax,ival,ctr, COUNTER, IVAL_ST
 integer                :: icandidate(npart),lst(npart),iprev(nclumpmax)
 real                   :: dist2tmp,sep,dv2,epot,etherm,ekin,dx,dy,dz,dr,menc,dxyz(3,nmembermax)
 real                   :: dist2(npart),xyz(4,nmembermax)
 logical                :: keep_searching, full_list

 !--Initialise values
 dist2      = 0.
 k          = 0
 icandidate = 0
 etherm     = 0.
 xyz        = 0.
 xyz(1,1)   = xi
 xyz(2,1)   = yi
 xyz(3,1)   = zi
 xyz(4,1)   = pmassi
 iprev      = 0
 imember    = 0
 if (isink > 0) then
    xyz(4,1) = xyzmh_ptmass(4,isink)
    if (idclumpsinkold(isink) > 0) iprev(idclumpsinkold(isink)) = 1
 endif
 COUNTER = 0

 !--Create list of particles within r=cloudmax of the peak
 !  this is done by looking through the array that is sorted in y
 dy = 0.
 jj = 1
 full_list = .false. ! To account for whether we have checked all available particles
 if (isink > 0) jj = 0 ! to account for the nearest particle
 do p = 1,2
    ctr = 0
    do while(abs(dy) < cloudmax .and. ctr <= ndense)
       ctr = ctr + 1 ! failsafe
       if (full_list) exit ! We do not want to exceed bounds of npart
       if (p==1) then
          ival = ii + jj
          if (ival > ndense) then
             ival = ival - ndense
             if (.not. periodic) ctr  = ndense + 1 ! stop the loop
          endif
       else
          ival = ii - jj
          if (ival < 1) then
             ival = ival + ndense
             if (.not. periodic) ctr  = ndense + 1 ! stop the loop
          endif
       endif

       if (ctr == 1) then
          IVAL_ST = ival ! save the first ival here
       else
          if (ival == IVAL_ST-1) full_list = .true.
       end if

       j = idense(lst_y(ival))
       if (j < 0) print*, j
       if (idclump(j)==0) then
          dx = xi-xyzh(1,j)
          dy = yi-xyzh(2,j)
          dz = zi-xyzh(3,j)
#ifdef PERIODIC
          if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
          if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
          if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
          dist2tmp = dx*dx + dy*dy + dz*dz
          if (dist2tmp < cloudmax2) then
             !print*, k, p, ival, j
             k = k + 1
             dist2(k) = dist2tmp
             icandidate(k) = j
             !if ((p == 2) .and. (k > npart)) print*, icandidate(k), SIZE(icandidate)/npart, p
             !read(*,*)
          endif
       endif
       jj = jj + 1
    enddo
    ! now search in the opposite direction
    dy = 0.
    jj = 1
 enddo

 !print*, npart

 !--Determine membership & clump properties
 if (k > 0) then
    call indexx(k,dist2,lst)
    keep_searching = .true.
    inotmember     = 0
    j              = 1
    nctr           = 1
    do while (keep_searching)
       i = icandidate(lst(j))
       !if (lst(j) > npart) print*, lst(j)
       ! properties of the candidate particle
       dv2 = (clump(iclump)%v(1) - vxyzu(1,i))**2 + (clump(iclump)%v(2) - vxyzu(2,i))**2 + (clump(iclump)%v(3) - vxyzu(3,i))**2
       nctr        = nctr + 1
       xyz(1,nctr) = xyzh(1,i)
       xyz(2,nctr) = xyzh(2,i)
       xyz(3,nctr) = xyzh(3,i)
       xyz(4,nctr) = pmassi
       epot        = 0.
       sep         = 0.
       iiorig      = minloc(abs(iorigold(1:npart) - iorig(i)), 1) ! the current particle's position in the array in the previous step
       do p = 1,nctr-1
          dx   = xyz(1,p) - xyz(1,nctr)
          dy   = xyz(2,p) - xyz(2,nctr)
          dz   = xyz(3,p) - xyz(3,nctr)
          dr   = sqrt(dx*dx + dy*dy + dz*dz)
          sep  = max(sep,dr)
          if (dr > 0.) epot = epot + xyz(4,p)/dr
       enddo
       ekin = clump(iclump)%kinetic   + 0.5*pmassi*dv2
       epot = clump(iclump)%potential - epot*pmassi
       if (maxvxyzu>=4) etherm = clump(iclump)%thermal + vxyzu(4,i)*pmassi

       ! if particle is bound, add contribution to clump
       if (ekin_coef*ekin + epot < 0. .or. clump(iclump)%nmember < nmin) then
          idclump(i) = idtmp
          if (idclumpold(iiorig) > 0) iprev(idclumpold(iiorig)) = iprev(idclumpold(iiorig)) + 1
          do p=1,3
             clump(iclump)%r(p)  = (clump(iclump)%mass*clump(iclump)%r(p) + pmassi*xyzh( p,i))/(clump(iclump)%mass + pmassi)
             clump(iclump)%v(p)  = (clump(iclump)%mass*clump(iclump)%v(p) + pmassi*vxyzu(p,i))/(clump(iclump)%mass + pmassi)
          enddo
          clump(iclump)%nmember   = nctr
          clump(iclump)%kinetic   = ekin
          clump(iclump)%potential = epot
          clump(iclump)%thermal   = etherm
          if (mhd) then
             clump(iclump)%magnetic = clump(iclump)%magnetic + 0.5*pmassi*dot_product(Bxyz(1:3,i),Bxyz(1:3,i))/rho(i)
          endif
          clump(iclump)%mass      = clump(iclump)%mass + pmassi
          clump(iclump)%size      = max(clump(iclump)%size,0.5*sep)
       else
          inotmember = inotmember + 1
          nctr       = nctr - 1
       endif
       if (inotmember >= clump(iclump)%nmember .or. j==k) keep_searching = .false.
       j = j + 1
    enddo
    ekin   = clump(iclump)%kinetic
    epot   = clump(iclump)%potential
    etherm = clump(iclump)%thermal

    !--This is a clump; perform additional bookkeeping
    if (ekin_coef*ekin + epot < 0. .and. ekin > 0.) then  ! ekin criterion is to prevent clumps from being formed from the unmoving initial background
       ! Determine clump id: if sink existed in previous dump, use its imember
       if (isink > 0) then
          if (idclumpsinkold(isink) > 0) then
             imember = idclumpsinkold(isink)
             idclumpsink(isink) = imember
          endif
       endif
       ! Determine clump id: for gas-only clouds, see which id from the previous dump is most used; use if not already in use
       if (imember==0) then
          imax = maxval(iprev,1)
          if (imax > 0) then
             imember = maxloc(iprev,1)
             do jj = 1,iclump-1
                if (clump(jj)%ID==imember) imember = 0  ! to avoid duplicating ID numbers if a cloud essentially split
             enddo
          endif
          if (imember==0) then
             idclumpmax = idclumpmax + 1
             imember    = idclumpmax
          endif
       endif
       do jj = 1,k
          i = icandidate(lst(jj))
          if (idclump(i) == idtmp) idclump(i) = imember
       enddo
       clump(iclump)%ID = imember
       ! calculate clump radius that contains 90% of the mass & update the membership number
       dist2 = 0.
       dxyz  = 0.
       clump(iclump)%dr(:) = 0.
       do j = 1,nctr
          dx = clump(iclump)%r(1)-xyz(1,j)
          dy = clump(iclump)%r(2)-xyz(2,j)
          dz = clump(iclump)%r(3)-xyz(3,j)
#ifdef PERIODIC
          if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
          if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
          if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
          dxyz(1,j) = dx*dx
          dxyz(2,j) = dy*dy
          dxyz(3,j) = dz*dz
          dist2( j) = dxyz(1,j) + dxyz(2,j) + dxyz(3,j)
          do p = 1,3
             clump(iclump)%dr(p) = max(clump(iclump)%dr(p),dxyz(p,j))
          enddo
       enddo
       if (isink > 0) then
          menc = 0.9*(clump(iclump)%mass-xyzmh_ptmass(4,isink))
       else
          menc = 0.9*clump(iclump)%mass
       endif
       ienc  = int(menc/pmassi + 1.0)
       call indexx(nctr,dist2,lst)
       clump(iclump)%size = sqrt(dist2(lst(ienc)))
       do p = 1,3
          clump(iclump)%dr(p) = sqrt(clump(iclump)%dr(p))
       enddo
    endif
    print*, 'eligible particles near candidate particle found, ', 'number of candidates =', k
 else
    COUNTER = COUNTER + 1
    print*, 'no eligible particles near the candidate particle', ii, npart
 endif

 !--Kill the clump if not bound
! do jj = 1, npart
!    print*, lst(1)
!    if (icandidate(jj) < 0) print*, icandidate(jj), jj
!    COUNTER = COUNTER + 1
! end do
 if (imember==0) then
    do jj = 1,k
       i = icandidate(lst(jj))
       if (lst(jj) > SIZE(icandidate)) then
            print*, 'BIGGER', jj, lst(jj), SIZE(icandidate)
            !read(*,*)
       end if
       !if (i < 0) print*, lst(jj)
       if (i > SIZE(idclump)) print*, i, SIZE(idclump), idclump(i)
       if (idclump(i)==idtmp) then
          if (jj < nmin) then
             idclump(i) = -1
          else
             idclump(i) = 0
          endif
       endif
       !COUNTER = COUNTER + 1
    enddo
    call reset_clump(iclump)
    if (isink > 0) then
       clump(iclump)%r(1:3)    = xyzmh_ptmass(1:3,isink)
       clump(iclump)%v(1:3)    = vxyz_ptmass(1:3,isink)
       clump(iclump)%mass      = xyzmh_ptmass(4,isink)
       clump(iclump)%size      = xyzmh_ptmass(ihacc,isink)
       clump(iclump)%pointmass = isink
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
 print*, 'done', imember

end subroutine obtain_membership
!--------------------------------------------------------------------------
!--Zero the clump
subroutine reset_clump(i)
 integer, intent(in) :: i

  clump(i)%r(:)      = 0.0
  clump(i)%v(:)      = 0.0
  clump(i)%dr(:)     = 0.0
  clump(i)%ID        = 0
  clump(i)%mass      = 0.0
  clump(i)%size      = 0.0
  clump(i)%nmember   = 0
  clump(i)%kinetic   = 0.0
  clump(i)%potential = 0.0
  clump(i)%thermal   = 0.0
  clump(i)%magnetic  = 0.0
  clump(i)%pointmass = 0

end subroutine reset_clump
!--------------------------------------------------------------------------
end module
