module analysis
!
! Analysis routine which compares the positions of particles between two dumps
!
 use dim, only: maxp_hard
 implicit none
 character(len=20), parameter, public :: analysistype = 'positions'
 character(len=50) :: dump1, dump2
 public :: do_analysis
 logical :: first_call = .true.
 real :: x1(maxp_hard), x2(maxp_hard), v1(maxp_hard), v2(maxp_hard), t1, t2 ! Have to define them up here so they are shared between dumps
 

contains

subroutine do_analysis(dumpfile, num, xyzh, vxyzu, particlemass, npart, time, iunit) ! These are required parameters when making do_analysis from the overarching phantomanalysis.f90 script
 use units, only: udist, utime
 character(len=*), intent(in) :: dumpfile
 integer, intent(in) :: npart, num, iunit
 real, intent(in) :: xyzh(:,:), vxyzu(:,:)
 real, intent(in) :: particlemass, time

 integer :: i, j, k
 real :: xi, vxi
 real :: delx(npart), delt, v_unit
 character(len=200) :: fileout

 v_unit = udist/utime
 print*, num

 if(first_call) then
    t1 = time
    dump1 = dumpfile
 else
    t2 = time
    dump2 = dumpfile
 endif

 do i = 1, npart
    xi = xyzh(1, i)
    vxi = vxyzu(1, i)
    if (first_call) then
        x1(i) = xi
        v1(i) = vxi
    else
        x2(i) = xi
        v2(i) = vxi
    endif
 enddo

 if (.not. first_call) then
    delt = abs(t2 - t1)
    print*, delt
    do i = 1, npart
        delx(i) = abs(x1(i) - x2(i))
        !print*, delx(i)
    enddo
 endif

 x2 = x2*udist
 delx = delx*udist
 delt = delt*utime

 if (.not. first_call) then
    write(fileout,'(5a)') 'dump_pos_comp_',dump1(index(dump1,'0'):index(dump1,'0')+4), &
    'vs',dump2(index(dump2,'0'):index(dump2,'0')+4),'.dat' ! xa tells you how many characters to expect
    fileout=trim(fileout)
    open(iunit,file=fileout)
    do i = 1, npart
        write(iunit,'(3(1pe18.10,1x))') x2(i), (delx(i)/delt), &
        abs(((v2(i) + v1(i))/2))*(udist/utime)
    enddo
 endif
 close(iunit)

 first_call = .false.

 end subroutine do_analysis

end module analysis

