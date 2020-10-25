program test_sph !a square phase water. movement unit g,cm,s

use imsl
use vars
use subs
use init
implicit none
integer::i,j
real::time_begin,time_end
    call init_output()
	call init_data()
	call cpu_time(time_begin)

do loop_time=1,int(t/dt)
	 print*,"loop_time=",loop_time,"=>",int(t/dt)

   call strain_rate()
 !  call pressure()            !
   call energy()
   call gruneisen()
   call stress()  
   call kinematic()
   call movement()
   call density()             !

end do

	call terminate_output()

	   
end program

