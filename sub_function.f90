module subs
use vars
use imsl
implicit none
contains

	subroutine density()                      !page114页 连续性密度方程
		integer::i,j
		real::d_rho
		real::rv(2)
		real::grad_w(2)
		rho_temp=rho
		do i=1,num
			d_rho=0.0
				do j=1,num
					rv(1)=v(1,i)-v(1,j)
					rv(2)=v(2,i)-v(2,j)
					call gradient_w(i,j,grad_w,h)
					d_rho=d_rho+m(j)*dot_product(rv,grad_w)*dt	 
				end do	
			rho(i)=rho_temp(i)+d_rho
		end do
	end subroutine density

	subroutine strain_rate()
		integer::i,j
		real::grad_w(2)
		real::rv_ji(2)
		    eps_rate=0.0		
			do i=1,num
				do j=1,num
					rv_ji(1)=v(1,j)-v(1,i)
					rv_ji(2)=v(2,j)-v(2,i)
					call gradient_w(i,j,grad_w,h)
					eps_rate(1,1,i)=eps_rate(1,1,i)+m(j)/rho(j)*rv_ji(1)*grad_w(1)+m(j)/rho(j)*rv_ji(1)*grad_w(1)-2/3*m(j)/rho(j)*dot_product(rv_ji,grad_w)
					eps_rate(1,2,i)=eps_rate(1,2,i)+m(j)/rho(j)*rv_ji(2)*grad_w(1)+m(j)/rho(j)*rv_ji(1)*grad_w(2)
					eps_rate(2,1,i)=eps_rate(2,1,i)+m(j)/rho(j)*rv_ji(1)*grad_w(2)+m(j)/rho(j)*rv_ji(2)*grad_w(1)
					eps_rate(2,2,i)=eps_rate(2,2,i)+m(j)/rho(j)*rv_ji(2)*grad_w(2)+m(j)/rho(j)*rv_ji(2)*grad_w(2)-2/3*m(j)/rho(j)*dot_product(rv_ji,grad_w)
				end do
			end do
	end subroutine strain_rate

	subroutine pressure()                                !EOS model
		integer::i
		do i=1,num	
	     !	p(i)=c**2/(rho(i)*7)*((rho(i)/rho_0(i))**7-1)                      !!!!!!
		    p(i)=c**2*rho(i)
        end do
	end subroutine pressure


	subroutine stress()
		integer::i,j
			do i=1,num
			tao(:,:,i)=miu*eps_rate(:,:,i)
			sig(:,:,i)=-p(i)*delta(:,:)+tao(:,:,i)
			end do
	end subroutine stress

    subroutine kinematic()
        integer::i,j
		real::grad_w(2),w
		real::q
		dv=0.0
		q=0.0
			do i=1,num
				do j=1,num
						call gradient_w(i,j,grad_w,h)
    			!		call viscosity(i,j,h,c,rho(i),rho(j),q)
					
							dv(1,i)=dv(1,i)+m(j)*((sig(1,1,i)/(rho(i)**2)+sig(1,1,j)/(rho(j)**2)+q)*grad_w(1))
							dv(2,i)=dv(2,i)+m(j)*((sig(2,2,i)/(rho(i)**2)+sig(2,2,j)/(rho(j)**2)+q)*grad_w(2))
					!       dv(1,i)=dv(1,i)+m(j)*((sig(1,1,i)/(rho(i)**2)+sig(1,1,j)/(rho(j)**2)+q)*grad_w(1)+(sig(1,2,i)/(rho(i)**2)+sig(1,2,j)/(rho(j)**2))*grad_w(2))
					!		dv(2,i)=dv(2,i)+m(j)*((sig(2,1,i)/(rho(i)**2)+sig(2,1,j)/(rho(j)**2))*grad_w(1)+(sig(2,2,i)/(rho(i)**2)+sig(2,2,j)/(rho(j)**2)+q)*grad_w(2))					
					
				end do			
			end do    
					dv=dv*dt
    end subroutine kinematic

	
	subroutine viscosity(i,j,h_ij,c_ij,rho_i,rho_j,q)
		integer::i,j
		real::alpha,beta,phi,theta
		real::h_i,h_j,h_ij					   !h_j=h_j=h_ij=constant
		real::c_i,c_j,c_ij                     !c_j=c_j=c_ij=constant
		real::rho_i,rho_j,rho_ij
		real::rv(2)
		real::vd(2)
		real::q

			rv(1)=v(1,i)-v(1,j)
			rv(2)=v(2,i)-v(2,j)
			vd(1)=xy(1,i)-xy(1,j)
			vd(2)=xy(2,i)-xy(2,j)	
					
		    alpha=0.000001
			beta=0.000001
			phi=0.1*h_ij
			rho_ij=(rho_i+rho_j)/2	         

			theta=(h_ij*dot_product(rv,vd))/(dot_product(vd,vd)**2+phi**2)		

		if(dot_product(rv,vd)<0)then

				q=(-alpha*c_ij*theta+beta*theta**2)/rho_ij
	!			q=0.0
			else 
				q=0.0
			end if
		return
	end subroutine viscosity

   subroutine energy()
	 integer::i,j
	 real::grad_w(2)
     real::de
	 real::rv(2)
	 de=0.0
	 	do i=1,num
	   	   do j=1,num
				rv(1)=v(1,i)-v(1,j)
				rv(2)=v(2,i)-v(2,j)
				call gradient_w(i,j,grad_w,h)
				de=de+0.5*m(j)*(((sig(1,1,i)/(rho(i)**2)+sig(1,1,j)/(rho(j)**2))*grad_w(1)+(sig(1,2,i)/(rho(i)**2)+sig(1,2,j)/(rho(j)**2))*grad_w(2))*(-rv(1))+((sig(2,1,i)/(rho(i)**2)+sig(2,1,j)/(rho(j)**2))*grad_w(1)+(sig(2,2,i)/(rho(i)**2)+sig(2,2,j)/(rho(j)**2))*grad_w(2))*(-rv(2)))
				e(i)=e(i)+de*dt			
			end do
        end do    


   end subroutine energy

   subroutine gruneisen()
		integer::i
		real::ph,eh
		real::v_0,v
		do i=1,num	
			v_0=rho_0(i)
			v=rho(i)
			ph=c**2*(v_0-v)/(v_0-s*(v_0-v))**2
			eh=0.5*ph*(v_0-v)
			p(i)=ph+g*rho(i)*(e(i)-eh)
			if(i==1)then
				 print*,ph
			end if
        end do
   end subroutine gruneisen


	subroutine movement()
	  integer::i,j
	  real::w
	  real::rv(2)
		do i=1,num
			do j=1,2
				call XSHP(i,j,w,h)					                                !movement part                                     !
				v_b(j,i)=v(j,i) 
				v(j,i)=v_b(j,i)+dv(j,i)                                          !
				dxy(j,i)=0.5*(v(j,i)+v_b(j,i))*dt                                         !
				xy(j,i)=xy(j,i)+dxy(j,i)
			end do
		end do   
	end subroutine movement


	function w(i,j,h)
		integer::i,j
		real::h	          !kernal distance
       
		real::pi	
		real::ad	
		real::w

		real::vd(2)       !vector  i->j
		real::s           !absolute value of vd_ij=distance of these pair particle
		real::R           !relative distance =s/h
		
		pi=3.1415         !constant of pi
 
		vd(1)=xy(1,i)-xy(1,j)
		vd(2)=xy(2,i)-xy(2,j)
		
		s=sqrt(dot_product(vd,vd))
		R=s/h
		
		ad=5/(pi*h**2)
		if(R>1) ad=0.0
		w=ad*(1+3*R)*(1-R)**3
		return
	end function w

	subroutine gradient_w_1(i,j,grad_w,h)
	
		integer::i,j
		real::h	          !kernal distance
		real::grad_w(2)   !grad vector of w

		real::pi	
		real::ad	
		real::w

		real::vd(2)       !vector  i->j
		real::s           !absolute value of vd_ij=distance of these pair particle
		real::R           !relative distance =s/h
		real::grad        !grad of kernal function w
		
		pi=3.1415         !constant of pi

		vd(1)=xy(1,i)-xy(1,j)
		vd(2)=xy(2,i)-xy(2,j)
		grad_w=0.0
		s=sqrt(dot_product(vd,vd))
	
		R=s/h

		ad=5/(pi*h**2)
		if(R>1) ad=0.0

		w=ad*(1+3*R)*(1-R)**3
		grad=ad*(3*(1-R)**3-3*(1-R)**2*(1+3*R))*(1/h)*(1/s)

		grad_w=grad*vd
		if(i==j) grad_w=0.0	
		return
	
	end subroutine gradient_w_1

	subroutine gradient_w(i,j,grad_w,h)
	
		integer::i,j
		real::h	          !kernal distance
		real::grad_w(2)   !grad vector of w

		real::pi	
		real::ad	
		real::w

		real::vd(2)       !vector  i->j
		real::s           !absolute value of vd_ij=distance of these pair particle
		real::R           !relative distance =s/h
		real::grad        !grad of kernal function w
		
		pi=3.1415         !constant of pi

		vd(1)=xy(1,i)-xy(1,j)
		vd(2)=xy(2,i)-xy(2,j)
		grad_w=0.0
		s=sqrt(dot_product(vd,vd))
	
		R=s/h

		ad=5/(pi*h**2)
		if(R>1) ad=0.0

		w=ad*(1+3*R)*(1-R)**3
		grad=ad*(3*(1-R)**3-3*(1-R)**2*(1+3*R))*(1/h)*(1/s)

		grad_w=grad*vd
		if(i==j) grad_w=0.0	
		return
	
	end subroutine gradient_w



    subroutine XSHP(i,j,w,h)
		integer::i,j
		real::h	          !kernal distance

		real::pi	
		real::ad	
		real::w

		real::vd(2)       !vector  i->j
		real::s           !absolute value of vd_ij=distance of these pair particle
		real::R           !relative distance =s/h
		real::grad        !grad of kernal function w
		
		pi=3.1415         !constant of pi

		vd(1)=xy(1,i)-xy(1,j)
		vd(2)=xy(2,i)-xy(2,j)
		s=sqrt(dot_product(vd,vd))
	
		R=s/h

		ad=5/(pi*h**2)
		if(R>1) ad=0.0
		w=ad*(1+3*R)*(1-R)**3

		return

	end subroutine

subroutine init_output()
integer::i,j
	open(100,file="ans.txt")
	do i=1,num
		write(100,*)xy(1,i),xy(2,i),v(1,i),p(i),tag(i)
	end do
end subroutine

subroutine terminate_output()
	integer::i,j
	open(101,file="ans_1.txt")
	do i=1,num
		write(101,*)xy(1,i),xy(2,i),rho(i),p(i),tag(i)
	end do
end subroutine terminate_output



end module subs