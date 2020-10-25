module init
use vars
use imsl
use subs
implicit none
contains

subroutine init_data()
	call allocate_number()
	call position_init()                             !用ansys读取个数据 两个半球mesh后节点,目前先用两个平面靶看看
    call velocity_init()
	call density_init()
	call pressure_e_init()
	call merge_init()
	call calculate_init()
end subroutine init_data

subroutine  allocate_number()
	
	dim=2
	num_x_1=20
	num_y_1=20
	num_x_2=20
	num_y_2=20

	num_x=num_x_1+num_x_2
	num_y=num_y_1
	
	num_1=num_x_1*num_y_1
	num_2=num_x_2*num_y_2
	num=num_x*num_y

	allocate(xy_1_0(dim,num_1))
	allocate(xy_1(dim,num_1))
	allocate(rho_1_0(num_1))
	allocate(rho_1(num_1))
	allocate(m_1(num_1))
	allocate(tag_1(num_1))
	allocate(v_1_0(dim,num_1))
	allocate(p_1(num_1))

	allocate(xy_2_0(dim,num_2))
	allocate(xy_2(dim,num_2))
	allocate(rho_2_0(num_2))
	allocate(rho_2(num_2))
	allocate(m_2(num_2))
	allocate(tag_2(num_2))
	allocate(v_2_0(dim,num_2))
	allocate(p_2(num_2))

	allocate(xy_0(2,num))
	allocate(xy(2,num))
	allocate(dxy(dim,num))
	allocate(v_0(dim,num))
	allocate(v(dim,num))
	allocate(v_b(dim,num))
	allocate(dv(dim,num))
	allocate(a(dim,num))
	allocate(dXSHP(dim,num))

	allocate(rho_0(num))
	allocate(rho(num))
	allocate(rho_temp(num))
	allocate(m(num))
	allocate(buck(num))
	allocate(p_0(num))
	allocate(p(num))
	allocate(e(num))



	allocate(tag(num))
	allocate(sig(dim,dim,num))
	allocate(eps_rate(dim,dim,num))	
	allocate(tao(dim,dim,num))
	allocate(delta(dim,dim))	
	allocate(mark(num))


end subroutine allocate_number


    
subroutine position_init()                   !两块矩形靶板，尺寸为1*1，中心位置为(0.5,0.5) and (3.5,0.5)
	integer::i,j

	L_x_1=0.0
	R_x_1=1.0
	D_y_1=0.0
	U_y_1=1.0

	L_x_2=1.1
	R_x_2=2.1
	D_y_2=0.0
	U_y_2=1.0



dx_1=abs(R_x_1-L_x_1)/num_x_1
dy_1=abs(U_y_1-D_y_1)/num_y_1

dx_2=abs(R_x_2-L_x_2)/num_x_2
dy_2=abs(U_y_2-D_y_2)/num_y_2


do i=1,num_x_1
	do j=1,num_y_1
		xy_1_0(1,(i-1)*num_y_1+j)=mod(i,num_x_1)*dx_1+L_x_1
		xy_1_0(2,(i-1)*num_y_1+j)=mod(j,num_y_1)*dy_1+D_y_1
	end do
end do

do i=1,num_x_2
	do j=1,num_y_2
		xy_2_0(1,(i-1)*num_y_2+j)=mod(i,num_x_2)*dx_2+L_x_2
		xy_2_0(2,(i-1)*num_y_2+j)=mod(j,num_y_2)*dy_2+D_y_2
	end do
end do

xy_1=xy_1_0
xy_2=xy_2_0

end subroutine position_init

subroutine velocity_init()
	integer::i,j
	do i=1,num_x_1
		do j=1,num_y_1
			v_1_0(1,(i-1)*num_x_1+j)=0.05         !100m/s
			v_1_0(2,(i-1)*num_x_1+j)=0.0 
			tag_1((i-1)*num_x_1+j)=1.0
		end do
	end do

	do i=1,num_x_2
		do j=1,num_y_2
			v_2_0(1,(i-1)*num_x_2+j)=0.0          !100m/s
			v_2_0(2,(i-1)*num_x_2+j)=0.00 
			tag_2((i-1)*num_x_2+j)=2.0
		end do
	end do
end subroutine velocity_init

subroutine density_init()
	rho_init=18.087
	rho_1_0=rho_init
	rho_2_0=rho_init
	rho_1=rho_1_0
	rho_2=rho_2_0
	m_1=dx_1*dy_1*rho_init
	m_2=dx_2*dy_2*rho_init
end subroutine density_init

subroutine pressure_e_init()
!c=0.251
c=0.0251
s=0.151
!s=1.51
g=2.03
miu=0.0
e=0.0
p_1=0.0
p_2=0.0

end subroutine pressure_e_init

subroutine merge_init()
	integer::i,j
	
	dx=dx_1
	dy=dy_1


	do i=1,num_x_1*num_y_1
		xy_0(:,i)=xy_1_0(:,i)
		v_0(:,i)=v_1_0(:,i)
		m(i)=m_1(i)
		rho_0(i)=rho_1_0(i)
		p_0(i)=p_1(i)
		tag(i)=tag_1(i)
	end do
	
	do i=num_x_1*num_y_1+1,num_x*num_y
		xy_0(:,i)=xy_2_0(:,i-num_x_1*num_y_1)
		v_0(:,i)=v_2_0(:,i-num_x_1*num_y_1)
		m(i)=m_2(i-num_x_1*num_y_1)
		rho_0(i)=rho_2_0(i-num_x_1*num_y_1)
		p_0(i)=p_2(i-num_x_1*num_y_1)
		tag(i)=tag_2(i-num_x_1*num_y_1)
	end do
	xy=xy_0
	v=v_0
	rho=rho_0
	p=p_0
end subroutine merge_init

subroutine calculate_init()
	h=1.5*dx                      ! smooth radius
	dt=0.02
	t=3.0
	num=num_x*num_y
	delta=0.0
	delta(1,1)=1.0
	delta(2,2)=1.0
end subroutine calculate_init


end module init