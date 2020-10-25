
module vars
real,allocatable,save::xy_1_0(:,:),xy_2_0(:,:),xy_1(:,:),xy_2(:,:),xy_0(:,:),xy(:,:),dxy(:,:),dXSHP(:,:)
real,save::L_x_1,R_x_1,D_y_1,U_y_1,L_x_2,R_x_2,D_y_2,U_y_2
real,save::dx_1,dy_1,dx_2,dy_2,dx,dy
integer,save::dim

integer,save::loop_time
real,save::t,dt
real,save::h

real,allocatable,save::rho_0(:),rho(:),rho_temp(:),m(:),tag(:)

real,allocatable,save::rho_1_0(:),rho_2_0(:),rho_1(:),rho_2(:),m_1(:),m_2(:),tag_1(:),tag_2(:)
real::c,s,g,miu,rho_init


real,allocatable,save::buck(:),p_0(:),p(:),e(:),v_0(:,:),v(:,:),v_b(:,:),dv(:,:),a(:,:)!!
real,allocatable,save::v_1_0(:,:),v_2_0(:,:),p_1(:),p_2(:)             
real,allocatable,save::sig(:,:,:),tao(:,:,:),eps_rate(:,:,:),delta(:,:)



integer,allocatable,save::mark(:)
integer,save::num_x_1,num_y_1,num_x_2,num_y_2,num_x,num_y,num_1,num_2,num
real,save::k

end module vars