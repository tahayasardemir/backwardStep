	Program FlowOverBackwardFacingStep
c..Taha Yaşar Demir / 1881978
c..CE-580 / Homework #9
	parameter(mx=1001)
	common/init/  rL,H,SL,SH,Um,vis,rho,gamma,dt
	common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >		  NS,MS
	common/flow/  u(mx,mx),v(mx,mx),p(mx,mx),S(mx,mx)
	common/xcoef/ ue(mx,mx),uw(mx,mx),un(mx,mx),us(mx,mx),conu(mx,mx),
     >              difu(mx,mx),F(mx,mx),gx(mx,mx),qxe(mx,mx),qxw(mx,mx)
     >              ,qxn(mx,mx),qxs(mx,mx)
	common/ycoef/ ve(mx,mx),vw(mx,mx),vn(mx,mx),vs(mx,mx),conv(mx,mx),
     >		  difv(mx,mx),G(mx,mx),gy(mx,mx),qye(mx,mx),qyw(mx,mx)
     >              ,qyn(mx,mx),qys(mx,mx)
      common/bndry/ Un0(mx),Vn0(mx),Unp1(mx),Vnp1(mx),Pn0(mx),Pnp1(mx),
     >              Um0(mx),Vm0(mx),Ump1(mx),Vmp1(mx),Pm0(mx),Pmp1(mx),
     >              Fn0(mx),Gm0(mx)
      common/Err/   E,res,tolerance

	open(11,file='mesh.tec',form="formatted")
	open(13,file="error.dat")
	open(14,file="Residual.dat")
	E = 1.
	tolerance = 1e-8 
	call Initialize
	l = 1
c	do l=1,10000
	do while(E.gt.tolerance.and.l.lt.60000)
		call Boundary
		call Coefficients
		call Poisson
		call UpdateVelocities
		call Error(l)
		l = l+1
		if(mod(l,1000).eq.0) print*, l,E,res/(N*M)
	enddo
	print*, l
	call Output
	close(11)
	close(13)
	close(14)
	stop
	end
c-----------------------------------------------------------------------
	subroutine Initialize
	parameter(mx=1001)
	common/init/  rL,H,SL,SH,Um,vis,rho,gamma,dt
	common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >		  NS,MS

	rL  = 2.  !m /total length
	H   = 0.1 !m /total height 
	SL  = 0.4 !m /step length
	SH  = 0.05!m /step height
	Um  = 2.  !m/s 0.2--1.0--2.0
	vis = 1e-4!m^2/s 
	rho = 1e3 !km/m^3

	gamma = 1.8 ! Over-Relaxation Parameter

	call GenerateGrid

	dt = 5.e-6

	call InitialCondition

	return
	end
c-----------------------------------------------------------------------
	subroutine GenerateGrid
	parameter(mx=1001)
	common/init/  rL,H,SL,SH,Um,vis,rho,gamma,dt
	common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >		  NS,MS

	N  = 300 ! Cell count in x-direction
	M  = 60  ! Cell count in y-direction
	dx = 2.0/N 
	dy = 0.1/M
	NS = ceiling(0.4/dx)  +1 
	MS = ceiling(0.05/dy) +1

	do j = 1,M+1
		x(1,j) = 0.
		do i=1,N
			x(i+1,j) = x(i,j) + dx
		enddo
		x(N+1,j) = rL
	enddo

	do i = 1,N+1
		y(i,1) = 0
		do j=1,M
			y(i,j+1) = y(i,j) + dy
		enddo
		y(i,M+1) = H
	enddo

	do j = 1,M
		do i = 1,N
			xc(i,j) = (x(i+1,j)+x(i,j))/2
			yc(i,j) = (y(i,j+1)+y(i,j))/2
		enddo
	enddo

	a = (2/(dx**2)) + (2/(dy**2))


	return
	end

c-----------------------------------------------------------------------
	subroutine InitialCondition
	parameter(mx=1001)
	common/init/  rL,H,SL,SH,Um,vis,rho,gamma,dt
	common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >		  NS,MS
	common/bndry/ Un0(mx),Vn0(mx),Unp1(mx),Vnp1(mx),Pn0(mx),Pnp1(mx),
     >              Um0(mx),Vm0(mx),Ump1(mx),Vmp1(mx),Pm0(mx),Pmp1(mx),
     >              Fn0(mx),Gm0(mx)
	common/flow/  u(mx,mx),v(mx,mx),p(mx,mx),S(mx,mx)
	real uniflow(mx)

	do i=1,N+2
		do j=1,M+2
			u(i,j) = 0.
			v(i,j) = 0.
			p(i,j) = 0.
		enddo
	enddo

	do i=1,N+1
	if(i.lt.NS) then
	k=1
	do j = MS+1,M+1
      	uniflow(k) = (-((0.1*Um/(0.025**2))*yc(1,k)**2)+
     +	 	       ((0.2*Um/0.025)*yc(1,k)))/0.1
      	u(i,j) = uniflow(k) ! wont updated in future iterations
      	v(1,j) = 0. ! wont updated in future iterations
      	k = k+1
      enddo
      else
      k=1
	do j = 1,M+1
      	uniflow(k) = (-((0.05*Um/(0.05**2))*yc(1,k)**2)+
     +                      ((0.1*Um/0.05)*yc(1,k)))/0.1
      	u(i,j) = uniflow(k) ! wont updated in future iterations
      	v(1,j) = 0. ! wont updated in future iterations
      	k = k+1
      enddo
      endif
      enddo

      write(11,*) ' variables="x","y","u","v" '
	write(11,*) ' zone i=',N+1, 'j=',M+1
	do j = 1,M+1
	do i = 1,N+1
		write(11,*) x(i,j),y(i,j),u(i,j),v(i,j)
	enddo
	enddo

	return
	end
c-----------------------------------------------------------------------
	subroutine Boundary
	parameter(mx=1001)
	common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >		  NS,MS
      common/flow/  u(mx,mx),v(mx,mx),p(mx,mx),S(mx,mx)
	common/xcoef/ ue(mx,mx),uw(mx,mx),un(mx,mx),us(mx,mx),conu(mx,mx),
     >              difu(mx,mx),F(mx,mx),gx(mx,mx),qxe(mx,mx),qxw(mx,mx)
     >              ,qxn(mx,mx),qxs(mx,mx)
	common/ycoef/ ve(mx,mx),vw(mx,mx),vn(mx,mx),vs(mx,mx),conv(mx,mx),
     >		  difv(mx,mx),G(mx,mx),gy(mx,mx),qye(mx,mx),qyw(mx,mx)
     >              ,qyn(mx,mx),qys(mx,mx)
	common/bndry/ Un0(mx),Vn0(mx),Unp1(mx),Vnp1(mx),Pn0(mx),Pnp1(mx),
     >              Um0(mx),Vm0(mx),Ump1(mx),Vmp1(mx),Pm0(mx),Pmp1(mx),
     >              Fn0(mx),Gm0(mx)

      ! Inflow Boundary Conditions
      do j = 1,M+2
      	p(1,j) = p(2,j)
      	F(1,j) = u(1,j)
      enddo
      ! Outflow Boundary Conditions
      do j = 1,M+1
      	v(N+2,j) = 0.
c      	Unp1(j) = u(N,j)
      	u(N+1,j) = u(N,j)
      	F(N+1,j) = u(N+1,j)
      	p(N+2,j) = p(N+1,j)
      enddo
            u(N+1,M+2) = u(N,M+2)
      	F(N+1,M+2) = u(N+1,M+2)
      	p(N+2,M+2) = p(N+1,M+2)
      !! Wall Boundaries !!
      ! Horizaontal Part of the Step
      do i=1,NS-1
      	v(i,MS)   = 0. ! on the wall
      	u(i,MS)   =-u(i,MS+1) ! ghost cell
      	G(i,MS)   = v(i,MS)
      	p(i,MS)   = p(i,MS+1)
      enddo
      	v(NS,MS)   = 0.
      	G(NS,MS)   = v(NS,MS)
      	p(NS,MS)   = p(NS,MS+1) !!!!!!!!!!!!!!!!!!!!!!

      ! Vertical Part of the Step
      do j=1,MS-1
      	v(NS,j)   =-v(NS+1,j) ! ghost cell
      	u(NS,j)   = 0.  ! on the wall
      	F(NS,j)   = u(NS,j)
      	p(NS,j)   = p(NS+1,j)
      enddo
      	p(NS,MS)   = p(NS+1,MS) !!!!!!!!!!!!!!!!!!!!!
      ! Lower Wall
      do i=NS,N+2 ! maybe add v(i,1) = 0  and G(i,1) too
      	v(i,1) = 0.
      	u(i,1) =-u(i,2)
      	G(i,1) = v(i,1)
      	p(i,1) = p(i,2)
      enddo
      ! Upper Wall
      do i=1,N+2
      	v(i,M+1) = 0.
      	G(i,M+1) = v(i,M+1)
      	p(i,M+2) = p(i,M+1)
      enddo
      do i=1,N+1
      	u(i,M+2) =-u(i,M+1)
      enddo

     	return
     	end

c-----------------------------------------------------------------------
	subroutine Coefficients
	parameter(mx=1001)
	common/init/  rL,H,SL,SH,Um,vis,rho,gamma,dt
	common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >		  NS,MS
	common/flow/  u(mx,mx),v(mx,mx),p(mx,mx),S(mx,mx)
	common/xcoef/ ue(mx,mx),uw(mx,mx),un(mx,mx),us(mx,mx),conu(mx,mx),
     >              difu(mx,mx),F(mx,mx),gx(mx,mx),qxe(mx,mx),qxw(mx,mx)
     >              ,qxn(mx,mx),qxs(mx,mx)
	common/ycoef/ ve(mx,mx),vw(mx,mx),vn(mx,mx),vs(mx,mx),conv(mx,mx),
     >		  difv(mx,mx),G(mx,mx),gy(mx,mx),qye(mx,mx),qyw(mx,mx)
     >              ,qyn(mx,mx),qys(mx,mx)
      common/bndry/ Un0(mx),Vn0(mx),Unp1(mx),Vnp1(mx),Pn0(mx),Pnp1(mx),
     >              Um0(mx),Vm0(mx),Ump1(mx),Vmp1(mx),Pm0(mx),Pmp1(mx),
     >              Fn0(mx),Gm0(mx)


     	do i=1,N+2
     		do j=1,M+2
     			gx(i,j) = 0.9
     			gy(i,j) = 0.9
     		enddo
     	enddo	
! u-coefficients
      DO i =2,N
      IF(i.le.NS) THEN ! Before Step
      	do j=MS+1,M+1
     		qxe(i,j) = 2*(u(i,j)+u(i+1,j))*dy
     		qxw(i,j) = 2*(u(i-1,j)+u(i,j))*dy
     		qxn(i,j) = 2*(v(i,j)+v(i+1,j))*dx
     		qxs(i,j) = 2*(v(i,j-1)+v(i+1,j-1))*dx
     		ue(i,j)  = (u(i,j)+u(i+1,j)+gx(i,j)*sign(1.,qxe(i,j))
     + 		    *(u(i,j)-u(i+1,j)))/2.
		un(i,j)  = (u(i,j)+u(i,j+1)+gx(i,j)*sign(1.,qxn(i,j))
     + 		    *(u(i,j)-u(i,j+1)))/2.
      	uw(i,j)  = (u(i-1,j)+u(i,j)+gx(i,j)*sign(1.,qxw(i,j))
     +                *(u(i-1,j)-u(i,j)))/2.
      	us(i,j)  = (u(i,j-1)+u(i,j)+gx(i,j)*sign(1.,qxs(i,j))
     +                *(u(i,j-1)-u(i,j)))/2.
		difu(i,j)= vis*( ( ((u(i+1,j)-u(i,j))/dx)
     +		           - ((u(i,j)-u(i-1,j))/dx) )*dy
     +		          +( ((u(i,j+1)-u(i,j))/dy)
     +		           - ((u(i,j)-u(i,j-1))/dy) )*dx )
     		conu(i,j)= ue(i,j)*qxe(i,j)-uw(i,j)*qxw(i,j)+
     +		     un(i,j)*qxn(i,j)-us(i,j)*qxs(i,j)
		F(i,j)   = u(i,j) + (dt/(dx*dy))*(difu(i,j)-conu(i,j))
      	enddo
      ELSE ! After Step
      	do j=2,M+1
     		qxe(i,j) = 2*(u(i,j)+u(i+1,j))*dy
     		qxw(i,j) = 2*(u(i-1,j)+u(i,j))*dy
     		qxn(i,j) = 2*(v(i,j)+v(i+1,j))*dx
     		qxs(i,j) = 2*(v(i,j-1)+v(i+1,j-1))*dx
     		ue(i,j)  = (u(i,j)+u(i+1,j)+gx(i,j)*sign(1.,qxe(i,j))
     +		    *(u(i,j)-u(i+1,j)))/2.
		un(i,j)  = (u(i,j)+u(i,j+1)+gx(i,j)*sign(1.,qxn(i,j))
     + 		    *(u(i,j)-u(i,j+1)))/2.
     		uw(i,j)  = (u(i-1,j)+u(i,j)+gx(i,j)*sign(1.,qxw(i,j))
     +                *(u(i-1,j)-u(i,j)))/2.
     		us(i,j)  = (u(i,j-1)+u(i,j)+gx(i,j)*sign(1.,qxs(i,j))
     +                *(u(i,j-1)-u(i,j)))/2.
     		difu(i,j)= vis*( ( ((u(i+1,j)-u(i,j))/dx)
     +		           - ((u(i,j)-u(i-1,j))/dx) )*dy
     +		          +( ((u(i,j+1)-u(i,j))/dy)
     +		           - ((u(i,j)-u(i,j-1))/dy) )*dx )
     		conu(i,j)= ue(i,j)*qxe(i,j)-uw(i,j)*qxw(i,j)+
     +		     un(i,j)*qxn(i,j)-us(i,j)*qxs(i,j)
     		F(i,j)   = u(i,j) + (dt/(dx*dy))*(difu(i,j)-conu(i,j))
      	enddo
      ENDIF
      ENDDO
! v-coefficients
	DO i=2,N+1
	IF (i.le.NS) THEN
		do j=MS+1,M ! Before Step
     		qye(i,j) = 2*(u(i,j)+u(i,j+1))*dy
     		qyw(i,j) = 2*(u(i-1,j)+u(i-1,j+1))*dy
     		qyn(i,j) = 2*(v(i,j)+v(i,j+1))*dx
     		qys(i,j) = 2*(v(i,j)+v(i,j-1))*dx
     		ve(i,j)  = (v(i,j)+v(i+1,j)+gy(i,j)*sign(1.,qye(i,j))
     +		    *(v(i,j)-v(i+1,j)))/2.
     		vw(i,j)  = (v(i-1,j)+v(i,j)+gy(i,j)*sign(1.,qyw(i,j))
     +  		    *(v(i-1,j)-v(i,j)))/2.
     		vn(i,j)  = (v(i,j)+v(i,j+1)+gy(i,j)*sign(1.,qyn(i,j))
     +		    *(v(i,j)-v(i,j+1)))/2.
     		vs(i,j)  = (v(i,j-1)+v(i,j)+gy(i,j)*sign(1.,qys(i,j))
     +		    *(v(i,j-1)-v(i,j)))/2.
     		difv(i,j)= vis*( ( ((v(i+1,j)-v(i,j))/dx)
     +		           - ((v(i,j)-v(i-1,j))/dx) )*dy
     +		          +( ((v(i,j+1)-v(i,j))/dy)
     +		           - ((v(i,j)-v(i,j-1))/dy) )*dx )
     		conv(i,j)= ve(i,j)*qye(i,j)-vw(i,j)*qyw(i,j)+
     +		     vn(i,j)*qyn(i,j)-vs(i,j)*qys(i,j)
     		G(i,j)   = v(i,j) + (dt/(dx*dy))*(difv(i,j)-conv(i,j))
     		enddo
     	ELSE
     		do j=2,M
     		qye(i,j) = 2*(u(i,j)+u(i,j+1))*dy
     		qyw(i,j) = 2*(u(i-1,j)+u(i-1,j+1))*dy
     		qyn(i,j) = 2*(v(i,j)+v(i,j+1))*dx
     		qys(i,j) = 2*(v(i,j)+v(i,j-1))*dx
     		ve(i,j)  = (v(i,j)+v(i+1,j)+gy(i,j)*sign(1.,qye(i,j))
     +		    *(v(i,j)-v(i+1,j)))/2.
     		vw(i,j)  = (v(i-1,j)+v(i,j)+gy(i,j)*sign(1.,qyw(i,j))
     +  		    *(v(i-1,j)-v(i,j)))/2.
     		vn(i,j)  = (v(i,j)+v(i,j+1)+gy(i,j)*sign(1.,qyn(i,j))
     +		    *(v(i,j)-v(i,j+1)))/2.
     		vs(i,j)  = (v(i,j-1)+v(i,j)+gy(i,j)*sign(1.,qys(i,j))
     +		    *(v(i,j-1)-v(i,j)))/2.
     		difv(i,j)= vis*( ( ((v(i+1,j)-v(i,j))/dx)
     +		           - ((v(i,j)-v(i-1,j))/dx) )*dy
     +		          +( ((v(i,j+1)-v(i,j))/dy)
     +		           - ((v(i,j)-v(i,j-1))/dy) )*dx )
     		conv(i,j)= ve(i,j)*qye(i,j)-vw(i,j)*qyw(i,j)+
     +	           vn(i,j)*qyn(i,j)-vs(i,j)*qys(i,j)
     		G(i,j)   = v(i,j) + (dt/(dx*dy))*(difv(i,j)-conv(i,j))     
     		enddo
     	ENDIF
     	ENDDO

	return
	end

c-----------------------------------------------------------------------
	subroutine Poisson
	parameter(mx=1001)
	common/init/  rL,H,SL,SH,Um,vis,rho,gamma,dt
	common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >		  NS,MS
	common/flow/  u(mx,mx),v(mx,mx),p(mx,mx),S(mx,mx)
	common/xcoef/ ue(mx,mx),uw(mx,mx),un(mx,mx),us(mx,mx),conu(mx,mx),
     >              difu(mx,mx),F(mx,mx),gx(mx,mx),qxe(mx,mx),qxw(mx,mx)
     >              ,qxn(mx,mx),qxs(mx,mx)
	common/ycoef/ ve(mx,mx),vw(mx,mx),vn(mx,mx),vs(mx,mx),conv(mx,mx),
     >		  difv(mx,mx),G(mx,mx),gy(mx,mx),qye(mx,mx),qyw(mx,mx)
     >              ,qyn(mx,mx),qys(mx,mx)
      common/bndry/ Un0(mx),Vn0(mx),Unp1(mx),Vnp1(mx),Pn0(mx),Pnp1(mx),
     >              Um0(mx),Vm0(mx),Ump1(mx),Vmp1(mx),Pm0(mx),Pmp1(mx),
     >              Fn0(mx),Gm0(mx)
      common/Err/   E,res,tolerance

      ! Source Term Calculation
      DO i=2,N+1
      IF (i.le.NS) THEN ! Before Step
      	do j=MS,M+1
      	S(i,j) = rho*((F(i,j)-F(i-1,j))/dx+(G(i,j)-G(i,j-1))/dy)/dt
       	enddo
      ELSE ! After Step
      	do j=2,M+1
      	S(i,j) = rho*((F(i,j)-F(i-1,j))/dx+(G(i,j)-G(i,j-1))/dy)/dt
      	enddo
      ENDIF
      ENDDO
      res = 0.
      ! Pressure Update 
      do k=1,30
     		
c     	res = 0. ! for determination of gamma	

      DO i=2,N+1
      IF (i.le.NS) THEN
      	do j=MS+1,M+1
            p(i,j) = p(i,j)+ gamma*(((((p(i+1,j)+p(i-1,j))/(dx**2) + 
     + 	        (p(i,j+1)+p(i,j-1))/(dy**2))-S(i,j)))/a - p(i,j))
            res = res + abs((((((p(i+1,j)+p(i-1,j))/(dx**2) + 
     + 	        (p(i,j+1)+p(i,j-1))/(dy**2))-S(i,j)))/a - p(i,j)))
      	enddo
      ELSE
      	do j=2,M+1
            p(i,j) = p(i,j)+ gamma*(((((p(i+1,j)+p(i-1,j))/dx**2 + 
     + 	        (p(i,j+1)+p(i,j-1))/dy**2)-S(i,j)))/a - p(i,j))
            res = res + abs((((((p(i+1,j)+p(i-1,j))/(dx**2) + 
     + 	        (p(i,j+1)+p(i,j-1))/(dy**2))-S(i,j)))/a - p(i,j)))
      	enddo
      ENDIF
      ENDDO
c      if(E.eq.1.) write(14,*) k,res/(N*M) ! for determination of gamma
      enddo


      return
      end


c-----------------------------------------------------------------------
	subroutine UpdateVelocities
	parameter(mx=1001)
	common/init/  rL,H,SL,SH,Um,vis,rho,gamma,dt
	common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >		  NS,MS
	common/flow/  u(mx,mx),v(mx,mx),p(mx,mx),S(mx,mx)
	common/xcoef/ ue(mx,mx),uw(mx,mx),un(mx,mx),us(mx,mx),conu(mx,mx),
     >              difu(mx,mx),F(mx,mx),gx(mx,mx),qxe(mx,mx),qxw(mx,mx)
     >              ,qxn(mx,mx),qxs(mx,mx)
	common/ycoef/ ve(mx,mx),vw(mx,mx),vn(mx,mx),vs(mx,mx),conv(mx,mx),
     >		  difv(mx,mx),G(mx,mx),gy(mx,mx),qye(mx,mx),qyw(mx,mx)
     >              ,qyn(mx,mx),qys(mx,mx)
      common/bndry/ Un0(mx),Vn0(mx),Unp1(mx),Vnp1(mx),Pn0(mx),Pnp1(mx),
     >              Um0(mx),Vm0(mx),Ump1(mx),Vmp1(mx),Pm0(mx),Pmp1(mx),
     >              Fn0(mx),Gm0(mx)

! u-component
      DO i=2,N
      IF(i.le.NS) THEN
      	do j=MS+1,M+1
      	u(i,j) = F(i,j) + dt*(p(i,j)-p(i+1,j))/(rho*dx)
      	enddo
      ELSE
      	do j=2,M+1
      	u(i,j) = F(i,j) + dt*(p(i,j)-p(i+1,j))/(rho*dx)
      	enddo
      ENDIF
      ENDDO
! v-component
      DO i=2,N+1
      IF(i.lt.NS) THEN
      	do j=MS+1,M
      	v(i,j) = G(i,j) + dt*(p(i,j)-p(i,j+1))/(rho*dy)
      	enddo
      ELSE
      	do j=2,M
      	v(i,j) = G(i,j) + dt*(p(i,j)-p(i,j+1))/(rho*dy)
      	enddo
      ENDIF
      ENDDO

      return
      end
c-----------------------------------------------------------------------
	subroutine Output
	parameter(mx=1001)
	common/init/  rL,H,SL,SH,Um,vis,rho,gamma,dt
	common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >		  NS,MS
	common/flow/  u(mx,mx),v(mx,mx),p(mx,mx),S(mx,mx)
	common/xcoef/ ue(mx,mx),uw(mx,mx),un(mx,mx),us(mx,mx),conu(mx,mx),
     >              difu(mx,mx),F(mx,mx),gx(mx,mx),qxe(mx,mx),qxw(mx,mx)
     >              ,qxn(mx,mx),qxs(mx,mx)
	common/ycoef/ ve(mx,mx),vw(mx,mx),vn(mx,mx),vs(mx,mx),conv(mx,mx),
     >		  difv(mx,mx),G(mx,mx),gy(mx,mx),qye(mx,mx),qyw(mx,mx)
     >              ,qyn(mx,mx),qys(mx,mx)
      common/bndry/ Un0(mx),Vn0(mx),Unp1(mx),Vnp1(mx),Pn0(mx),Pnp1(mx),
     >              Um0(mx),Vm0(mx),Ump1(mx),Vmp1(mx),Pm0(mx),Pmp1(mx),
     >              Fn0(mx),Gm0(mx)
     	real disc_in,disc_out,u_out(mx,mx),v_out(mx,mx)
     	real stream(mx,mx)


     	do i=2,N+1
     	do j=2,M+1
     		u_out(i,j) = (u(i,j)+u(i-1,j))/2
     		v_out(i,j) = (v(i,j)+v(i-1,j))/2
     	enddo
     	enddo

     	do i=2,N+1
     		stream(i,1) = 0.
     		stream(i,2) = -u_out(i,j)*dy/2 + stream(i,1)
     	enddo

     	do i=2,N+1
     	do j=2,M
     		stream(i,j+1) = -u_out(i,j+1)*dy + stream(i,j)
     	enddo
     	enddo


      open(12,file="Velocities.tec",form='formatted')
      write(12,*) ' variables="x","y","u","v","p","stream" '
      write(12,*) 'zone i=',N, ' j=',M
      do j=2,M+1
      do i=2,N+1
      if((i.le.NS.and.j.le.MS).or.j.eq.1) then
      write(12,'(8E12.4)') xc(i-1,j-1),yc(i-1,j-1),0.,0.,p(i,j),0.
      else
      write(12,'(8E12.4)') xc(i-1,j-1),yc(i-1,j-1),
     + 			   u_out(i,j),v_out(i,j),p(i,j),stream(i,j)
      endif
      enddo
      enddo
      close(12)
      disc_in=0.
      disc_out=0.

      do j=MS+1,M+1
      	disc_in  = disc_in+ u(1,j)*dy
      enddo
      do j=2,M+1
      	disc_out = disc_out + u(N,j)*dy
      enddo
      print*, "Discharge In", disc_in, "Discharge Out", disc_out
      print*, "Percent Difference" ,abs(100*(disc_out-disc_in)/disc_in)



     	return
     	end
c-----------------------------------------------------------------------
	subroutine Error(iter)
	parameter(mx=1001)
	common/init/  rL,H,SL,SH,Um,vis,rho,gamma,dt
	common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >		  NS,MS
	common/flow/  u(mx,mx),v(mx,mx),p(mx,mx),S(mx,mx)
	common/xcoef/ ue(mx,mx),uw(mx,mx),un(mx,mx),us(mx,mx),conu(mx,mx),
     >              difu(mx,mx),F(mx,mx),gx(mx,mx),qxe(mx,mx),qxw(mx,mx)
     >              ,qxn(mx,mx),qxs(mx,mx)
	common/ycoef/ ve(mx,mx),vw(mx,mx),vn(mx,mx),vs(mx,mx),conv(mx,mx),
     >		  difv(mx,mx),G(mx,mx),gy(mx,mx),qye(mx,mx),qyw(mx,mx)
     >              ,qyn(mx,mx),qys(mx,mx)
      common/bndry/ Un0(mx),Vn0(mx),Unp1(mx),Vnp1(mx),Pn0(mx),Pnp1(mx),
     >              Um0(mx),Vm0(mx),Ump1(mx),Vmp1(mx),Pm0(mx),Pmp1(mx),
     >              Fn0(mx),Gm0(mx)
     	common/Err/   E,res,tolerance
      real u_old(mx,mx),v_old(mx,mx), sum


      if(iter.eq.1) then
      	sum = N*M
      	do i=2,N
      	do j=2,M
      		u_old(i,j) = u(i,j)
      		v_old(i,j) = v(i,j)
      	enddo
      	enddo
      else
      	sum = 0. 
      	do i=2,N
      	do j=2,M
      		sum = sum + abs(u_old(i,j)-u(i,j))
      		sum = sum + abs(v_old(i,j)-v(i,j))
      		u_old(i,j) = u(i,j)
      		v_old(i,j) = v(i,j)
      	enddo
      	enddo
      endif
      E = sum/(N*M*Um)
      write(13,*) iter, E , res/(N*M)
c      print*, 'velocity Res', E


     	return
     	end

c-----------------------------------------------------------------------
