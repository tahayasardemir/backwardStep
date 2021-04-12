      Program FlowOverBackwardFacingStep
c..Taha Yaşar Demir / 1881978
c..CE-580 / Homework #9
      parameter(mx=1001)
      common/init/  rL,H,SL,SH,Um,vis,rho,gamma,dt
      common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >              NS,MS
      common/flow/  u(mx,mx),v(mx,mx),p(mx,mx),S(mx,mx)
      common/xcoef/ ue(mx,mx),uw(mx,mx),un(mx,mx),us(mx,mx),conu(mx,mx),
     >              difu(mx,mx),F(mx,mx),gx(mx,mx),qxe(mx,mx),qxw(mx,mx)
     >              ,qxn(mx,mx),qxs(mx,mx)
      common/ycoef/ ve(mx,mx),vw(mx,mx),vn(mx,mx),vs(mx,mx),conv(mx,mx),
     >              difv(mx,mx),G(mx,mx),gy(mx,mx),qye(mx,mx),qyw(mx,mx)
     >              ,qyn(mx,mx),qys(mx,mx)
      common/bndry/ Un0(mx),Vn0(mx),Unp1(mx),Vnp1(mx),Pn0(mx),Pnp1(mx),
     >              Um0(mx),Vm0(mx),Ump1(mx),Vmp1(mx),Pm0(mx),Pmp1(mx),
     >              Fn0(mx),Gm0(mx)

      open(11,file='mesh.tec',form="formatted")
      call Initialize
      do j=1,2000
            call Boundary
            call Coefficients
            call Poisson
            call UpdateVelocities
      print*, j ,"----------------------------------------------------"
      enddo
      call Output

      stop
      end
c-----------------------------------------------------------------------
      subroutine Initialize
      parameter(mx=1001)
      common/init/  rL,H,SL,SH,Um,vis,rho,gamma,dt
      common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >              NS,MS

      rL  = 2.  !m /total length
      H   = 0.1 !m /total height 
      SL  = 0.4 !m /step length
      SH  = 0.05!m /step height
      Um  = 2.  !m/s 0.2--1.0--2.0
      vis = 1e-4!m^2/s 
      rho = 1e3 !km/m^3

      gamma = 1. ! Over-Relaxation Parameter

      call GenerateGrid

      dt = 1e-7

      call InitialCondition

      return
      end
c-----------------------------------------------------------------------
      subroutine GenerateGrid
      parameter(mx=1001)
      common/init/  rL,H,SL,SH,Um,vis,rho,gamma,dt
      common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >              NS,MS

      N  = 300 ! Cell count in x-direction
      M  = 60  ! Cell count in y-direction
      dx = 2.0/N 
      dy = 0.1/M
      NS = ceiling(0.4/dx) + 1
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

      write(11,*) ' variables="x","y" '
      write(11,*) ' zone i=',N+1, 'j=',M+1
      do j = 1,M+1
      do i = 1,N+1
            write(11,*) x(i,j),y(i,j)
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine InitialCondition
      parameter(mx=1001)
      common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >              NS,MS
      common/bndry/ Un0(mx),Vn0(mx),Unp1(mx),Vnp1(mx),Pn0(mx),Pnp1(mx),
     >              Um0(mx),Vm0(mx),Ump1(mx),Vmp1(mx),Pm0(mx),Pmp1(mx),
     >              Fn0(mx),Gm0(mx)
      common/flow/  u(mx,mx),v(mx,mx),p(mx,mx),S(mx,mx)
      real uniflow(30)

      do i=1,N
            do j=1,M
                  u(i,j) = 0.
                  v(i,j) = 0.
                  p(i,j) = 0.
            enddo
      enddo
      k=1
      do j = MS,M
            uniflow(k) = (-(320*yc(1,k)**2)+(16*yc(1,k)))/0.1
            Un0(j) = uniflow(k) ! wont updated in future iterations
            Vn0(j) = 0. ! wont updated in future iterations
            k = k+1
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine Boundary
      parameter(mx=1001)
      common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >              NS,MS
      common/flow/  u(mx,mx),v(mx,mx),p(mx,mx),S(mx,mx)
      common/xcoef/ ue(mx,mx),uw(mx,mx),un(mx,mx),us(mx,mx),conu(mx,mx),
     >              difu(mx,mx),F(mx,mx),gx(mx,mx),qxe(mx,mx),qxw(mx,mx)
     >              ,qxn(mx,mx),qxs(mx,mx)
      common/ycoef/ ve(mx,mx),vw(mx,mx),vn(mx,mx),vs(mx,mx),conv(mx,mx),
     >              difv(mx,mx),G(mx,mx),gy(mx,mx),qye(mx,mx),qyw(mx,mx)
     >              ,qyn(mx,mx),qys(mx,mx)
      common/bndry/ Un0(mx),Vn0(mx),Unp1(mx),Vnp1(mx),Pn0(mx),Pnp1(mx),
     >              Um0(mx),Vm0(mx),Ump1(mx),Vmp1(mx),Pm0(mx),Pmp1(mx),
     >              Fn0(mx),Gm0(mx)

      ! Inflow Boundary Conditions
      do j = 1,M
            Pn0(j) = p(1,j)
            Fn0(j) = Un0(j)
      enddo
      ! Outflow Boundary Conditions
      do j = 1,M
            Vnp1(j) = 0.
c           Unp1(j) = u(N,j)
            u(N,j)  = u(N-1,j)
            F(N,j)  = u(N,j)
            Pnp1(j) = p(N,j)
      enddo
      !! Wall Boundaries !!
      ! Horizaontal Part of the Step
      do i=1,NS
            v(i,MS-1) = 0.
            u(i,MS-1) =-u(i,MS)
            G(i,MS-1) = v(i,MS-1)
            p(i,MS-1) = p(i,MS)
      enddo
      ! Vertical Part of the Step
      do j=1,MS-1
            v(NS-1,j) =-v(NS,j)
            u(NS-1,j) = 0.
            F(NS-1,j) = u(NS-1,j)
            p(NS-1,j) = p(NS,j)
      enddo
      ! Lower Wall
      do i=NS,N
            Vm0(i) = 0.
            Um0(i) =-u(i,1)
            Gm0(i) = Vm0(i)
            Pm0(i) = p(i,1)
      enddo
      ! Upper Wall
      do i=1,N
            v(i,M) = 0.
            Ump1(i)=-u(i,M)
            G(i,M) = v(i,M)
            Pmp1(i)= p(i,M)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine Coefficients
      parameter(mx=1001)
      common/init/  rL,H,SL,SH,Um,vis,rho,gamma,dt
      common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >              NS,MS
      common/flow/  u(mx,mx),v(mx,mx),p(mx,mx),S(mx,mx)
      common/xcoef/ ue(mx,mx),uw(mx,mx),un(mx,mx),us(mx,mx),conu(mx,mx),
     >              difu(mx,mx),F(mx,mx),gx(mx,mx),qxe(mx,mx),qxw(mx,mx)
     >              ,qxn(mx,mx),qxs(mx,mx)
      common/ycoef/ ve(mx,mx),vw(mx,mx),vn(mx,mx),vs(mx,mx),conv(mx,mx),
     >              difv(mx,mx),G(mx,mx),gy(mx,mx),qye(mx,mx),qyw(mx,mx)
     >              ,qyn(mx,mx),qys(mx,mx)
      common/bndry/ Un0(mx),Vn0(mx),Unp1(mx),Vnp1(mx),Pn0(mx),Pnp1(mx),
     >              Um0(mx),Vm0(mx),Ump1(mx),Vmp1(mx),Pm0(mx),Pmp1(mx),
     >              Fn0(mx),Gm0(mx)


      do i=1,N
            do j=1,M
                  if(u(i,j).ge.0.) then
                  gx(i,j) = 1.
                  else
                  gx(i,j) =-1.
                  endif
                  if(v(i,j).ge.0.) then
                  gy(i,j) = 1.
                  else
                  gy(i,j) =-1.
                  endif
            enddo
      enddo 

      DO i =1,N
            if(i.le.NS) then ! Before Step
            do j=MS,M
                  if(j.lt.M.and.i.eq.1) then ! Inflow
                  ! u-convection constants
                        qxe(i,j) = 2*(u(i,j)+u(i+1,j))*dy
                        qxw(i,j) = 2*(Un0(j)+u(i,j))*dy
                        qxn(i,j) = 2*(v(i,j)+v(i+1,j))*dx
                        qxs(i,j) = 2*(v(i,j-1)+v(i+1,j-1))*dx
                  ue(i,j)= (u(i,j)+u(i+1,j)+gx(i,j)*sign(1.,qxe(i,j))
     +                    *(u(i,j)-u(i+1,j)))/2.
                  un(i,j)= (u(i,j)+u(i,j+1)+gx(i,j)*sign(1.,qxn(i,j))
     +                    *(u(i,j)-u(i,j+1)))/2.
                  uw(i,j)= (Un0(j)+u(i,j)+gx(i,j)*sign(1.,qxw(i,j))
     +                    *(Un0(j)-u(i,j)))/2.
                  us(i,j)= (u(i,j-1)+u(i,j)+gx(i,j)*sign(1.,qxs(i,j))
     +                    *(u(i,j-1)-u(i,j)))/2.
                  ! v-convection constants
                        qye(i,j) = 2*(u(i,j)+u(i,j+1))*dy
                        qyw(i,j) = 2*(Un0(j)+Un0(j+1))*dy
                        qyn(i,j) = 2*(v(i,j)+v(i,j+1))*dx
                        qys(i,j) = 2*(v(i,j)+v(i,j-1))*dx
                  ve(i,j)= (v(i,j)+v(i+1,j)+gy(i,j)*sign(1.,qye(i,j))
     +                    *(v(i,j)-v(i+1,j)))/2.
                  vw(i,j)= (Vn0(j)+v(i,j)+gy(i,j)*sign(1.,qyw(i,j))
     +                    *(Vn0(j)-v(i,j)))/2.
                  vn(i,j)= (v(i,j)+v(i,j+1)+gy(i,j)*sign(1.,qyn(i,j))
     +                    *(v(i,j)-v(i,j+1)))/2.
                  vs(i,j)= (v(i,j-1)+v(i,j)*gy(i,j)*sign(1.,qys(i,j))
     +                    *(v(i,j-1)-v(i,j)))/2.
                  ! u-diffusion constant
                  difu(i,j) = vis*( ( ((u(i+1,j)-u(i,j))/dx)
     +                              - ((u(i,j)-Un0(j))/dx) )*dy
     +                             +( ((u(i,j+1)-u(i,j))/dy)
     +                              - ((u(i,j)-u(i,j-1))/dy) )*dx )
                  ! v-diffusion constant
                  difv(i,j) = vis*( ( ((v(i+1,j)-v(i,j))/dx)
     +                              - ((v(i,j)-Vn0(j))/dx) )*dy
     +                             +( ((v(i,j+1)-v(i,j))/dy)
     +                              - ((v(i,j)-v(i,j-1))/dy) )*dx )                       

                  elseif(j.eq.M.and.i.gt.1) then ! Upper Wall Before Step
                  ! u-convection constants
                        qxe(i,j)=2*(u(i,j)+u(i+1,j))*dy
                        qxw(i,j)=2*(u(i-1,j)+u(i,j))*dy
                        qxn(i,j)=2*(v(i,j)+v(i+1,j))*dx
                        qxs(i,j)=2*(v(i,j-1)+v(i+1,j-1))*dx
                  ue(i,j)= (u(i,j)+u(i+1,j)+gx(i,j)*sign(1.,qxe(i,j))
     +                    *(u(i,j)-u(i+1,j)))/2.
                  un(i,j)= (u(i,j)+Ump1(i)+gx(i,j)*sign(1.,qxn(i,j))
     +                    *(u(i,j)-Ump1(i)))/2.
                  uw(i,j)= (u(i-1,j)+u(i,j)+gx(i,j)*sign(1.,qxw(i,j))
     +                    *(u(i-1,j)-u(i,j)))/2.
                  us(i,j)= (u(i,j-1)+u(i,j)+gx(i,j)*sign(1.,qxs(i,j))
     +                    *(u(i,j-1)-u(i,j)))/2.
                  ! v-convection constants
                        qye(i,j) = 2*(u(i,j)+Ump1(i))*dy
                        qyw(i,j) = 2*(u(i-1,j)+Ump1(i-1))*dy
                        qyn(i,j) = 2*(v(i,j)+Vmp1(i))*dx
                        qys(i,j) = 2*(v(i,j)+v(i,j-1))*dx
                  ve(i,j)= (v(i,j)+v(i+1,j)+gy(i,j)*sign(1.,qye(i,j))
     +                    *(v(i,j)-v(i+1,j)))/2.
                  vw(i,j)= (v(i-1,j)+v(i,j)+gy(i,j)*sign(1.,qyw(i,j))
     +                    *(v(i-1,j)-v(i,j)))/2.
                  vn(i,j)= (v(i,j)+Vmp1(i)+gy(i,j)*sign(1.,qyn(i,j))
     +                    *(v(i,j)-Vmp1(i)))/2.
                  vs(i,j)= (v(i,j-1)+v(i,j)*gy(i,j)*sign(1.,qys(i,j))
     +                    *(v(i,j-1)-v(i,j)))/2.
                  ! u-diffusion constant
                  difu(i,j) = vis*( ( ((u(i+1,j)-u(i,j))/dx)
     +                              - ((u(i,j)-u(i-1,j))/dx) )*dy
     +                             +( ((Ump1(i)-u(i,j))/dy)
     +                              - ((u(i,j)-u(i,j-1))/dy) )*dx )
                  ! v-diffusion constant
                  difv(i,j) = vis*( ( ((v(i+1,j)-v(i,j))/dx)
     +                              - ((v(i,j)-v(i-1,j))/dx) )*dy
     +                             +( ((Vmp1(i)-v(i,j))/dy)
     +                              - ((v(i,j)-v(i,j-1))/dy) )*dx )
                  elseif(i.eq.1.and.j.eq.M) then ! Upper Left Corner
                  ! u-convection constants
                        qxe(i,j) = 2*(u(i,j)+u(i+1,j))*dy
                        qxw(i,j) = 2*(Un0(j)+u(i,j))*dy
                        qxn(i,j) = 2*(v(i,j)+v(i+1,j))*dx
                        qxs(i,j) = 2*(v(i,j-1)+v(i+1,j-1))*dx
                  ue(i,j)= (u(i,j)+u(i+1,j)+gx(i,j)*sign(1.,qxe(i,j))
     +                    *(u(i,j)-u(i+1,j)))/2.
                  un(i,j)= (u(i,j)+Ump1(i)+gx(i,j)*sign(1.,qxn(i,j))
     +                    *(u(i,j)-Ump1(i)))/2.
                  uw(i,j)= (Un0(j)+u(i,j)+gx(i,j)*sign(1.,qxw(i,j))
     +                    *(Un0(j)-u(i,j)))/2.
                  us(i,j)= (u(i,j-1)+u(i,j)+gx(i,j)*sign(1.,qxs(i,j))
     +                    *(u(i,j-1)-u(i,j)))/2.
                  !print*, us(i,j),(u(i,j-1)),u(i,j),qxs(i,j) !---------------------------------------------
                  ! v-convection constants
                        qye(i,j) = 2*(u(i,j)+Ump1(i))*dy
                        qyw(i,j) = 2*(Un0(j)+Ump1(i))*dy 
                        qyn(i,j) = 2*(v(i,j)+Vmp1(i))*dx
                        qys(i,j) = 2*(v(i,j)+v(i,j-1))*dx
                  ve(i,j) = (v(i,j)+v(i+1,j)+gy(i,j)*sign(1.,qye(i,j))
     +                     *(v(i,j)-v(i+1,j)))/2.
                  vw(i,j) = (Vn0(j)+v(i,j)+gy(i,j)*sign(1.,qyw(i,j))
     +                     *(Vn0(j)-v(i,j)))/2.
                  vn(i,j) = (v(i,j)+Vmp1(i)+gy(i,j)*sign(1.,qyn(i,j))
     +                     *(v(i,j)-Vmp1(i)))/2.
                  vs(i,j) = (v(i,j-1)+v(i,j)+gy(i,j)*sign(1.,qys(i,j))
     +                     *(v(i,j-1)-v(i,j)))/2.
                  ! u-diffusion constant
                  difu(i,j) = vis*( ( ((u(i+1,j)-u(i,j))/dx)
     +                              - ((u(i,j)-Un0(j))/dx) )*dy
     +                             +( ((Ump1(i)-u(i,j))/dy)
     +                              - ((u(i,j)-u(i,j-1))/dy) )*dx )
                  ! v-diffusion constant
                  difv(i,j) = vis*( ( ((v(i+1,j)-v(i,j))/dx)
     +                              - ((v(i,j)-Vn0(j))/dx) )*dy
     +                             +( ((Vmp1(i)-v(i,j))/dy)
     +                              - ((v(i,j)-v(i,j-1))/dy) )*dx )                 
                  else ! Everywhere else before step
                  ! u-convection constants
                        qxe(i,j) = 2*(u(i,j)+u(i+1,j))*dy
                        qxw(i,j) = 2*(u(i-1,j)+u(i,j))*dy
                        qxn(i,j) = 2*(v(i,j)+v(i+1,j))*dx
                        qxs(i,j) = 2*(v(i,j-1)+v(i+1,j-1))*dx
                  ue(i,j)= (u(i,j)+u(i+1,j)+gx(i,j)*sign(1.,qxe(i,j))
     +                    *(u(i,j)-u(i+1,j)))/2.
                  un(i,j)= (u(i,j)+u(i,j+1)+gx(i,j)*sign(1.,qxn(i,j))
     +                    *(u(i,j)-u(i,j+1)))/2.
                  uw(i,j)= (u(i-1,j)+u(i,j)+gx(i,j)*sign(1.,qxw(i,j))
     +                    *(u(i-1,j)-u(i,j)))/2.
                  us(i,j)= (u(i,j-1)+u(i,j)+gx(i,j)*sign(1.,qxs(i,j))
     +                    *(u(i,j-1)-u(i,j)))/2.
                  ! v-convection constants
                        qye(i,j) = 2*(u(i,j)+u(i,j+1))*dy
                        qyw(i,j) = 2*(u(i-1,j)+u(i-1,j+1))*dy
                        qyn(i,j) = 2*(v(i,j)+v(i,j+1))*dx
                        qys(i,j) = 2*(v(i,j)+v(i,j-1))*dx
                  ve(i,j) = (v(i,j)+v(i+1,j)+gy(i,j)*sign(1.,qye(i,j))
     +                     *(v(i,j)-v(i+1,j)))/2.
                  vw(i,j) = (v(i-1,j)+v(i,j)+gy(i,j)*sign(1.,qyw(i,j))
     +                     *(v(i-1,j)-v(i,j)))/2.
                  vn(i,j) = (v(i,j)+v(i,j+1)+gy(i,j)*sign(1.,qyn(i,j))
     +                     *(v(i,j)-v(i,j+1)))/2.
                  vs(i,j) = (v(i,j-1)+v(i,j)+gy(i,j)*sign(1.,qys(i,j))
     +                     *(v(i,j-1)-v(i,j)))/2.
                  ! u-diffusion constant
                  difu(i,j) = vis*( ( ((u(i+1,j)-u(i,j))/dx)
     +                              - ((u(i,j)-u(i-1,j))/dx) )*dy
     +                             +( ((u(i,j+1)-u(i,j))/dy)
     +                              - ((u(i,j)-u(i,j-1))/dy) )*dx )
                  ! v-diffusion constant
                  difv(i,j) = vis*( ( ((v(i+1,j)-v(i,j))/dx)
     +                              - ((v(i,j)-v(i-1,j))/dx) )*dy
     +                             +( ((v(i,j+1)-v(i,j))/dy)
     +                              - ((v(i,j)-v(i,j-1))/dy) )*dx )
                  endif
            conu(i,j) = ue(i,j)*qxe(i,j)-uw(i,j)*qxw(i,j)+
     +                  un(i,j)*qxn(i,j)-us(i,j)*qxs(i,j)
            conv(i,j) = ve(i,j)*qye(i,j)-vw(i,j)*qyw(i,j)+
     +                  vn(i,j)*qyn(i,j)-vs(i,j)*qys(i,j)
            enddo
            else ! After Step
            do j=1,M
                  if(j.eq.1.and.i.lt.N) then ! Bottom Wall After Step
                  ! u-convection constants
                        qxe(i,j)=2*(u(i,j)+u(i+1,j))*dy
                        qxw(i,j)=2*(u(i-1,j)+u(i,j))*dy
                        qxn(i,j)=2*(v(i,j)+v(i+1,j))*dx
                        qxs(i,j)=2*(Vm0(i)+Vm0(i+1))*dx
                  ue(i,j)= (u(i,j)+u(i+1,j)+gx(i,j)*sign(1.,qxe(i,j))
     +                    *(u(i,j)-u(i+1,j)))/2.
                  un(i,j)= (u(i,j)+u(i,j+1)+gx(i,j)*sign(1.,qxn(i,j))
     +                    *(u(i,j)-u(i,j+1)))/2.
                  uw(i,j)= (u(i-1,j)+u(i,j)+gx(i,j)*sign(1.,qxw(i,j))
     +                    *(u(i-1,j)-u(i,j)))/2.
                  us(i,j)= (Um0(i)+u(i,j)+gx(i,j)*sign(1.,qxs(i,j))
     +                    *(Um0(i)-u(i,j)))/2.  
                  ! v-convection constants
                        qye(i,j) = 2*(u(i,j)+u(i,j+1))*dy
                        qyw(i,j) = 2*(u(i-1,j)+u(i-1,j+1))*dy
                        qyn(i,j) = 2*(v(i,j)+v(i,j+1))*dx
                        qys(i,j) = 2*(v(i,j)+Vm0(i))*dx
                  ve(i,j) = (v(i,j)+v(i+1,j)+gy(i,j)*sign(1.,qye(i,j))
     +                     *(v(i,j)-v(i+1,j)))/2.
                  vw(i,j) = (v(i-1,j)+v(i,j)+gy(i,j)*sign(1.,qyw(i,j))
     +                     *(v(i-1,j)-v(i,j)))/2.
                  vn(i,j) = (v(i,j)+v(i,j+1)+gy(i,j)*sign(1.,qyn(i,j))
     +                     *(v(i,j)-v(i,j+1)))/2.
                  vs(i,j) = (Vm0(i)+v(i,j)+gy(i,j)*sign(1.,qys(i,j))
     +                     *(Vm0(i)-v(i,j)))/2.
                  ! u-diffusion constant
                  difu(i,j) = vis*( ( ((u(i+1,j)-u(i,j))/dx)
     +                              - ((u(i,j)-u(i-1,j))/dx) )*dy
     +                             +( ((u(i,j+1)-u(i,j))/dy)
     +                              - ((u(i,j)-Um0(i))/dy) )*dx )
                  ! v-diffusion constant
                  difv(i,j) = vis*( ( ((v(i+1,j)-v(i,j))/dx)
     +                              - ((v(i,j)-v(i-1,j))/dx) )*dy
     +                             +( ((v(i,j+1)-v(i,j))/dy)
     +                              - ((v(i,j)-Vm0(i))/dy) )*dx )                               
                  elseif(j.gt.1.and.i.eq.N) then ! Outflow
                  ! u-convection constants
                        qxe(i,j)=2*(u(i,j)+Unp1(j))*dy
                        qxw(i,j)=2*(u(i-1,j)+u(i,j))*dy
                        qxn(i,j)=2*(v(i,j)+Vnp1(j))*dx
                        qxs(i,j)=2*(v(i,j-1)+Vnp1(j-1))*dx
                  ue(i,j)= (u(i,j)+Unp1(j)+gx(i,j)*sign(1.,qxe(i,j))
     +                    *(u(i,j)-Unp1(j)))/2.
                  un(i,j)= (u(i,j)+u(i,j+1)+gx(i,j)*sign(1.,qxn(i,j))
     +                    *(u(i,j)-u(i,j+1)))/2.
                  uw(i,j)= (u(i-1,j)+u(i,j)+gx(i,j)*sign(1.,qxw(i,j))
     +                    *(u(i-1,j)-u(i,j)))/2.
                  us(i,j)= (u(i,j-1)+u(i,j)+gx(i,j)*sign(1.,qxs(i,j))
     +                    *(u(i,j-1)-u(i,j)))/2.
                  ! v-convection constants
                        qye(i,j) = 2*(u(i,j)+u(i,j+1))*dy
                        qyw(i,j) = 2*(u(i-1,j)+u(i-1,j+1))*dy
                        qyn(i,j) = 2*(v(i,j)+v(i,j+1))*dx
                        qys(i,j) = 2*(v(i,j)+v(i,j-1))*dx
                  ve(i,j) = (v(i,j)+Vnp1(j)+gy(i,j)*sign(1.,qye(i,j))
     +                     *(v(i,j)-Vnp1(j)))/2.
                  vw(i,j) = (v(i-1,j)+v(i,j)+gy(i,j)*sign(1.,qyw(i,j))
     +                     *(v(i-1,j)-v(i,j)))/2.
                  vn(i,j) = (v(i,j)+v(i,j+1)+gy(i,j)*sign(1.,qyn(i,j))
     +                     *(v(i,j)-v(i,j+1)))/2.
                  vs(i,j) = (v(i,j-1)+v(i,j)+gy(i,j)*sign(1.,qys(i,j))
     +                     *(v(i,j-1)-v(i,j)))/2.
                  ! u-diffusion constant
                  difu(i,j) = vis*( ( ((Ump1(j)-u(i,j))/dx)
     +                              - ((u(i,j)-u(i-1,j))/dx) )*dy
     +                             +( ((u(i,j+1)-u(i,j))/dy)
     +                              - ((u(i,j)-u(i,j-1))/dy) )*dx )
                  ! v-diffusion constant
                  difv(i,j) = vis*( ( ((Vnp1(j)-v(i,j))/dx)
     +                              - ((v(i,j)-v(i-1,j))/dx) )*dy
     +                             +( ((v(i,j+1)-v(i,j))/dy)
     +                              - ((v(i,j)-v(i,j-1))/dy) )*dx )
                  elseif(j.eq.M.and.i.lt.N) then ! Upper Wall After Step
                  ! u-convection constants
                        qxe(i,j)=2*(u(i,j)+u(i+1,j))*dy
                        qxw(i,j)=2*(u(i-1,j)+u(i,j))*dy
                        qxn(i,j)=2*(v(i,j)+v(i+1,j))*dx
                        qxs(i,j)=2*(v(i,j-1)+v(i+1,j-1))*dx
                  ue(i,j)= (u(i,j)+u(i+1,j)+gx(i,j)*sign(1.,qxe(i,j))
     +                    *(u(i,j)-u(i+1,j)))/2.
                  un(i,j)= (u(i,j)+Ump1(i)+gx(i,j)*sign(1.,qxn(i,j))
     +                    *(u(i,j)-Ump1(i)))/2.
                  uw(i,j)= (u(i-1,j)+u(i,j)+gx(i,j)*sign(1.,qxw(i,j))
     +                    *(u(i-1,j)-u(i,j)))/2.
                  us(i,j)= (u(i,j-1)+u(i,j)+gx(i,j)*sign(1.,qxs(i,j))
     +                    *(u(i,j-1)-u(i,j)))/2.
                  ! v-convection constants
                        qye(i,j) = 2*(u(i,j)+Ump1(i))*dy
                        qyw(i,j) = 2*(u(i-1,j)+Ump1(i-1))*dy
                        qyn(i,j) = 2*(v(i,j)+Vmp1(i))*dx
                        qys(i,j) = 2*(v(i,j)+v(i,j-1))*dx
                  ve(i,j) = (v(i,j)+v(i+1,j)+gy(i,j)*sign(1.,qye(i,j))
     +                     *(v(i,j)-v(i+1,j)))/2.
                  vw(i,j) = (v(i-1,j)+v(i,j)+gy(i,j)*sign(1.,qyw(i,j))
     +                     *(v(i-1,j)-v(i,j)))/2.
                  vn(i,j) = (v(i,j)+Vmp1(i)+gy(i,j)*sign(1.,qyn(i,j))
     +                     *(v(i,j)-v(i,j+1)))/2.
                  vs(i,j) = (v(i,j-1)+v(i,j)+gy(i,j)*sign(1.,qys(i,j))
     +                     *(v(i,j-1)-v(i,j)))/2.
                  ! u-diffusion constant
                  difu(i,j) = vis*( ( ((u(i+1,j)-u(i,j))/dx)
     +                              - ((u(i,j)-u(i-1,j))/dx) )*dy
     +                             +( ((Ump1(i)-u(i,j))/dy)
     +                              - ((u(i,j)-u(i,j-1))/dy) )*dx )
                  ! v-diffusion constant
                  difv(i,j) = vis*( ( ((v(i+1,j)-v(i,j))/dx)
     +                              - ((v(i,j)-v(i-1,j))/dx) )*dy
     +                             +( ((Vmp1(i)-v(i,j))/dy)
     +                              - ((v(i,j)-v(i,j-1))/dy) )*dx )           
                  else ! Interior After Step
                        qxe(i,j)=2*(u(i,j)+u(i+1,j))*dy
                        qxw(i,j)=2*(u(i-1,j)+u(i,j))*dy
                        qxn(i,j)=2*(v(i,j)+v(i+1,j))*dx
                        qxs(i,j)=2*(v(i,j-1)+v(i+1,j-1))*dx
                  ue(i,j)= (u(i,j)+u(i+1,j)+gx(i,j)*sign(1.,qxe(i,j))
     +                    *(u(i,j)-u(i+1,j)))/2.
                  un(i,j)= (u(i,j)+u(i,j+1)+gx(i,j)*sign(1.,qxn(i,j))
     +                    *(u(i,j)-u(i,j+1)))/2.
                  uw(i,j)= (u(i-1,j)+u(i,j)+gx(i,j)*sign(1.,qxw(i,j))
     +                    *(u(i-1,j)-u(i,j)))/2.
                  us(i,j)= (u(i,j-1)+u(i,j)+gx(i,j)*sign(1.,qxs(i,j))
     +                    *(u(i,j-1)-u(i,j)))/2.
                        qye(i,j) = 2*(u(i,j)+u(i,j+1))*dy
                        qyw(i,j) = 2*(u(i-1,j)+u(i-1,j+1))*dy
                        qyn(i,j) = 2*(v(i,j)+v(i,j+1))*dx
                        qys(i,j) = 2*(v(i,j)+v(i,j-1))*dx
                  ve(i,j) = (v(i,j)+v(i+1,j)+gy(i,j)*sign(1.,qye(i,j))
     +                     *(v(i,j)-v(i+1,j)))/2.
                  vw(i,j) = (v(i-1,j)+v(i,j)+gy(i,j)*sign(1.,qyw(i,j))
     +                     *(v(i-1,j)-v(i,j)))/2.
                  vn(i,j) = (v(i,j)+v(i,j+1)+gy(i,j)*sign(1.,qyn(i,j))
     +                     *(v(i,j)-v(i,j+1)))/2.
                  vs(i,j) = (v(i,j-1)+v(i,j)+gy(i,j)*sign(1.,qys(i,j))
     +                     *(v(i,j-1)-v(i,j)))/2.
                  ! u-diffusion constant
                  difu(i,j) = vis*( ( ((u(i+1,j)-u(i,j))/dx)
     +                              - ((u(i,j)-u(i-1,j))/dx) )*dy
     +                             +( ((u(i,j+1)-u(i,j))/dy)
     +                              - ((u(i,j)-u(i,j-1))/dy) )*dx )
                  ! v-diffusion constant
                  difv(i,j) = vis*( ( ((v(i+1,j)-v(i,j))/dx)
     +                              - ((v(i,j)-v(i-1,j))/dx) )*dy
     +                             +( ((v(i,j+1)-v(i,j))/dy)
     +                              - ((v(i,j)-v(i,j-1))/dy) )*dx )
                  endif
            conu(i,j) = ue(i,j)*qxe(i,j)-uw(i,j)*qxw(i,j)+
     +                  un(i,j)*qxn(i,j)-us(i,j)*qxs(i,j)
            conv(i,j) = ve(i,j)*qye(i,j)-vw(i,j)*qyw(i,j)+
     +                  vn(i,j)*qyn(i,j)-vs(i,j)*qys(i,j)
            enddo
            endif
      ENDDO

      DO i=1,N
      IF (i.le.NS) THEN
            do j=MS,M-1
            F(i,j)    = u(i,j) + (dt/(dx*dy))*(difu(i,j)-conu(i,j))
            G(i,j)    = v(i,j) + (dt/(dx*dy))*(difv(i,j)-conv(i,j))
            enddo
      ELSEIF(i.gt.NS) THEN
            do j=1,M-1
            F(i,j)    = u(i,j) + (dt/(dx*dy))*(difu(i,j)-conu(i,j))
            G(i,j)    = v(i,j) + (dt/(dx*dy))*(difv(i,j)-conv(i,j))
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
     >              NS,MS
      common/flow/  u(mx,mx),v(mx,mx),p(mx,mx),S(mx,mx)
      common/xcoef/ ue(mx,mx),uw(mx,mx),un(mx,mx),us(mx,mx),conu(mx,mx),
     >              difu(mx,mx),F(mx,mx),gx(mx,mx),qxe(mx,mx),qxw(mx,mx)
     >              ,qxn(mx,mx),qxs(mx,mx)
      common/ycoef/ ve(mx,mx),vw(mx,mx),vn(mx,mx),vs(mx,mx),conv(mx,mx),
     >              difv(mx,mx),G(mx,mx),gy(mx,mx),qye(mx,mx),qyw(mx,mx)
     >              ,qyn(mx,mx),qys(mx,mx)
      common/bndry/ Un0(mx),Vn0(mx),Unp1(mx),Vnp1(mx),Pn0(mx),Pnp1(mx),
     >              Um0(mx),Vm0(mx),Ump1(mx),Vmp1(mx),Pm0(mx),Pmp1(mx),
     >              Fn0(mx),Gm0(mx)

      ! Source Term Calculation
      DO i=1,N
      IF (i.le.NS) THEN ! Before Step
            do j=MS,M
            if (i.eq.1) then
            S(i,j) = rho*((F(i,j)-Fn0(j))/dx+(G(i,j)-G(i,j-1))/dy)/dt
c           print*, i,j,S(i,j),F(i,j),Fn0(j),G(i,j),G(i,j-1)
            else
            S(i,j) = rho*((F(i,j)-F(i-1,j))/dx+(G(i,j)-G(i,j-1))/dy)/dt
c           if(i.eq.2) print*, i,j,F(i,j)-F(i-1,j),G(i,j)-G(i,j-1)
            endif
            enddo
      ELSE ! After Step
            do j=1,M
            if(j.eq.1) then
            S(i,j) = rho*((F(i,j)-F(i-1,j))/dx+(G(i,j)-Gm0(i))/dy)/dt
            else
            S(i,j) = rho*((F(i,j)-F(i-1,j))/dx+(G(i,j)-G(i,j-1))/dy)/dt
            endif
            enddo
      ENDIF
      ENDDO

      ! Pressure Update 
      do k=1,10
            !print*, "--------------------------------------------------"
      DO i=1,N
      IF (i.le.NS) THEN
            do j=MS,M
            if (i.eq.1.and.j.lt.M) then
            p(i,j) = p(i,j)+ gamma*(((((p(i+1,j)+Pn0(j))/(dx**2) + 
     +              (p(i,j+1)+p(i,j-1))/(dy**2))-S(i,j)))/a - p(i,j))
            !print*, "Pn0", Pn0(j), p(i,j)
            elseif(i.eq.1.and.j.eq.M) then
            p(i,j) = p(i,j)+ gamma*(((((p(i+1,j)+Pn0(j))/(dx**2) + 
     +              (Pmp1(i)+p(i,j-1))/(dy**2))-S(i,j)))/a - p(i,j))
            elseif(i.gt.1.and.j.eq.M) then
            p(i,j) = p(i,j)+ gamma*(((((p(i+1,j)+p(i-1,j))/(dx**2) + 
     +              (Pmp1(i)+p(i,j-1))/(dy**2))-S(i,j)))/a - p(i,j))
            else
            !if (i.eq.2) print*, i,j,p(i,j),p(i,j+1),p(i,j-1),S(i,j)
            !if (i.eq.2) print*, p(i-1,j),p(i+1,j),gamma,a,dx**2,dy**2
            p(i,j) = p(i,j)+ gamma*(((((p(i+1,j)+p(i-1,j))/(dx**2) + 
     +              (p(i,j+1)+p(i,j-1))/(dy**2))-S(i,j)))/a - p(i,j))

            endif
            enddo
      ELSE
            do j=1,M
            if (i.lt.N.and.j.eq.1) then 
            p(i,j) = p(i,j)+ gamma*(((((p(i+1,j)+p(i-1,j))/dx**2 + 
     +              (p(i,j+1)+Pm0(i))/dy**2)-S(i,j)))/a - p(i,j))
            elseif (i.lt.N.and.j.eq.M) then
            p(i,j) = p(i,j)+ gamma*(((((p(i+1,j)+p(i-1,j))/dx**2 + 
     +              (Pmp1(i)+p(i,j-1))/dy**2)-S(i,j)))/a - p(i,j))
            elseif (i.eq.N) then
            p(i,j) = p(i,j)+ gamma*(((((Pnp1(j)+p(i-1,j))/dx**2 + 
     +              (p(i,j+1)+p(i,j-1))/dy**2)-S(i,j)))/a - p(i,j))
            else
            p(i,j) = p(i,j)+ gamma*(((((p(i+1,j)+p(i-1,j))/dx**2 + 
     +              (p(i,j+1)+p(i,j-1))/dy**2)-S(i,j)))/a - p(i,j))
            endif
            enddo
      ENDIF
      ENDDO
      enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine UpdateVelocities
      parameter(mx=1001)
      common/init/  rL,H,SL,SH,Um,vis,rho,gamma,dt
      common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >              NS,MS
      common/flow/  u(mx,mx),v(mx,mx),p(mx,mx),S(mx,mx)
      common/xcoef/ ue(mx,mx),uw(mx,mx),un(mx,mx),us(mx,mx),conu(mx,mx),
     >              difu(mx,mx),F(mx,mx),gx(mx,mx),qxe(mx,mx),qxw(mx,mx)
     >              ,qxn(mx,mx),qxs(mx,mx)
      common/ycoef/ ve(mx,mx),vw(mx,mx),vn(mx,mx),vs(mx,mx),conv(mx,mx),
     >              difv(mx,mx),G(mx,mx),gy(mx,mx),qye(mx,mx),qyw(mx,mx)
     >              ,qyn(mx,mx),qys(mx,mx)
      common/bndry/ Un0(mx),Vn0(mx),Unp1(mx),Vnp1(mx),Pn0(mx),Pnp1(mx),
     >              Um0(mx),Vm0(mx),Ump1(mx),Vmp1(mx),Pm0(mx),Pmp1(mx),
     >              Fn0(mx),Gm0(mx)
      real sumu,sumv

      sumu = 0.
      sumv = 0.

      DO i=1,N
      IF(i.le.NS) THEN
            do j=MS,M
                  if (i.eq.1.or.j.eq.M) then
                  u(i,j) = F(i,j) + dt*(p(i,j)-p(i+1,j))/(rho*dx)
                  v(i,j) = 0.
c                 print*, i,j,u(i,j),F(i,j),Un0(j)
                  else
                  u(i,j) = F(i,j) + dt*(p(i,j)-p(i+1,j))/(rho*dx)
                  v(i,j) = G(i,j) + dt*(p(i,j)-p(i,j+1))/(rho*dy)
c                 if(u(i,j).gt.8.) print*, i,j,p(i,j),p(i,j+1)
                  endif
            sumu = sumu + u(i,j)
            sumv = sumv + v(i,j)
            enddo
      ELSE
            do j=1,M
                  if (j.eq.M) then
                  u(i,j) = F(i,j) + dt*(p(i,j)-p(i+1,j))/(rho*dx)
                  v(i,j) = 0.
                  elseif (i.eq.N) then
                  u(i,j) = u(i-1,j)
                  v(i,j) = 0.
                  else
                  u(i,j) = F(i,j) + dt*(p(i,j)-p(i+1,j))/(rho*dx)
                  v(i,j) = G(i,j) + dt*(p(i,j)-p(i,j+1))/(rho*dy)
                  endif
            enddo
      ENDIF
      ENDDO

      print*, sumu,sumv

      return
      end
c-----------------------------------------------------------------------
      subroutine Output
      parameter(mx=1001)
      common/init/  rL,H,SL,SH,Um,vis,rho,gamma,dt
      common/grid/  N,M,dx,dy,x(mx,mx),xc(mx,mx),y(mx,mx),yc(mx,mx),a,
     >              NS,MS
      common/flow/  u(mx,mx),v(mx,mx),p(mx,mx),S(mx,mx)
      common/xcoef/ ue(mx,mx),uw(mx,mx),un(mx,mx),us(mx,mx),conu(mx,mx),
     >              difu(mx,mx),F(mx,mx),gx(mx,mx),qxe(mx,mx),qxw(mx,mx)
     >              ,qxn(mx,mx),qxs(mx,mx)
      common/ycoef/ ve(mx,mx),vw(mx,mx),vn(mx,mx),vs(mx,mx),conv(mx,mx),
     >              difv(mx,mx),G(mx,mx),gy(mx,mx),qye(mx,mx),qyw(mx,mx)
     >              ,qyn(mx,mx),qys(mx,mx)
      common/bndry/ Un0(mx),Vn0(mx),Unp1(mx),Vnp1(mx),Pn0(mx),Pnp1(mx),
     >              Um0(mx),Vm0(mx),Ump1(mx),Vmp1(mx),Pm0(mx),Pmp1(mx),
     >              Fn0(mx),Gm0(mx)

      open(11,file="Velocities.tec",form='formatted')
      write(11,*) ' variables="x","y","u","v","p" '
      write(11,*) 'zone i=',N, ' j=',M
      do j=1,M
            do i=1,N
            write(11,'(8E12.4)') x(i,j),y(i,j),u(i,j),v(i,j),p(i,j)
            enddo
      enddo
      close(11)


      return
      end
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
     
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------