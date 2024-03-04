subroutine funRrz(N_node,N_face,rho,N_GAUSS,Matrix_P0,F1,N_edge,C1,Fac_Ed_loc,Aed,Led,N_thread,R)
implicit none
integer*8 N_node,N_face,N_edge,N_thread,N_GAUSS
real*8 Matrix_P0(N_node,3)
integer*8 G2(8,N_edge), F1(4,N_face), C1(4,N_face), Fac_Ed_loc(4,N_face)
real*8 Aed(N_edge), Led(N_edge), rho(N_face)
real*8 R(N_edge,N_edge)
!!
real*8 pi
real*8 xabsc(N_GAUSS), weig(N_GAUSS)
integer*8 NP
real*8 xi(N_GAUSS**2), eta(N_GAUSS**2), zet(N_GAUSS**2), WG_loc(N_GAUSS**2)
integer*8 ii, kk, hh, qq, ee, jj
real*8 Hexa_master(8,3), PH1(3,N_GAUSS**2)
real*8 wf_loc(3,N_GAUSS**2,4)
integer*8 indii, ind_e1s, ind_e1, ind_loc_wk, ind_e2s, ind_e2
integer*8 ind_loc_wh
real*8 Face(4,3), Hexa(8,3), xyz(3)
real*8 PP(N_GAUSS**2,3), J(3,3), dett, detJ(N_GAUSS**2)
real*8 we_glo(3,N_GAUSS**2,4), tmp3(3), Rcoeff
!!
pi=4.d0*datan(1.d0)
call gauleg(N_GAUSS,xabsc(1:N_GAUSS),weig(1:N_GAUSS))
NP=N_GAUSS**2
xi(1:NP)=0.d0
eta(1:NP)=-1.d0
zet(1:NP)=0.d0
WG_loc(1:NP)=0.d0
hh=1;
do ii = 1,N_GAUSS 
    do kk = 1,N_GAUSS
           xi(hh) =xabsc(ii)
           zet(hh)=xabsc(kk)
           WG_loc(hh)= weig(ii)*weig(kk)
           hh=hh+1
    enddo 
enddo
!
Hexa_master(1,1:3)=[-1.d0,-1.d0, 1.d0] 
Hexa_master(2,1:3)=[-1.d0,-1.d0,-1.d0] 
Hexa_master(3,1:3)=[+1.d0,-1.d0,-1.d0]  
Hexa_master(4,1:3)=[+1.d0,-1.d0, 1.d0]  
Hexa_master(5,1:3)=[-1.d0,+1.d0, 1.d0]  
Hexa_master(6,1:3)=[-1.d0,+1.d0,-1.d0]  
Hexa_master(7,1:3)=[+1.d0,+1.d0,-1.d0]  
Hexa_master(8,1:3)=[+1.d0,+1.d0, 1.d0]
!!
PH1(1,1:NP)=xi
PH1(2,1:NP)=eta
PH1(3,1:NP)=zet
call fun_whitney_face_hexa_6_quad(Hexa_master,PH1,NP,wf_loc)
wf_loc=2.d0*wf_loc
!!
R(1:N_edge,1:N_edge)=0.d0
!!
call omp_set_num_threads(N_thread)
!!
do ii = 1,N_face
    Face=Matrix_P0(F1(1:4,ii),1:3)
    Hexa(1:4,[1,3])=Face(1:4,[1,3])
    Hexa(5:8,[1,3])=Face(1:4,[1,3])
    Hexa(1:4,2)=-1.d0
    Hexa(5:8,2)=1.d0
!!  punti di gauss in globale
	do kk = 1,NP
        call funTrilinear(xi(kk),0.0,zet(kk),Hexa,xyz)
		PP(kk,1:3)=xyz
    enddo
    !!
    !do jj = 1!:NP
	jj=1
        J(1:3,1)=(-(1.0d0-eta(jj))*(1.0d0+zet(jj))*Hexa(1,1:3)-(1.0d0-eta(jj))*(1.0d0-zet(jj))*Hexa(2,1:3) &
                  +(1.0d0-eta(jj))*(1.0d0-zet(jj))*Hexa(3,1:3)+(1.0d0-eta(jj))*(1.0d0+zet(jj))*Hexa(4,1:3) &
                  -(1.0d0+eta(jj))*(1.0d0+zet(jj))*Hexa(5,1:3)-(1.0d0+eta(jj))*(1.0d0-zet(jj))*Hexa(6,1:3) &
                  +(1.0d0+eta(jj))*(1.0d0-zet(jj))*Hexa(7,1:3)+(1.0d0+eta(jj))*(1.0d0+zet(jj))*Hexa(8,1:3))*0.125d0

        J(1:3,2)=(-(1.0d0- xi(jj))*(1.0d0+zet(jj))*Hexa(1,1:3)-(1.0d0- xi(jj))*(1.0d0-zet(jj))*Hexa(2,1:3) &
                  -(1.0d0+ xi(jj))*(1.0d0-zet(jj))*Hexa(3,1:3)-(1.0d0+ xi(jj))*(1.0d0+zet(jj))*Hexa(4,1:3) &
                  +(1.0d0- xi(jj))*(1.0d0+zet(jj))*Hexa(5,1:3)+(1.0d0- xi(jj))*(1.0d0-zet(jj))*Hexa(6,1:3) &
                  +(1.0d0+ xi(jj))*(1.0d0-zet(jj))*Hexa(7,1:3)+(1.0d0+ xi(jj))*(1.0d0+zet(jj))*Hexa(8,1:3))*0.125d0

        J(1:3,3)=(+(1.0d0- xi(jj))*(1.0d0-eta(jj))*Hexa(1,1:3)-(1.0d0- xi(jj))*(1.0d0-eta(jj))*Hexa(2,1:3) &
                  -(1.0d0+ xi(jj))*(1.0d0-eta(jj))*Hexa(3,1:3)+(1.0d0+ xi(jj))*(1.0d0-eta(jj))*Hexa(4,1:3) &
                  +(1.0d0- xi(jj))*(1.0d0+eta(jj))*Hexa(5,1:3)-(1.0d0- xi(jj))*(1.0d0+eta(jj))*Hexa(6,1:3) &
                  -(1.0d0+ xi(jj))*(1.0d0+eta(jj))*Hexa(7,1:3)+(1.0d0+ xi(jj))*(1.0d0+eta(jj))*Hexa(8,1:3))*0.125d0
        call det3(J,dett)
	    detJ(1:NP)=dett     
    !enddo
!     Vol=Vol*2*pi
!! 
    do jj=1,4
        do kk = 1,NP
		    call mat33_x_vec3(J(1:3,1:3),wf_loc(1:3,kk,jj),tmp3(1:3))
            we_glo(1:3,kk,jj)=(1.d0/detJ(kk))*tmp3(1:3)
        enddo
    enddo

    do kk = 1,4 
        ind_e1s=C1(kk,ii)
        ind_e1=abs(ind_e1s)
        ind_loc_wk=Fac_Ed_loc(kk,ii)
        do hh = 1,4
            ind_e2s=(C1(hh,ii))
            ind_e2=abs(ind_e2s)
            ind_loc_wh=Fac_Ed_loc(hh,ii)
            Rcoeff=0.0;
            do ee = 1,NP
                Rcoeff=Rcoeff+2.d0*pi*PP(ee,1)&
                              *rho(ii)&
                              *dot_product(we_glo(1:3,ee,ind_loc_wk)*Led(ind_e1)/Aed(ind_e1),we_glo(1:3,ee,ind_loc_wh)*Led(ind_e2)/Aed(ind_e2))&
                              *WG_loc(ee)*detJ(ee)&
                              *sign(1.d0,real(ind_e1s,8))*sign(1.d0,real(ind_e2s,8))
            enddo
            R(ind_e1,ind_e2)=R(ind_e1,ind_e2)+Rcoeff
        enddo
    enddo    
enddo
end subroutine funRrz
!!
subroutine mat33_x_vec3(a,b,c) 
implicit none
real(kind=8), dimension(3,3) :: a
real(kind=8), dimension(3) :: b
real(kind=8), dimension(3) :: c
c(1)=a(1,1)*b(1)+a(1,2)*b(2)+a(1,3)*b(3)
c(2)=a(2,1)*b(1)+a(2,2)*b(2)+a(2,3)*b(3)
c(3)=a(3,1)*b(1)+a(3,2)*b(2)+a(3,3)*b(3)
end subroutine mat33_x_vec3
!!
subroutine my_ellipke(m,k,e)
implicit none 
real*8 m, k, e, tol, pi, c1, a1, b1
real*8 a0, b0, c0, s0,  i1, mm, w1
!ELLIPKE Complete elliptic integral.
!   [K,E] = ELLIPKE(M) returns the value of the complete elliptic
!   integrals of the first and second kinds, evaluated do each
!   element of M.  As currently implemented, M is limited to 0 <= M <= 1.
!   
!   [K,E] = ELLIPKE(M,TOL) computes the complete elliptic integrals to
!   the accuracy TOL instead of the default TOL = EPS(CLASS(M)).
!
!   Some definitions of the complete elliptic integrals use the modulus
!   k instead of the parameter M.  They are related by M = k**2.
!
!   Class support do input M:
!      float: double, single
!
!   See also ELLIPJ.

!   Modified to include the second kind by Bjorn Bonnevier
!   from the Alfven Laboratory, KTH, Stockholm, Sweden
!   Copyright 1984-2013 The MathWorks, Inc. 

!   ELLIPKE uses the method of the arithmetic-geometric mean
!   described in [1].

!   References:
!   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
!       Functions" Dover Publications", 1965, 17.6.
tol=1.0d-20!2.2204d-16
pi=4.d0*datan(1.0d0)
a0 = 1.0d0
b0 = dsqrt(1.0d0-m)
c0 = 0.0d0
s0 = m
i1 = 0.0d0 
mm = 666.0d6
do while (mm>=tol) 
    a1 = (a0+b0)/2.0d0
    b1 = sqrt(a0*b0)
    c1 = (a0-b0)/2.0d0
    i1 = i1 + 1.0d0
    w1 = (2.0d0**i1)*(c1**2.0d0)
    mm = w1
    
    ! test do stagnation (may happen do TOL < machine precision)
    !if (c0==c1) then
	    !write(66,*) 'my_ellipke fails'
    !    stop 
    !end if
    
    s0 = s0 + w1  
    a0 = a1  
	b0 = b1  
	c0 = c1
end do
k = pi/(2.0*a1)
e = k*(1.0-s0/2.0)
end subroutine my_ellipke
!!
subroutine gauleg(ngp,xabsc,weig)
implicit none
integer i,j,m
real(kind=8) p1,p2,p3,pp,z,z1
integer, intent(IN) :: ngp !# of Gauss Points
real(kind=8), intent(OUT) :: xabsc(ngp),weig(ngp)
real(kind=8) :: eps2, pi
eps2=3.0d-15
pi=3.141592653589793d0
m=(ngp+1)/2
!Roots are symmetric in the interval so only need to find half of them
do i=1,m
z=cos(pi*(i-0.25d0)/(ngp+0.5d0)) !starting approximation
!Newton's method
100     p1 = 1.0d0
p2 = 0.0d0
!Loop up the recurrence relation to get the Legendre
!polynomial evaluated at z
do j=1,ngp
p3=p2
p2=p1
p1=((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
end do
!p1 is now the desired Legendre polynomial. We next compute pp,
!its derivative, by a standard relation involving also p2, the
!polynomial of one lower order.
pp=ngp*(z*p1-p2)/(z*z-1.0d0)
z1=z
z=z1-p1/pp             !Newton's Method
if (dabs(z-z1) .gt. EPS2) GOTO  100
xabsc(i)=-z                    		! Roots will be bewteen -1.0 & 1.0
xabsc(ngp+1-i)=+z             		! and symmetric about the origin
weig(i)=2.0d0/((1.0d0-z*z)*pp*pp)	! Compute the weight and its
weig(ngp+1-i)=weig(i)               ! symmetric counterpart
end do
end subroutine gauleg
!!
subroutine  fun_whitney_face_hexa_6_quad(Hexa,P,NP,wf)
implicit none
real*8 Hexa(8,3)
integer*8 NP
real*8 P(3,NP)
real*8 wf(3,NP,4), tmp3(3)
real*8 xi(NP), eta(NP), zet(NP)
real*8 grad_snF(6,3), J(3,3), sn(8), wf_loc(3,4), dett
real*8 gF64(3), gF45(3), gF52(3), gF26(3), gF61(3), gF15(3), gF53(3), gF36(3), gF62(3), gF25(3), gF41(3), gF32(3)
real*8 gF54(3), gF46(3), gF63(3), gF35(3), gF51(3), gF16(3), gF14(3), gF43(3), gF21(3), gF12(3), gF23(3), gF34(3)
!!
integer*8 ii, jj
xi(1:NP) =P(1,1:NP)
eta(1:NP)=P(2,1:NP)
zet(1:NP)=P(3,1:NP)
grad_snF(1,1:3)=[0.0d0,-0.5d0,0.0d0]
grad_snF(2,1:3)=[0.5d0,0.0d0,0.0d0]
grad_snF(3,1:3)=[0.0d0,0.5d0,0.0d0]
grad_snF(4,1:3)=[-0.5d0,0.0d0,0.0d0]
grad_snF(5,1:3)=[0.0d0,0.0d0,-0.5d0]
grad_snF(6,1:3)=[0.0d0,0.0d0,0.5d0]
!!
call fun_my_cross(grad_snF(6,1:3),grad_snF(4,1:3),gF64)
call fun_my_cross(grad_snF(4,1:3),grad_snF(5,1:3),gF45)
call fun_my_cross(grad_snF(5,1:3),grad_snF(2,1:3),gF52)
call fun_my_cross(grad_snF(2,1:3),grad_snF(6,1:3),gF26)
call fun_my_cross(grad_snF(6,1:3),grad_snF(1,1:3),gF61)
call fun_my_cross(grad_snF(1,1:3),grad_snF(5,1:3),gF15)
call fun_my_cross(grad_snF(5,1:3),grad_snF(3,1:3),gF53)
call fun_my_cross(grad_snF(3,1:3),grad_snF(6,1:3),gF36)
call fun_my_cross(grad_snF(6,1:3),grad_snF(2,1:3),gF62)
call fun_my_cross(grad_snF(2,1:3),grad_snF(5,1:3),gF25)
call fun_my_cross(grad_snF(5,1:3),grad_snF(4,1:3),gF54)
call fun_my_cross(grad_snF(4,1:3),grad_snF(6,1:3),gF46)
call fun_my_cross(grad_snF(6,1:3),grad_snF(3,1:3),gF63)
call fun_my_cross(grad_snF(3,1:3),grad_snF(5,1:3),gF35)
call fun_my_cross(grad_snF(5,1:3),grad_snF(1,1:3),gF51)
call fun_my_cross(grad_snF(1,1:3),grad_snF(6,1:3),gF16)
call fun_my_cross(grad_snF(1,1:3),grad_snF(4,1:3),gF14)
call fun_my_cross(grad_snF(4,1:3),grad_snF(3,1:3),gF43)
call fun_my_cross(grad_snF(3,1:3),grad_snF(2,1:3),gF32)
call fun_my_cross(grad_snF(2,1:3),grad_snF(1,1:3),gF21)
call fun_my_cross(grad_snF(1,1:3),grad_snF(2,1:3),gF12)
call fun_my_cross(grad_snF(2,1:3),grad_snF(3,1:3),gF23)
call fun_my_cross(grad_snF(3,1:3),grad_snF(4,1:3),gF34)
call fun_my_cross(grad_snF(4,1:3),grad_snF(1,1:3),gF41)
!!
do ii = 1,NP
J(1:3,1)=(-(1.0d0-eta(ii))*(1.0d0+zet(ii))*Hexa(1,1:3)-(1.0d0-eta(ii))*(1.0d0-zet(ii))*Hexa(2,1:3) &
		  +(1.0d0-eta(ii))*(1.0d0-zet(ii))*Hexa(3,1:3)+(1.0d0-eta(ii))*(1.0d0+zet(ii))*Hexa(4,1:3) &
		  -(1.0d0+eta(ii))*(1.0d0+zet(ii))*Hexa(5,1:3)-(1.0d0+eta(ii))*(1.0d0-zet(ii))*Hexa(6,1:3) &
		  +(1.0d0+eta(ii))*(1.0d0-zet(ii))*Hexa(7,1:3)+(1.0d0+eta(ii))*(1.0d0+zet(ii))*Hexa(8,1:3))*0.125d0

J(1:3,2)=(-(1.0d0- xi(ii))*(1.0d0+zet(ii))*Hexa(1,1:3)-(1.0d0- xi(ii))*(1.0d0-zet(ii))*Hexa(2,1:3) &
		  -(1.0d0+ xi(ii))*(1.0d0-zet(ii))*Hexa(3,1:3)-(1.0d0+ xi(ii))*(1.0d0+zet(ii))*Hexa(4,1:3) &
		  +(1.0d0- xi(ii))*(1.0d0+zet(ii))*Hexa(5,1:3)+(1.0d0- xi(ii))*(1.0d0-zet(ii))*Hexa(6,1:3) &
		  +(1.0d0+ xi(ii))*(1.0d0-zet(ii))*Hexa(7,1:3)+(1.0d0+ xi(ii))*(1.0d0+zet(ii))*Hexa(8,1:3))*0.125d0

J(1:3,3)=(+(1.0d0- xi(ii))*(1.0d0-eta(ii))*Hexa(1,1:3)-(1.0d0- xi(ii))*(1.0d0-eta(ii))*Hexa(2,1:3) &
		  -(1.0d0+ xi(ii))*(1.0d0-eta(ii))*Hexa(3,1:3)+(1.0d0+ xi(ii))*(1.0d0-eta(ii))*Hexa(4,1:3) &
		  +(1.0d0- xi(ii))*(1.0d0+eta(ii))*Hexa(5,1:3)-(1.0d0- xi(ii))*(1.0d0+eta(ii))*Hexa(6,1:3) &
		  -(1.0d0+ xi(ii))*(1.0d0+eta(ii))*Hexa(7,1:3)+(1.0d0+ xi(ii))*(1.0d0+eta(ii))*Hexa(8,1:3))*0.125d0
sn(1)=1.0d0/8.0d0*(1.0d0-xi(ii))*(1.0d0-eta(ii))*(1.0d0+zet(ii))
sn(2)=1.0d0/8.0d0*(1.0d0-xi(ii))*(1.0d0-eta(ii))*(1.0d0-zet(ii))
sn(3)=1.0d0/8.0d0*(1.0d0+xi(ii))*(1.0d0-eta(ii))*(1.0d0-zet(ii))
sn(4)=1.0d0/8.0d0*(1.0d0+xi(ii))*(1.0d0-eta(ii))*(1.0d0+zet(ii))
sn(5)=1.0d0/8.0d0*(1.0d0-xi(ii))*(1.0d0+eta(ii))*(1.0d0+zet(ii))
sn(6)=1.0d0/8.0d0*(1.0d0-xi(ii))*(1.0d0+eta(ii))*(1.0d0-zet(ii))
sn(7)=1.0d0/8.0d0*(1.0d0+xi(ii))*(1.0d0+eta(ii))*(1.0d0-zet(ii))
sn(8)=1.0d0/8.0d0*(1.0d0+xi(ii))*(1.0d0+eta(ii))*(1.0d0+zet(ii))
wf_loc(1:3,3)=sn(4)*gF61+sn(3)*gF15+sn(7)*gF53+sn(8)*gF36
wf_loc(1:3,1)=sn(5)*gF63+sn(6)*gF35+sn(2)*gF51+sn(1)*gF16
wf_loc(1:3,2)=sn(2)*gF14+sn(6)*gF43+sn(7)*gF32+sn(3)*gF21
wf_loc(1:3,4)=sn(4)*gF12+sn(8)*gF23+sn(5)*gF34+sn(1)*gF41
call det3(J,dett)
do jj = 1,4
	call mat33_x_vec3(J(1:3,1:3),wf_loc(1:3,jj),tmp3(1:3))
	wf(1:3,ii,jj)=(1.d0/dett)*tmp3(1:3)
enddo
enddo
!!
end subroutine fun_whitney_face_hexa_6_quad
!!
subroutine fun_my_cross(a, b, c)
real*8, DIMENSION(3) :: c
real*8, DIMENSION(3), INTENT(IN) :: a, b
c(1) = a(2) * b(3) - a(3) * b(2)
c(2) = a(3) * b(1) - a(1) * b(3)
c(3) = a(1) * b(2) - a(2) * b(1)
end subroutine fun_my_cross
!!
function CP(a,b) result(c)
implicit none
real(kind=8), dimension(3), intent(in) :: a,b
real(kind=8), dimension(3)    :: c
c(1)=a(2)*b(3)-a(3)*b(2)
c(2)=a(3)*b(1)-a(1)*b(3)
c(3)=a(1)*b(2)-a(2)*b(1)
end function
!!
subroutine det3(matrix,det)
	implicit none
	real(kind=8), dimension(3,3) :: matrix
	real(kind=8)				 :: det
	det=( matrix(3,3)*matrix(1,1)*matrix(2,2)-matrix(1,1)*matrix(2,3)*matrix(3,2) &
	 	 & -matrix(3,3)*matrix(2,1)*matrix(1,2)+matrix(2,1)*matrix(1,3)*matrix(3,2) &
		 & +matrix(3,1)*matrix(1,2)*matrix(2,3)-matrix(3,1)*matrix(1,3)*matrix(2,2) )
end subroutine	

subroutine funTrilinear(xi,eta,zet,QQ,PP)
implicit none
integer*8 ii
real*8 xi, eta, zet, QQ(8,3), PP(3), sn(8)
!! orientazione Paolo
! kcsi=[-1,+1,+1,-1,-1,+1,+1,-1]
! keta=[-1,-1,+1,+1,-1,-1,+1,+1]
! kzeta=[-1,-1,-1,-1,+1,+1,+1,+1]
! sn=0.125*(1+kcsi*csi).*(1+keta*eta).*(1+kzeta*zeta)
! PP(1)=sn*QQ(:,1)
! PP(2)=sn*QQ(:,2)
! PP(3)=sn*QQ(:,3)
!! orientazione Riccardo
sn(1)=1.0d0/8.0d0*(1.0d0-xi)*(1.0d0-eta)*(1.0d0+zet)
sn(2)=1.0d0/8.0d0*(1.0d0-xi)*(1.0d0-eta)*(1.0d0-zet)
sn(3)=1.0d0/8.0d0*(1.0d0+xi)*(1.0d0-eta)*(1.0d0-zet)
sn(4)=1.0d0/8.0d0*(1.0d0+xi)*(1.0d0-eta)*(1.0d0+zet)
sn(5)=1.0d0/8.0d0*(1.0d0-xi)*(1.0d0+eta)*(1.0d0+zet)
sn(6)=1.0d0/8.0d0*(1.0d0-xi)*(1.0d0+eta)*(1.0d0-zet)
sn(7)=1.0d0/8.0d0*(1.0d0+xi)*(1.0d0+eta)*(1.0d0-zet)
sn(8)=1.0d0/8.0d0*(1.0d0+xi)*(1.0d0+eta)*(1.0d0+zet)
PP(1:3)=0.0
do ii = 1,8
PP(1)=PP(1)+sn(ii)*QQ(ii,1)
PP(2)=PP(2)+sn(ii)*QQ(ii,2)
PP(3)=PP(3)+sn(ii)*QQ(ii,3)
enddo 
end subroutine funTrilinear