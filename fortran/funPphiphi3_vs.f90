subroutine funPphiphi3_vs(N_node,N_edge,N_GAUSS1,N_GAUSS2,N_thread,N_face,Matrix_P0,G1,F1,Area,Area_vol,P)
implicit none
integer*8 N_node, N_edge, N_face, N_GAUSS1, N_GAUSS2, N_thread
real*8 Matrix_P0(N_node,3), Area(N_face), Area_vol(N_face)
integer*8 G1(2,N_edge)
integer*8 hh, ii, kk, qq, ee, jj
integer*8 NP1, NP2, F1(4,N_face)
real*8 P(N_edge,N_face)
real*8 Kg2, k2, Kg, k, e1, e2, G2Daxi, pi
real*8 xabsc1(N_GAUSS1), weig1(N_GAUSS1), xi1(N_GAUSS1), WG_loc1(N_GAUSS1), PP1(N_GAUSS1,3)
real*8 xabsc2(N_GAUSS2), weig2(N_GAUSS2), xi2(N_GAUSS2**2), WG_loc2(N_GAUSS2**2), PP2(N_GAUSS2**2,3)
real*8 zet2(N_GAUSS2**2), eta2(N_GAUSS2**2)
real*8 Ed1(2,3), Cent1(3), vec1(3), lung1
real*8 Face2(4,3), dett, Hexa2(8,3), xyz(3)
real*8 J2(3,3), detJ2(N_GAUSS2**2), Pint
real*8 ree, rqq, h
!!
pi=4.d0*datan(1.d0)
!!
call gauleg(int(N_GAUSS1,8),xabsc1(1:N_GAUSS1),weig1(1:N_GAUSS1))
NP1=N_GAUSS1; 
hh=1;
do ii = 1,N_GAUSS1 ! 
   xi1(hh) =xabsc1(ii)
   WG_loc1(hh)= weig1(ii)
   hh=hh+1
enddo
!!
call gauleg(int(N_GAUSS2,8),xabsc2(1:N_GAUSS2),weig2(1:N_GAUSS2))
NP2=N_GAUSS2**2; 
hh=1;
eta2(1:NP2)=0.0d0
do ii = 1,N_GAUSS2  
    do jj = 1,N_GAUSS2
           xi2(hh) =xabsc2(ii)
           zet2(hh)=xabsc2(jj)
           WG_loc2(hh)= weig2(ii)*weig2(jj)
           hh=hh+1
    enddo 
enddo
!!
P(1:N_edge,1:N_face)=0.0
call omp_set_num_threads(N_thread)
!$OMP PARALLEL SHARED(N_face,N_edge,Matrix_P0,G1,F1,NP1,xi1,NP2,xi2,eta2,zet2,pi,Area,Area_vol,P)
!$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(ii,Ed1,Cent1,vec1,lung1,kk,PP1,jj, & 
!$OMP                                PP2,Pint,qq,ee,ree,rqq,h,Kg2,k2,Kg,k,e1,e2,G2Daxi,Face2, &
!$OMP                                     Hexa2,xyz,J2,dett,detJ2)
do ii = 1,N_edge ! ciclo sui lati target
!! lato target
    Ed1=Matrix_P0(G1(1:2,ii),1:3);
    Cent1=0.5*(Ed1(1,1:3)+Ed1(2,1:3));
    vec1=(Ed1(2,1:3)-Ed1(1,1:3))*0.5;
    lung1=norm2(vec1)*2.0;
!! punti di gauss in globale
    do kk = 1,NP1 
        PP1(kk,:)=Cent1+vec1*xi1(kk);
    enddo
!! 
    do jj = 1,N_face
        Face2=Matrix_P0(F1(1:4,jj),1:3);
!         Hexa2=zeros(8,3);
        Hexa2(1:4,[1,3])=Face2(1:4,[1,3]);
        Hexa2(5:8,[1,3])=Face2(1:4,[1,3]);
        Hexa2(1:4,2)=-1.0;
        Hexa2(5:8,2)=1.0;
        do kk = 1,NP2 ! punti di gauss in globale
            call funTrilinear(xi2(kk),0.0d0,zet2(kk),Hexa2,xyz);
			PP2(kk,1:3)=xyz
        enddo

		kk = 1
        !do kk = 1,NP2
            J2(1:3,1)=(-(1.0d0-eta2(kk))*(1.0d0+zet2(kk))*Hexa2(1,1:3)-(1.0d0-eta2(kk))*(1.0d0-zet2(kk))*Hexa2(2,1:3) &
                       +(1.0d0-eta2(kk))*(1.0d0-zet2(kk))*Hexa2(3,1:3)+(1.0d0-eta2(kk))*(1.0d0+zet2(kk))*Hexa2(4,1:3) &
                       -(1.0d0+eta2(kk))*(1.0d0+zet2(kk))*Hexa2(5,1:3)-(1.0d0+eta2(kk))*(1.0d0-zet2(kk))*Hexa2(6,1:3) &
                       +(1.0d0+eta2(kk))*(1.0d0-zet2(kk))*Hexa2(7,1:3)+(1.0d0+eta2(kk))*(1.0d0+zet2(kk))*Hexa2(8,1:3))*0.125d0;

            J2(1:3,2)=(-(1.0d0- xi2(kk))*(1.0d0+zet2(kk))*Hexa2(1,1:3)-(1.0d0- xi2(kk))*(1.0d0-zet2(kk))*Hexa2(2,1:3) &
                       -(1.0d0+ xi2(kk))*(1.0d0-zet2(kk))*Hexa2(3,1:3)-(1.0d0+ xi2(kk))*(1.0d0+zet2(kk))*Hexa2(4,1:3) &
                       +(1.0d0- xi2(kk))*(1.0d0+zet2(kk))*Hexa2(5,1:3)+(1.0d0- xi2(kk))*(1.0d0-zet2(kk))*Hexa2(6,1:3) &
                       +(1.0d0+ xi2(kk))*(1.0d0-zet2(kk))*Hexa2(7,1:3)+(1.0d0+ xi2(kk))*(1.0d0+zet2(kk))*Hexa2(8,1:3))*0.125d0;

            J2(1:3,3)=(+(1.0d0- xi2(kk))*(1.0d0-eta2(kk))*Hexa2(1,1:3)-(1.0d0- xi2(kk))*(1.0d0-eta2(kk))*Hexa2(2,1:3) &
                       -(1.0d0+ xi2(kk))*(1.0d0-eta2(kk))*Hexa2(3,1:3)+(1.0d0+ xi2(kk))*(1.0d0-eta2(kk))*Hexa2(4,1:3) &
                       +(1.0d0- xi2(kk))*(1.0d0+eta2(kk))*Hexa2(5,1:3)-(1.0d0- xi2(kk))*(1.0d0+eta2(kk))*Hexa2(6,1:3) &
                       -(1.0d0+ xi2(kk))*(1.0d0+eta2(kk))*Hexa2(7,1:3)+(1.0d0+ xi2(kk))*(1.0d0+eta2(kk))*Hexa2(8,1:3))*0.125d0;
           call det3(J2,dett);
		   detJ2(:)=dett		   
        !enddo 
        ! 
        Pint=0.0d0;
        do qq = 1,NP1
            do ee = 1,NP2
                ree = PP2(ee,1); !source point r
                rqq = PP1(qq,1); !target point r
                h = PP2(ee,3)-PP1(qq,3);
                Kg2=(ree+rqq)**2.0+h**2.0;
                k2=4.0*ree*rqq/Kg2;
                Kg=sqrt(Kg2);
                k=sqrt(k2);
                call my_ellipke(k2,e1,e2);
                G2Daxi=(ree)*(4.d0*e1/Kg)/(4.d0*pi);				
				!
                Pint=Pint+  + G2Daxi*WG_loc1(qq)*WG_loc2(ee)*detJ2(qq)*lung1*0.5;
            enddo
        enddo
        Pint=Pint/(lung1*Area_vol(jj));
        P(ii,jj)=Pint;
    enddo
enddo
!$OMP END DO NOWAIT 
!$OMP END PARALLEL 
end subroutine funPphiphi3_vs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
tol=1d-20!2.2204d-16
pi=4.d0*datan(1.d0)
a0 = 1.0;
b0 = sqrt(1.0-m);
c0 = 0.0;
s0 = m;
i1 = 0.0; 
mm = 666.0d6;
do while (mm >= tol) 
    a1 = (a0+b0)/2.0;
    b1 = sqrt(a0*b0);
    c1 = (a0-b0)/2.0;
    i1 = i1 + 1.0;
    w1 = (2.0**i1)*(c1**2.0);
    mm = w1
    
    ! test do stagnation (may happen do TOL < machine precision)
    if (c0==c1) then
	    write(66,*) 'my_ellipke fails'
        stop 
    end if
    
    s0 = s0 + w1;  
    a0 = a1;  
	b0 = b1;  
	c0 = c1;
end do
k = pi/(2.0*a1);
e = k*(1.0-s0/2.0);
end subroutine my_ellipke
!!
subroutine L_self_Axial_rect_cross(rm,dr,dz,L)
implicit none
real*8 rm,dr,dz,L, pi
pi=4.d0*datan(1.d0)
!!
L = rm*(&
       log(8*rm/dr) + 1.0/12.0 - pi*dz/(3.0*dr) - 0.5*log(1+(dz**2)/(dr**2))&
       +(1.0/12.0)*(dr**2.0/dz**2.0)*log(1.0+dz**2/dr**2)+(1.0/12.0)*dz**2.0/dr**2.0*log(1.0+dr**2.0/dz**2.0)&
       +(2.0/3.0)*(dz/dr-dr/dz)*datan2(dz,dr)+dr**2.0/(96.0*rm**2.0)*((log(8.0*rm/dr)&
       -0.5*log(1.0+dz**2.0/dr**2.0))*(1.0+3.0*dz**2.0/dr**2.0)+3.45*dz**2.0/dr**2.0 &
       +221.0/60.0 -1.6*pi*dz**3.0/dr**3.0 +3.2*dz**3.0/dr**3.0*datan2(dz,dr)-0.1&
       *dr**2.0/dz**2.0*log(1.0+dz**2.0/dr**2.0)+0.5*dz**4.0/dr**4.0*log(1.0+dr**2.0/dz**2.0)));
end subroutine L_self_Axial_rect_cross
!--------------------------------------------------------------------------------
! Calculation of GAUSS-LEGENDRE abscissas and weights do Gaussian Quadrature
! integration of polynomial functions.
!      For normalized lower and upper limits of integration -1.0 & 1.0, and
! given n, this routine calculates, arrays xabsc(1:n) and  weig(1:n) of length n,
! containing the abscissas and weights of the Gauss-Legendre n-point quadrature
! formula.
!--------------------------------------------------------------------------------
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
subroutine funTrilinear(xi,eta,zet,QQ,PP)
implicit none
integer*8 ii
real*8 xi, eta, zet, QQ(8,3), PP(3), sn(8)
!! orientazione Paolo
! kcsi=[-1,+1,+1,-1,-1,+1,+1,-1];
! keta=[-1,-1,+1,+1,-1,-1,+1,+1];
! kzeta=[-1,-1,-1,-1,+1,+1,+1,+1];
! sn=0.125*(1+kcsi*csi).*(1+keta*eta).*(1+kzeta*zeta);
! PP(1)=sn*QQ(:,1);
! PP(2)=sn*QQ(:,2);
! PP(3)=sn*QQ(:,3);
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
!!
subroutine det3(matrix,det)
	implicit none
	real(kind=8), dimension(3,3) :: matrix
	real(kind=8)				 :: det
	det=( matrix(3,3)*matrix(1,1)*matrix(2,2)-matrix(1,1)*matrix(2,3)*matrix(3,2) &
	 	 & -matrix(3,3)*matrix(2,1)*matrix(1,2)+matrix(2,1)*matrix(1,3)*matrix(3,2) &
		 & +matrix(3,1)*matrix(1,2)*matrix(2,3)-matrix(3,1)*matrix(1,3)*matrix(2,2) )
end subroutine	