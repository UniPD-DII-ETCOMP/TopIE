function [R] = funRrz_face2(N_face,rho,N_GAUSS,Matrix_P0,F1,N_edge,C1,Fac_Ed_loc,Aed,Led)
%%
[xabsc,weig]=lgwt(N_GAUSS,-1,1);
NP=N_GAUSS^2;
xi=zeros(NP,1);
eta=-ones(NP,1);
zet=zeros(NP,1);
WG_loc=zeros(NP,1);
hh=1;
for ii = 1:N_GAUSS  
    for jj = 1:N_GAUSS
           xi(hh) =xabsc(ii);
           zet(hh)=xabsc(jj);
           WG_loc(hh)= weig(ii)*weig(jj);
           hh=hh+1;
    end 
end
%%
Hexa_master=[-1,-1, 1;... 
             -1,-1,-1;...
              1,-1,-1;...
              1,-1, 1;...
             -1,1, 1;...
             -1,1,-1;...
              1,1,-1;...
              1,1, 1];
[wf_loc] = 2*fun_whitney_face_hexa_6_quad(Hexa_master,[xi,eta,zet].',NP);
%%
R=zeros(N_edge,N_edge);
Hexa=zeros(8,3);
PP=zeros(NP,3);
detJ=zeros(NP,1);
we_glo=zeros(3,NP,4);
for ii = 1:N_face 
    Face=Matrix_P0(F1(1:4,ii),:);
    Hexa(1:4,[1,3])=Face(1:4,[1,3]);
    Hexa(5:8,[1,3])=Face(1:4,[1,3]);
    Hexa(1:4,2)=-1;
    Hexa(5:8,2)=1;
    for jj = 1:NP 
        [PP(jj,:)]=funTrilinear(xi(jj),0.0,zet(jj),Hexa);
    end
    for jj = 1
        J(1:3,1)=(-(1.0d0-eta(jj))*(1+zet(jj))*Hexa(1,1:3)-(1-eta(jj))*(1-zet(jj))*Hexa(2,1:3) ...
                  +(1.0d0-eta(jj))*(1-zet(jj))*Hexa(3,1:3)+(1-eta(jj))*(1+zet(jj))*Hexa(4,1:3) ...
                  -(1.0d0+eta(jj))*(1+zet(jj))*Hexa(5,1:3)-(1+eta(jj))*(1-zet(jj))*Hexa(6,1:3) ...
                  +(1.0d0+eta(jj))*(1-zet(jj))*Hexa(7,1:3)+(1+eta(jj))*(1+zet(jj))*Hexa(8,1:3))*0.125;

        J(1:3,2)=(-(1.0d0- xi(jj))*(1+zet(jj))*Hexa(1,1:3)-(1- xi(jj))*(1-zet(jj))*Hexa(2,1:3) ...
                  -(1.0d0+ xi(jj))*(1-zet(jj))*Hexa(3,1:3)-(1+ xi(jj))*(1+zet(jj))*Hexa(4,1:3) ...
                  +(1.0d0- xi(jj))*(1+zet(jj))*Hexa(5,1:3)+(1- xi(jj))*(1-zet(jj))*Hexa(6,1:3) ...
                  +(1.0d0+ xi(jj))*(1-zet(jj))*Hexa(7,1:3)+(1+ xi(jj))*(1+zet(jj))*Hexa(8,1:3))*0.125;

        J(1:3,3)=(+(1.0d0- xi(jj))*(1-eta(jj))*Hexa(1,1:3)-(1- xi(jj))*(1-eta(jj))*Hexa(2,1:3) ...
                  -(1.0d0+ xi(jj))*(1-eta(jj))*Hexa(3,1:3)+(1+ xi(jj))*(1-eta(jj))*Hexa(4,1:3) ...
                  +(1.0d0- xi(jj))*(1+eta(jj))*Hexa(5,1:3)-(1- xi(jj))*(1+eta(jj))*Hexa(6,1:3) ...
                  -(1.0d0+ xi(jj))*(1+eta(jj))*Hexa(7,1:3)+(1+ xi(jj))*(1+eta(jj))*Hexa(8,1:3))*0.125;
       detJ(:)=det(J); 
    end
    for jj=1:4
        for kk = 1:NP
            we_glo(1:3,kk,jj)=(1/detJ(kk))*(J(1:3,1:3))*wf_loc(1:3,kk,jj);
        end
    end
    for kk = 1:4 
        ind_e1s=C1(kk,ii);
        ind_e1=abs(ind_e1s); 
        ind_loc_wk=Fac_Ed_loc(kk,ii);
        for hh = 1:4
            ind_e2s=(C1(hh,ii));
            ind_e2=abs(ind_e2s);
            ind_loc_wh=Fac_Ed_loc(hh,ii);
            Rcoeff=0.0;
            for ee = 1:NP
                Rcoeff=Rcoeff+2*pi*PP(ee,1)...
                              *rho(ii)...
                              *dot(we_glo(1:3,ee,ind_loc_wk)*Led(ind_e1)/Aed(ind_e1),we_glo(1:3,ee,ind_loc_wh)*Led(ind_e2)/Aed(ind_e2))...
                              *WG_loc(ee)*detJ(ee)...
                              *sign(ind_e1s)*sign(ind_e2s);
            end
            R(ind_e1,ind_e2)=R(ind_e1,ind_e2)+Rcoeff;
        end
    end       
end
end