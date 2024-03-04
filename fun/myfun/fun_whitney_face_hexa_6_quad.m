function [wf] = fun_whitney_face_hexa_6_quad(Hexa,P,NP)
%%
xi =P(1,1:NP);
eta=P(2,1:NP);
zet=P(3,1:NP);
%%
grad_snF=zeros(6,3);
grad_snF(1,1:3)=[0,-0.5,0];
grad_snF(2,1:3)=[0.5,0,0];
grad_snF(3,1:3)=[0,0.5d0,0];
grad_snF(4,1:3)=[-0.5,0,0];
grad_snF(5,1:3)=[0,0,-0.5];
grad_snF(6,1:3)=[0,0,0.5];
%%
wf=zeros(3,NP,4);
for ii = 1:NP
J(1:3,1)=(-(1-eta(ii))*(1+zet(ii))*Hexa(1,1:3)-(1-eta(ii))*(1-zet(ii))*Hexa(2,1:3) ...
		  +(1-eta(ii))*(1-zet(ii))*Hexa(3,1:3)+(1-eta(ii))*(1+zet(ii))*Hexa(4,1:3) ...
		  -(1+eta(ii))*(1+zet(ii))*Hexa(5,1:3)-(1+eta(ii))*(1-zet(ii))*Hexa(6,1:3) ...
		  +(1+eta(ii))*(1-zet(ii))*Hexa(7,1:3)+(1+eta(ii))*(1+zet(ii))*Hexa(8,1:3))*0.125;

J(1:3,2)=(-(1- xi(ii))*(1+zet(ii))*Hexa(1,1:3)-(1- xi(ii))*(1-zet(ii))*Hexa(2,1:3) ...
		  -(1+ xi(ii))*(1-zet(ii))*Hexa(3,1:3)-(1+ xi(ii))*(1+zet(ii))*Hexa(4,1:3) ...
		  +(1- xi(ii))*(1+zet(ii))*Hexa(5,1:3)+(1- xi(ii))*(1-zet(ii))*Hexa(6,1:3) ...
		  +(1+ xi(ii))*(1-zet(ii))*Hexa(7,1:3)+(1+ xi(ii))*(1+zet(ii))*Hexa(8,1:3))*0.125;

J(1:3,3)=(+(1- xi(ii))*(1-eta(ii))*Hexa(1,1:3)-(1- xi(ii))*(1-eta(ii))*Hexa(2,1:3) ...
		  -(1+ xi(ii))*(1-eta(ii))*Hexa(3,1:3)+(1+ xi(ii))*(1-eta(ii))*Hexa(4,1:3) ...
		  +(1- xi(ii))*(1+eta(ii))*Hexa(5,1:3)-(1- xi(ii))*(1+eta(ii))*Hexa(6,1:3) ...
		  -(1+ xi(ii))*(1+eta(ii))*Hexa(7,1:3)+(1+ xi(ii))*(1+eta(ii))*Hexa(8,1:3))*0.125;
%%
sn(1)=0.125*(1-xi(ii))*(1-eta(ii))*(1+zet(ii));
sn(2)=0.125*(1-xi(ii))*(1-eta(ii))*(1-zet(ii));
sn(3)=0.125*(1+xi(ii))*(1-eta(ii))*(1-zet(ii));
sn(4)=0.125*(1+xi(ii))*(1-eta(ii))*(1+zet(ii));
sn(5)=0.125*(1-xi(ii))*(1+eta(ii))*(1+zet(ii));
sn(6)=0.125*(1-xi(ii))*(1+eta(ii))*(1-zet(ii));
sn(7)=0.125*(1+xi(ii))*(1+eta(ii))*(1-zet(ii));
sn(8)=0.125*(1+xi(ii))*(1+eta(ii))*(1+zet(ii)); 
%%
wf_loc(1:3,3)=sn(4)*cross(grad_snF(6,1:3),grad_snF(1,1:3))+sn(3)*cross(grad_snF(1,1:3),...
    grad_snF(5,1:3))+sn(7)*cross(grad_snF(5,1:3),grad_snF(3,1:3))+sn(8)*cross(grad_snF(3,1:3),grad_snF(6,1:3));
wf_loc(1:3,1)=sn(5)*cross(grad_snF(6,1:3),grad_snF(3,1:3))+sn(6)*cross(grad_snF(3,1:3),...
    grad_snF(5,1:3))+sn(2)*cross(grad_snF(5,1:3),grad_snF(1,1:3))+sn(1)*cross(grad_snF(1,1:3),grad_snF(6,1:3));
wf_loc(1:3,2)=sn(2)*cross(grad_snF(1,1:3),grad_snF(4,1:3))+sn(6)*cross(grad_snF(4,1:3),...
    grad_snF(3,1:3))+sn(7)*cross(grad_snF(3,1:3),grad_snF(2,1:3))+sn(3)*cross(grad_snF(2,1:3),grad_snF(1,1:3));
wf_loc(1:3,4)=sn(4)*cross(grad_snF(1,1:3),grad_snF(2,1:3))+sn(8)*cross(grad_snF(2,1:3),...
    grad_snF(3,1:3))+sn(5)*cross(grad_snF(3,1:3),grad_snF(4,1:3))+sn(1)*cross(grad_snF(4,1:3),grad_snF(1,1:3));
%% 
for jj = 1:4
wf(1:3,ii,jj)=(1/det(J))*(J(1:3,1:3))*wf_loc(1:3,jj);
end 
end
end

