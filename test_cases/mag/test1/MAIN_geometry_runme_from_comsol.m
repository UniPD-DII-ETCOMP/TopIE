%%
clear 
close all
clc
%% BEGIN USER SETTINGS
C = textread('data.mphtxt', '%s','delimiter', '\n'); %#ok<DTXTRD>
%% END USER SETTINGS
dad=pwd; cd ..; cd ..; cd ..; cd('fun'); cd('myfun'); cd('GCD'); addpath(pwd); cd(dad)
[G1,C1,D1,F1,VP,Matrix_P0,Matrix_C,Matrix_G,Matrix_D] = ...
    fun_extract_GCD(C);
%% barycenter
N.face=size(F1,2);
barf=zeros(3,N.face);
for ii = 1:N.face
    barf(1:3,ii)=(sum(Matrix_P0(F1(1:4,ii),:))/4).';
end
%%
tol=1e-3;
indf=find(abs(barf(2,:))<tol & barf(1,:)>0);
%%
F1=F1(1:4,indf);
C1=C1(:,indf);
ed=unique(abs(reshape(C1,size(C1,1)*size(C1,2),1)));
%%
ed_ori=1:size(G1,2);
ed_map=zeros(size(G1,2),1);
ed_map(ed)=1:length(ed);
for ii = 1:size(C1,2)
   C1(:,ii)=ed_map(abs(C1(:,ii))).*sign(C1(:,ii));
end
G1=G1(:,ed);
%%
nod=unique(abs(reshape(G1,size(G1,1)*size(G1,2),1)));
nd_map=zeros(size(Matrix_P0,1),1);
nd_map(nod)=1:length(nod);
for ii = 1:size(G1,2)
   G1(:,ii)=nd_map(abs(G1(:,ii))).*sign(G1(:,ii));
end
for ii = 1:size(F1,2)
   F1(:,ii)=nd_map(abs(F1(:,ii))).*sign(F1(:,ii)); 
end
Matrix_P0=Matrix_P0(nod,:);
%%
Matrix_G=Matrix_G(ed,nod);
Matrix_C=Matrix_C(indf,ed);
%%
N.edge=size(G1,2);
N.face=size(C1,2);
%%
figure
patch('Faces',F1(1:4,:).','Vertices',Matrix_P0,'Facecolor','b','FaceAlpha',0.2)
axis equal
view(0,0)
axis([min(Matrix_P0(:,1)) max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
%%
G2=fun_G2(F1,C1,N,G1);
%%
figure
patch('Faces',F1(1:4,:).','Vertices',Matrix_P0,'Facecolor','g','FaceAlpha',0.2)
axis equal
view(0,0)
axis([min(Matrix_P0(:,1)) max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
%% check normal vector
vec=zeros(N.face,3);
for ii = 1:N.face
   vec(ii,1:3)=cross(Matrix_P0(F1(2,ii),:)-Matrix_P0(F1(1,ii),:),...
                     Matrix_P0(F1(3,ii),:)-Matrix_P0(F1(2,ii),:));
   if vec(ii,2)>0
      F1(:,ii)=flip(F1(:,ii)); 
      vec(ii,2)=-vec(ii,2);
   end
end
%%
C1_m=C1;
G1_m=G1;
Matrix_G_m=Matrix_G;
Matrix_C_m=Matrix_C;
Matrix_P0_m=Matrix_P0;
F1_m=F1;
G2_m=G2;
%%
figure
hold on
patch('Faces',F1(1:4,:).','Vertices',Matrix_P0,'Facecolor','b','FaceAlpha',0.2)
quiver3(barf(1,indf),barf(2,indf),barf(3,indf),...
        vec(:,1).',vec(:,2).',vec(:,3).')
axis equal
view(0,0)
axis([min(Matrix_P0(:,1))*0 max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
%%
save data_mag.mat Matrix_G_m G2_m C1_m F1_m G1_m Matrix_C_m  Matrix_P0_m  




