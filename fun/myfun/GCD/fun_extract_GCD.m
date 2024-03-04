function [G1,C1,D1,F1,VP,Matrix_P0,Matrix_C,Matrix_G,Matrix_D] = ...
    fun_extract_GCD(C)
for ii = 1:size(C,1)
    CC = char(C(ii));
   if strcmp(CC,'# Mesh vertex coordinates')
       ind_start_point = ii+1;
   end
   if strcmp(CC,'4 # number of vertices per element')
       ind_start_q = ii+3;
   end
end
Npoint= C(18)
Ntetra = C(ind_start_q-2)
CC=char(C(18));
ind=find(CC=='#');
Np=str2double(CC(1:ind-1))
CC=char(C(ind_start_q-2));
ind=find(CC=='#');
Nv=str2double(CC(1:ind-1))
a = str2num(char(C(ind_start_point:ind_start_point+Np)));
Matrix_P0(:,1)=a(:,1);
Matrix_P0(:,3)=a(:,2);
VP = str2num(char(C(ind_start_q:ind_start_q+Nv))); 
VP = VP+1;
check = max(max(VP)) == Np
%%
nn=max(max(VP));
VP=[VP,VP+nn];
Matrix_P0=[Matrix_P0;Matrix_P0];
Matrix_P0(nn+1:end,2)=Matrix_P0(nn+1:end,2)+1;
 VP = VP.';

%%

reorder = 1;
permute_Mat=[5 1 2 6 7 3 4 8];

%% riordinamento punti
if reorder == 1 
    for ii=1:size(VP,2)
        VP(:,ii)=VP(permute_Mat,ii);
    end
end
%%

[G1,C1,D1,F1]=gcd_mexed_HEXA(VP);

%% number of nodes edges faces volumes

nn = max(max(VP))
ne = size(G1,2)
nf = size(C1,2)
nv = size(D1,2)

%% G edges xnodes

 Matrix_G = sparse(1:ne,G1(1,1:ne),-ones(ne,1),ne,nn);
 Matrix_G = Matrix_G+sparse(1:ne,G1(2,1:ne),ones(ne,1),ne,nn);
%  Gnew = full(Gnew);

 %% C faces x edges

 Cpos = (C1+abs(C1))/2;
 Cneg = (C1-abs(C1))/2;
 [r_pos_c,c_pos_c,val_pos_c] = find(Cpos);
 [r_neg_c,c_neg_c,val_neg_c] = find(Cneg);
 Matrix_C = sparse(c_pos_c,val_pos_c,ones(size(c_pos_c,1),1),nf,ne);
 Matrix_C = Matrix_C+sparse(c_neg_c,abs(val_neg_c),-ones(size(c_neg_c,1),1),nf,ne);
%  Cnew = full(Cnew);
 
 %% check rot(grad())
 disp('check rot(grad()) = 0...')
 CxG=Matrix_C*Matrix_G;
 find(CxG)
 
 %% D volumes x faces

 Dpos = (D1+abs(D1))/2;
 Dneg = (D1-abs(D1))/2;
 [r_pos_d,c_pos_d,val_pos_d] = find(Dpos);
 [r_neg_d,c_neg_d,val_neg_d] = find(Dneg);
 Matrix_D = sparse(c_pos_d,val_pos_d,ones(size(c_pos_d,1),1),nv,nf);
 Matrix_D = Matrix_D+sparse(c_neg_d,abs(val_neg_d),-ones(size(c_neg_d,1),1),nv,nf);
%  Dnew = full(Dnew);
 
 %% check rot(grad())
 disp('check div(rot()) = 0...')
 DxC=Matrix_D*Matrix_C;
 find(DxC)
 %% check
 
 %% save
end