function [h_ext,j_ext] = fun_compute_ext_field(N,Matrix_Projm,Matrix_G_m,Matrix_Nme,bar_e,Area_e,J_ext)
j_ext=zeros(N.face_ext,1);
for ii = 1:N.face_ext
    jextii=J_ext(bar_e(1:3,ii));
    j_ext(ii)=jextii(2)*Area_e(ii);
end
h_ext=-Matrix_Projm*Matrix_G_m*Matrix_Nme*j_ext;