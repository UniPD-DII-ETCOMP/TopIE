function [out]=fun_solve_system(x,...
                                freq,...
                                fun_Matrix_Rm,...
                                Matrix_Ca_m,...
                                Matrix_Pm_vv,Matrix_Pm_vs,Matrix_Pm_ss,...
                                Matrix_Projm,Matrix_G_m,Matrix_Nme,...
                                h_ext,j_ext,...
                                opt)
%%
mu_0=4*pi*1e-7;
w=2*pi*freq;
%%
F_mur=opt.murmin+(opt.murmax-opt.murmin)*x;
param=1./(mu_0*(F_mur-1));
Matrix_Rm=fun_Matrix_Rm(param);
S=(1/(1j*w))*Matrix_Rm+...
    (1/(1j*w))*(Matrix_Ca_m.'*[Matrix_Pm_vv,Matrix_Pm_vs.';Matrix_Pm_vs,Matrix_Pm_ss]*Matrix_Ca_m);
Xm=S\h_ext;
%%
objf=fun_eval_objective(Xm,Matrix_Projm,Matrix_G_m,Matrix_Nme,j_ext,freq);
%%
out.u=Xm;
out.rhs=h_ext;
out.curr=j_ext/((sum(j_ext)/sqrt(2))^2);
out.X_K=(-1/(1j*w))*(Matrix_Projm*Matrix_G_m*Matrix_Nme);
out.S=S;
out.Rm=Matrix_Rm;
out.f=real(objf);
end