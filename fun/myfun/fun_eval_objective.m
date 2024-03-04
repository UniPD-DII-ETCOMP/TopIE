function [f] = fun_eval_objective(u,Matrix_Projm,Matrix_G_m,Matrix_Nme,jext,freq)
w=2*pi*freq;
a=(-1/(1j*w))*(Matrix_Projm*Matrix_G_m*Matrix_Nme).'*u;
Wm=0.5*jext'*a/((sum(jext)/sqrt(2))^2);
f=Wm;
end