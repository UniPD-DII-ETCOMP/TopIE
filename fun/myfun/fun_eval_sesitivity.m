function [dfdx] = fun_eval_sesitivity(x,u,lambda,lamDF_Ru,opt)
%%
mu_0=4*pi*1e-7;
alphaSIMP=1;
%%
F_mur=opt.murmin+(opt.murmax-opt.murmin)*(x.^alphaSIMP);
dF_mur=alphaSIMP*(opt.murmax-opt.murmin)*(x.^(alphaSIMP-1));
dparam=-dF_mur./(mu_0*(F_mur-1).^2);
%%
dfdx=lamDF_Ru(dparam,lambda,imag(u));
end