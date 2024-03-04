%% TopIE
close all
clear global
clear
clc
restoredefaultpath
format short
warning off
delete log_lamdRu_axi.txt 
delete logL_axi.txt 
delete logL_axi_ss.txt
delete logRrz_axi.txt
warning on
dad=pwd;
cd('fun'); addpath(genpath(pwd)); cd(dad)
cd('fortran'); addpath(pwd); cd(dad)
%%
%% BEGIN USER SETTINGS
%%
%% Directory
test_case_dir.mag='test1'; % Source directory for magnetic media
test_case_dir.ext='test1'; % Source directory for external coil
%% Current density in the external coil
J_ext=@(r,phi,z) [0,1/7.5e-4,0];
%% Frequency
freq=85e3; 
%% Optimization settings
opt.murmax=100; % Maximum relative permeability
opt.murmin=1.01; % Minumum relative permeability
opt.volfrac=0.6; % Final volume fraction
%% TOBS settings
epsilons=0.008; % Constraint relaxation parameter 
flip_limits=0.01; % Flip limits parameter (beta)
%% Settings
opt.nThread=22; % Number of threads
opt.nGauss=2; % Gauss points numerical integration
%% Flag
show_geometry_flag=true;
%%
%% END USER SETTINGS
%%
%% Checks
if opt.murmin<=1
   disp('murmin must be greater than one!')
   opt.murmin=1.01;
   disp('murmin set to 1.01');
end
if epsilons>flip_limits
   disp('epsilon must be less than flip_limits!')
   epsilon=flip_limits/2;
   disp('epsilon set to flip_limits/2');
end
%% Preprocessing of data
disp('======================================')
disp('Pre-processing...')
mytic_pre=tic;
[N,Matrix_P0_m,Matrix_P0_e,F1_m,F1_e,C1_m,...
 Fac_Ed_loc_m,Aed_m,Led_m,Area_m,Area_e,Area_vol_m,...
 ind_edge_free_mag,...
 bar_e,bar_m,Matrix_Cb_m,Matrix_Ca_m,...
 G1_m,Matrix_G_m,Matrix_C_m]=fun_pre_processing(test_case_dir);
mytoc_pre=toc(mytic_pre);
disp(['Time for preprocessing : ' ,num2str(mytoc_pre),' s']);
disp(' ')
%% Boundary plot
if show_geometry_flag
    fun_boundary_plot(F1_m,F1_e,Matrix_P0_m,Matrix_P0_e,1,1:size(F1_m,2),[],1);
end
%% Computing matrices
disp('======================================')
disp('Computing matrices...')
mytic_mat=tic;
[Matrix_Pm_vv,Matrix_Pm_vs,Matrix_Pm_ss,...
 Matrix_Projm,Matrix_Nme]=fun_compute_matrices(N,Matrix_P0_m,...
                                                Matrix_P0_e,...
                                                G1_m,F1_m,F1_e,C1_m,...
                                                Fac_Ed_loc_m,Aed_m,Led_m,...
                                                Area_e,Area_m,Area_vol_m,...
                                                ind_edge_free_mag,...
                                                opt.nGauss,opt.nThread);
mytoc_mat=toc(mytic_mat);
disp(['Time for computing matrices : ' ,num2str(mytoc_mat),' s']);
disp(' ')
%% External field
[h_ext,j_ext]=fun_compute_ext_field(N,Matrix_Projm,Matrix_G_m,Matrix_Nme,bar_e,Area_e,J_ext);
%% Rm(rho_m) constructor
fun_Matrix_Rm=@(x) fun_construct_Rm(x,N,Matrix_P0_m,...
                                F1_m,C1_m,Fac_Ed_loc_m,...
                                Aed_m,Led_m,...
                                opt.nGauss,opt.nThread);
%% dRm(rho_m)/drho constructor
lamDF_Ru=@(x,lambda,u) fun_vec_lamDRu(x,lambda,u,...
                                      N,Matrix_P0_m,F1_m,...
                                      C1_m,Fac_Ed_loc_m,...
                                      Aed_m,Led_m,...
                                      opt.nThread,...
                                      opt.nGauss);
%% Solver
fun_solve_EM=@(x) fun_solve_system(x,...
                                freq,...
                                fun_Matrix_Rm,...
                                Matrix_Ca_m,...
                                Matrix_Pm_vv,Matrix_Pm_vs,Matrix_Pm_ss,...
                                Matrix_Projm,Matrix_G_m,Matrix_Nme,...
                                h_ext,j_ext,opt);
%% TOBS settings
nDesign=size(F1_m,2);
id_fill=1:nDesign;
tobs=myTOBS(opt.volfrac,epsilons,flip_limits,nDesign,id_fill);
%%
disp('======================================')
disp('Running TopIE...')
mytic_opt=tic; 
[xbest,fxbest,gbest,stuff,tobs]=TopIE_optimization(fun_solve_EM,...
                                                  freq,...
                                                  fun_Matrix_Rm,lamDF_Ru,...
                                                  Matrix_Ca_m,Matrix_Pm_vv,Matrix_Pm_vs,Matrix_Pm_ss,...
                                                  Matrix_Projm,Matrix_G_m,Matrix_Nme,...
                                                  F1_m,F1_e,C1_m,Fac_Ed_loc_m,...
                                                  Matrix_P0_m,Matrix_P0_e,...
                                                  bar_m,Area_vol_m,...
                                                  Aed_m,Led_m,Area_m,...
                                                  opt,tobs);
mytoc_opt=toc(mytic_opt);
disp(['Time for computing matrices : ' ,num2str(mytoc_opt),' s']);
disp(' ')
%% 
id_fill=find(xbest);
id_empty=setdiff([1:nDesign],id_fill);
fun_boundary_plot(F1_m,F1_e,Matrix_P0_m,Matrix_P0_e,1,...
                   id_fill,id_empty,2)
title('Final topology')