function [xbest,fxbest,gbest,...
            stuff,tobs] = TopIE_optimization(fun_solve_EM,...
                                            freq,...
                                            fun_Matrix_Rm,lamDF_Ru,...
                                            Matrix_Ca_m,Matrix_Pm_vv,Matrix_Pm_vs,Matrix_Pm_ss,...
                                            Matrix_Projm,Matrix_G_m,Matrix_Nme,...
                                            F1_m,F1_e,C1_m,Fac_Ed_loc_m,...
                                            Matrix_P0_m,Matrix_P0_e,...
                                            bar_m,Area_vol_m,...
                                            Aed_m,Led_m,Area_m,...
                                            opt,tobs,...
                                            max_it,tau,min_change,...
                                            flag_correct_first,...
                                            flag_plot_geo,flag_plot_conv)
%% Check input arguments
if (nargin<24)
    disp('Error, function requires at least 24 input');
    return;
end
if (nargin < 25 || isempty(max_it))
    max_it=500; % maximum number of iterations
end
if (nargin < 26 || isempty(tau))
    tau=1e-3; % convergence tolerance 1 
end
if (nargin < 27 || isempty(min_change))
    min_change=1e-4; % convergence tolerance 2
end
if (nargin < 28 || isempty(flag_correct_first))
    flag_correct_first=true; % correct first topology according to sensitivity
end
if (nargin < 29 || isempty(flag_plot_geo))
    flag_plot_geo=true; % show layouts during optimization
end
if (nargin < 30 || isempty(flag_plot_geo))
    flag_plot_conv=true; % plot convergence during optimization
end
%% Settings
idx_fig_geo=120;
idx_fig_conv=121;
%% Initialization
flag_vol=0;
loop=0;
x=tobs.design_variables;
nDesign=length(x);
idx_design=1:nDesign;
tobs.history=zeros(max_it,3);
Vtot=sum(Area_vol_m);
%% Prepare density filter
filter_radius=2*(min(Area_vol_m))^1/3;
H=my_prepareDensityFilter(bar_m',filter_radius);
%% Assign material
id_fill=idx_design(x==1);
id_empty=setdiff([1:nDesign],id_fill);
%% Plot initial topology
if flag_plot_geo
    fun_boundary_plot(F1_m,F1_e,Matrix_P0_m,Matrix_P0_e,1,id_fill,id_empty,idx_fig_geo);
    title('Initial topology')
    drawnow
end
%% Solution and objective 
[out]=fun_solve_EM(x);
tobs.objective=out.f;
tobs.objective=-tobs.objective;
tobs.constraints=(x.'*Area_vol_m)/(sum(Area_vol_m));
xbest=x;
fxbest=out.f;
gbest=tobs.constraints;
if (tobs.constraints < 1.05*tobs.constraints_limits && ...
        tobs.constraints > 0.95*tobs.constraints_limits)
    flag_vol=1;
end
%% 
dpfdpu=out.X_K*out.curr;
lambda=out.S.'\(-dpfdpu);
[dfdx]=fun_eval_sesitivity(x,out.u,lambda,lamDF_Ru,opt);
dfdx=-dfdx;
dgdx=Area_vol_m/(sum(Area_vol_m)); 
%% Filtering sensitivities
dfdx=(H*dfdx)./sum(H,2);
tobs.objective_sensitivities=dfdx;
tobs.constraints_sensitivities=dgdx;
%%
sensitivities_previous=tobs.objective_sensitivities;
%% Storing optimization history
tobs.history(loop+1,:)=[abs(tobs.objective) tobs.constraints fxbest];
%% Modify first x by using the first information of sensitivity
if flag_correct_first
    [~,id]=sort(dfdx);
    idf=[];
    vvol=0;
    indii=1;
    while vvol<opt.volfrac*Vtot
        idf=[idf;id(indii)];
        vvol=vvol+Area_vol_m(id(indii));
        indii=indii+1;
    end
    x=zeros(nDesign,1);
    x(idf)=1;
    tobs.design_variables=x;
    xbest=x;
    [out]=fun_solve_EM(xbest);
    %
    tobs.objective=out.f;
    tobs.objective=-tobs.objective;
    tobs.constraints=(x.'*Area_vol_m)/(sum(Area_vol_m));
    dpfdpu=out.X_K*out.curr;
    lambda=out.S.'\(-dpfdpu);
    [dfdx]=fun_eval_sesitivity(x,out.u,lambda,lamDF_Ru,opt);
    dfdx=-dfdx;
    dgdx=Area_vol_m/(sum(Area_vol_m)); 
    dfdx=(H*dfdx)./sum(H,2);
    tobs.objective_sensitivities=dfdx;
    tobs.constraints_sensitivities=dgdx;
    %
    fxbest=out.f;
    gbest=tobs.constraints;
end
%% Plot objective and constraint trend
if flag_plot_conv
    figure(idx_fig_conv)
    grid on
    hold on
    subplot(2,1,1)
    plot(1,tobs.history(1,1),'-dm','linewidth',2);
    hold on
    plot(1,tobs.history(1,3),'-sr','linewidth',2,...
            'MarkerSize',8,'MarkerFaceColor',"#A2142F");
    grid on
    ylabel('Objective function')
    subplot(2,1,2)
    plot(1,tobs.history(1,2),'-ob','linewidth',2,...
            'MarkerSize',8,'MarkerFaceColor','c');
    hold on
    grid on
    yl=yline(opt.volfrac,'-','Vfrac','linewidth',1.5);
    yl.LabelHorizontalAlignment = 'left';
    ylim([0 1])
    ylabel('Volume fraction')
    xlabel('Iterations');
    sgtitle('Convergence trend')
    set(findobj(gcf,'type','axes'),'FontName','Century','FontSize',16);
    drawnow
end
%%
disp([' It.: ',num2str(loop), ' Obj.: ', num2str(abs(tobs.objective)),...
    '  Vol. Frac.: ',num2str(abs(tobs.constraints))])
if flag_correct_first
disp('--- reinitialize with maximum sensitivity variables ---')
disp([' Obj.: ', num2str(abs(fxbest)),'  Vol. Frac.: ',num2str(abs(gbest))])
disp('-------------------------------------------------------')
end
%% Iterations
is_converged=0;
change=1;
while (is_converged==0) && (loop<=max_it) && (change>min_change)
    % Iteration counter update
    loop=loop+1;
    %% Solve sequential approximate subproblem with ILP
    warning off
    tobs=mySolveWithILP(tobs,'Minimize');
    warning on
    x=tobs.design_variables;
    %% Solution of EM problem and objective evaluation
    [out]=fun_solve_EM(x);
    tobs.objective=out.f;
    tobs.objective=-tobs.objective;
    tobs.constraints=(x.'*Area_vol_m)/(sum(Area_vol_m));
    %% Check maximum that satisfied constraint
    if (tobs.constraints < 1.05*tobs.constraints_limits && ...
            tobs.constraints > 0.95*tobs.constraints_limits)
        flag_vol=1;
    else
        flag_vol=0;
    end
    %
    if out.f > fxbest || flag_vol==1
        xbest=x;
        fxbest=out.f;
        gbest=tobs.constraints;
    end
    %% Construct material structure
    id_fill=idx_design(x==1);
    id_empty=setdiff([1:nDesign],id_fill);
    %% Plot current topology
    if flag_plot_geo
        clf(idx_fig_geo);
        fun_boundary_plot(F1_m,F1_e,Matrix_P0_m,Matrix_P0_e,1,id_fill,id_empty,idx_fig_geo);
        title(['Topology iteration ',num2str(loop)])
        drawnow
    end
    %%
    dpfdpu=out.X_K*out.curr;
    lambda=out.S.'\(-dpfdpu);
    [dfdx]=fun_eval_sesitivity(x,out.u,lambda,lamDF_Ru,opt);
    dfdx=-dfdx;
    dgdx=Area_vol_m/(sum(Area_vol_m));   
    %% Filtering sensitivities
    dfdx=(H*dfdx)./sum(H,2);
    tobs.objective_sensitivities=dfdx;
    tobs.constraints_sensitivities=dgdx;
    %% Stabilization technique (average of the sensitivities history)
    tobs.objective_sensitivities=(tobs.objective_sensitivities+sensitivities_previous)/2;
    sensitivities_previous=tobs.objective_sensitivities;
    %% Storing optimization history
    tobs.history(loop+1,:)=[abs(tobs.objective) tobs.constraints fxbest];
    %%
    disp([' It.: ',num2str(loop), '  Obj.: ', num2str(abs(tobs.objective)),'  Vol. Frac.: ',num2str(abs(tobs.constraints))])
    %% Plot objective and constraint trend
    if flag_plot_conv
        figure(idx_fig_conv)
        subplot(2,1,1)
        plot(1:loop,tobs.history(1:loop,1),'-dm','linewidth',2);
        hold on
        plot(1:loop,tobs.history(1:loop,3),'-sr','linewidth',2,...
            'MarkerSize',8,'MarkerFaceColor',"#A2142F");
        grid on
        ylabel('Objective function')
        subplot(2,1,2)
        plot(1:loop,tobs.history(1:loop,2),'-ob','linewidth',2,...
            'MarkerSize',8,'MarkerFaceColor','c');
        hold on
        grid on
        yl=yline(opt.volfrac,'-','Vfrac','linewidth',1.5);
        yl.LabelHorizontalAlignment = 'left';
        ylim([0 1])
        ylabel('Volume fraction')
        xlabel('Iterations');
        sgtitle('Convergence trend')
        set(findobj(gcf,'type','axes'),'FontName','Century','FontSize',16);
        drawnow
    end
    %% Convergence analysis TOBS 101
    if loop>10
        change=abs(sum(tobs.history(loop-9:loop-5,1))-sum(tobs.history(loop-4:loop,1)))...
            /sum(tobs.history(loop-4:loop,1));
        if change<=min_change
            disp('-----------------------------')
            disp(' - exit for change <= change_min')
            disp('-----------------------------')
        end
    end    
    %% Convergence analysis [Huang and Xie, 2007]
    N=3;
    if (loop>=2*N) % Analyzing 2*N consecutive iterations
        error_1=zeros(N,1);
        error_2=zeros(N,1);
        for i=1:N
            error_1(i)=tobs.history(loop-i+1,1)-tobs.history(loop-N-i+1,1);
            error_2(i)=tobs.history(loop-i+1,1);
        end
        % Evaluating convergence
        difference=abs(sum(error_1))/abs(sum(error_2));
        % Verifying error tolerance and if constraint is satisfied
        if ((difference <= tau) && (tobs.history(loop,2) <= 1.002*tobs.constraints_limits))
            is_converged=1;
            disp('-----------------------------')
            disp(' - exit for Huang and Xie convergence')
            disp('-----------------------------')
        end
    end
end
%%
stuff.iter=loop;
disp('-----------------------------')
disp([' - TopIE terminated at iteration : ',num2str(loop)])
disp([' - Objective value : ',num2str(fxbest)])
disp([' - Constraint value : ',num2str(gbest)])
disp('-----------------------------')
end