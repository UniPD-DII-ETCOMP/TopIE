%% --------------------------------------------------------------------- %%
%                           ** TOBS class **                              %
% Modified by LF
%-------------------------------------------------------------------------%

classdef myTOBS
    
    %% Properties
    properties
        
        % TOBS parameters
        epsilons
        flip_limits
        symmetry
        iter
        
        % Design variables {0,1}
        design_variables
        
        % Objective and sensitivities
        objective
        objective_sensitivities
        
        % Constraints and sensitivities
        % [ constraint_1, constraint_2, ..., constraint_n ]
        constraints
        constraints_limits
        constraints_sensitivities
        
        % Optimization history
        % [ objective, constraint_1, constraint_2, ..., constraint_n ]
        history
        
    end
    
    %% Methods
    methods
        
        %% Constructor
        function tobs = myTOBS(constraints_limits_in, epsilons_in, flip_limits_in, number_of_variables,id_filled)
            
%             disp([' '])
%             disp(['         Preparing TOBS.'])
            
            % Optimization constraint limits
            tobs.constraints_limits = constraints_limits_in;
            
            % Optimization parameters
            tobs.epsilons = epsilons_in;
            tobs.flip_limits = flip_limits_in;
            
            % Initial design variables
            tobs.design_variables=zeros(number_of_variables,1);
            %
            tobs.design_variables(id_filled)=1;
            
            % Default symmetry condition
            tobs.symmetry = 0;
            
        end % end Constructor
        
        %% Function to solve optimization problem with ILP
        function tobs = mySolveWithILP(tobs,opt_dir) % added external opt_dir to select optimization type
            
            % Add CPLEX library.
%             addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64');
%             addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\examples\src\matlab');
%             addpath('/opt/ibm/ILOG/CPLEX_Studio1271/cplex/matlab/x86-64_linux');
%             addpath('/opt/ibm/ILOG/CPLEX_Studio1271/cplex/examples/src/');

            % Prepare ILP problem.
            clear options
            options.Symmetry = tobs.symmetry;
%             options.Optimizer = 'cplex';
            options.Optimizer = 'intlinprog';
            COptimize = myILP(tobs.epsilons', tobs.constraints_limits', tobs.constraints', tobs.design_variables,tobs.flip_limits,opt_dir);
            tobs.design_variables = COptimize.Optimize (tobs.objective_sensitivities, tobs.constraints_sensitivities, options);
            
        end % end SolveWithILP
    
    end % end methods
end