addpath(genpath('../src'))
addpath(genpath("./"))
addpath(genpath("../../matlab_tol/samples/hiv/"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='HIVModel'; % Folder to keep results
inputs.pathd.short_name='HIV';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=5;                                  % Number of states
inputs.model.n_par=10;                                 % Number of model parameters
%inputs.model.n_stimulus=0;                            % Number of inputs, stimuli or control variables
inputs.model.st_names=char('x','yy','vv','w','z');           % Names of the states
inputs.model.par_names=char('lm', 'd', 'beta', 'a', 'k', 'uu', 'c', 'q', 'b', 'h');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char('dx = lm - d * x - beta * x * vv;','dyy = beta * x * vv - a * yy;','dvv = k * yy - uu * vv;','dw = c * x * yy * w - c * q * yy * w - b * w;','dz = c * q * yy * w - h * z;');                                 % Equations describing system dynamics.
                            %Time derivatives are regarded 'd'st_name''
inputs.model.par = [0.1 0.11 0.1 0.1 0.09 0.1 0.1 0.1 0.1 0.1];         % Nominal value for the parameters
% inputs.model.AMIGOsensrhs = 1;                       % Generate the sensitivity equations for exact
%                                                      % Jacobian computation
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[1.0 1.0 1.0 1.0 1.0];        % Initial conditions
inputs.exps.t_f{1}=20;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('Y1', 'Y2', 'Y3', 'Y4'); % Names of the observables
inputs.exps.obs{1}=char('Y1=x', 'Y2=z', 'Y3=w', 'Y4=yy+vv');
inputs.exps.t_con{1}=[0 10];                 % Input swithching times including:
inputs.exps.n_s{1}=10;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=readmatrix(sprintf('hiv-%i.csv', inputs.exps.n_s{1})); % read sample data
% inputs.exps.t_s{1} = inputs.exps.exp_data{1}(:, 1);
inputs.exps.exp_data{1}(:, 1) = [];
inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=1.*ones(1,10);
inputs.PEsol.global_theta_min=0.0001.*ones(1,10);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=2*ones(1,5);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=0.001*ones(1,5);                % Minimum allowed values for the initial conditions
%=============================================================
% COST FUNCTION RELATED DATA
% SOLVING THE PROBLEM WITH WEIGHTED LEAST SQUARES FUNCTION
%=============================================================
inputs.PEsol.PEcost_type='lsq';          % 'lsq' (weighted least squares default)
inputs.PEsol.lsq_type='Q_I';             % Weights:
                                         % Q_I: identity matrix; Q_expmax: maximum experimental data
                                         % Q_expmean: mean experimental data;
                                         % Q_mat: user selected weighting matrix
% OPTIMIZATION
%inputs.nlpsol.nlpsolver='local_lsqnonlin';  % In this case the problem will be solved with
                                         % a local non linear least squares
                                         % method.AMIGO_Prep(inputs);
% %
 inputs.nlpsol.nlpsolver='eSS';                      % Solver used for optimization
%  inputs.nlpsol.eSS.log_var = 1:3;                    % Index of parameters to be considered in log scale
 inputs.nlpsol.eSS.maxeval = 20000;                  % Maximum number of cost function evaluations
 inputs.nlpsol.eSS.maxtime = 600;                    % Maximum time spent for optimization
 inputs.nlpsol.eSS.local.solver = 'nl2sol';
 inputs.nlpsol.eSS.local.finish = 'nl2sol';
%  inputs.nlpsol.eSS.local.nl2sol.maxiter = 150;       % Parameters for local solver
%  inputs.nlpsol.eSS.local.nl2sol.maxfeval = 200;
  inputs.nlpsol.eSS.local.nl2sol.display = 1;
%  inputs.nlpsol.eSS.local.nl2sol.objrtol = 1e-6;
%  inputs.nlpsol.eSS.local.nl2sol.tolrfun = 1e-5;
% %
% inputs.exps.u_interp{1}='sustained';          % Stimuli definition for experiment 1
                                              % Initial and final time
%inputs.exps.u{1}=1;                           % Values of the inputs for exp 1
inputs.ivpsol.rtol=1.0e-10;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-10;
AMIGO_Prep(inputs);
AMIGO_PE(inputs);

