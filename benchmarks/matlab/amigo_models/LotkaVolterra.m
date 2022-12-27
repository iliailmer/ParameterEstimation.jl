addpath(genpath('~/parameter-estimation/matlab/src'))
addpath(genpath("./"))
% k1 : 1.8322e-02  +-  1.6087e-03 (    8.78%);
% k2 : 2.9989e-02  +-  1.2393e-05 (  0.0413%);
% k3 : 4.9062e-02  +-  9.8952e-04 (    2.02%);
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='LVModel'; % Folder to keep results
inputs.pathd.short_name='LV';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=2;                                  % Number of states
inputs.model.n_par=3;                                 % Number of model parameters
%inputs.model.n_stimulus=0;                            % Number of inputs, stimuli or control variables
inputs.model.st_names=char('r','w');           % Names of the states
inputs.model.par_names=char('k1','k2','k3');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char('dr = k1*r-k2*r*w;', 'dw = k2*r*w-k3*w;');                                 % Equations describing system dynamics.
                            %Time derivatives are regarded 'd'st_name''
inputs.model.par = [0.02 0.03 0.05];         % Nominal value for the parameters
% inputs.model.AMIGOsensrhs = 1;                       % Generate the sensitivity equations for exact
%                                                      % Jacobian computation
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[100 100];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=1;                       % Number of observables
inputs.exps.obs_names{1}=char('Y'); % Names of the observables
inputs.exps.obs{1}=char('Y=r');
inputs.exps.t_con{1}=[0 1];                 % Input swithching times including:
inputs.exps.n_s{1}=16; % 8 or 16
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[  100.0
80.39902767817165
62.2690571399069
46.6561358533665
34.013260966178855
24.27292673393346
17.050518686913666
11.843883191532147
8.164810029856213
5.600773144479524
3.8302834137534165
2.6150536131895934
1.7840373850972975
1.2169728866377594
0.8304294265528263
0.5670214463059633];
inputs.PEsol.id_global_theta=char('k1', 'k2', 'k3');
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=[110 110];                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=[90 90];
inputs.PEsol.global_theta_max=1.*ones(1,3);
inputs.PEsol.global_theta_min=0.0001.*ones(1,3);
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
AMIGO_Prep(inputs);
AMIGO_PE(inputs);
