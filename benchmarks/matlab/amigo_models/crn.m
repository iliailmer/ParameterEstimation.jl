addpath(genpath('../src'))
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='CRNModel'; % Folder to keep results
inputs.pathd.short_name='CRN';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=6;                                  % Number of states
inputs.model.n_par=6;                                 % Number of model parameters
%inputs.model.n_stimulus=0;                            % Number of inputs, stimuli or control variables
inputs.model.st_names=char('x1','x2', 'x3', 'x4','x5','x6');    %x1=V, x2=R        % Names of the states
inputs.model.par_names=char('k1', 'k2', 'k3', 'k4', 'k5', 'k6');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char('dx1 = -k1 * x1 * x2 + k2 * x4 + k4 * x6;', 'dx2 = -k1 * x1 * x2 + k2 * x4 + k3 * x4;', 'dx3 = k3 * x4 + k5 * x6 - k6 * x3 * x5;', 'dx4 = k1 * x1 * x2 - k2 * x4 - k3 * x4;', 'dx5 = k4 * x6 + k5 * x6 - k6 * x3 * x5;', 'dx6 = -k4 * x6 - k5 * x6 + k6 * x3 * x5;');                                 % Equations describing system dynamics.
inputs.model.par = [0.03, 0.02, 0.05, 0.03, 0.02, 0.05];         % Nominal value for the parameters
% inputs.model.AMIGOsensrhs = 1;                       % Generate the sensitivity equations for exact
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0];        % Initial conditions
inputs.exps.t_f{1}=30;                       % Experiments duration
inputs.exps.n_obs{1}=2;                       % Number of observables
inputs.exps.obs_names{1}=char('y1','y2'); % Names of the observables
inputs.exps.obs{1}=char('y1=x2','y2=x3');
inputs.exps.t_con{1}=[0 30];                 % Input swithching times including:
inputs.exps.n_s{1}=10;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[1.0 1.0
1.1105980863382912 1.0524829527692867
1.1854409551877314 1.0845025566173467
1.2359391619599875 1.104348979975541
1.270046091741425 1.1169904599859937
1.2931807136842328 1.1253639594375038
1.308986032668638 1.1311917824224922
1.3198928706217408 1.1354791380888003
1.327517268253678 1.138811327070273
1.3329323448933883 1.1415292371942407
];
inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=1.*ones(1,6);
inputs.PEsol.global_theta_min=0.0001.*ones(1,6);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=1.*ones(1,6);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=0.00001*ones(1,6);
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
