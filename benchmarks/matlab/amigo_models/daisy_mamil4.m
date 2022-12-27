addpath(genpath('~/parameter-estimation/matlab/src'))
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
inputs.model.n_st=4;                                  % Number of states
inputs.model.n_par=7;                                 % Number of model parameters
%inputs.model.n_stimulus=0;                            % Number of inputs, stimuli or control variables
inputs.model.st_names=char('x1','x2','x3', 'x4');    %x1=V, x2=R        % Names of the states
inputs.model.par_names=char('k01', 'k12', 'k13', 'k14', 'k21', 'k31', 'k41');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char('dx1 = -k01 * x1 + k12 * x2 + k13 * x3 + k14 * x4 - k21 * x1 - k31 * x1 - k41 * x1;', 'dx2 = -k12 * x2 + k21 * x1;', 'dx3 = -k13 * x3 + k31 * x1;', 'dx4 = -k14 * x4 + k41 * x1;');                                 % Equations describing system dynamics.
inputs.model.par = [0.18, 0.21, 0.58, 0.64, -0.12, 0.9, 0.02];         % Nominal value for the parameters
% inputs.model.AMIGOsensrhs = 1;                       % Generate the sensitivity equations for exact
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[1.0, 2.2, 0.8, -1.0];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1','y2','y3','y4'); % Names of the observables
inputs.exps.obs{1}=char('y1=x1+x2', 'y2=x3+x4', 'y3=x1+x3', 'y4=x2+x4');
inputs.exps.t_con{1}=[0 1];                 % Input swithching times including:
inputs.exps.n_s{1}=10;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[3.0 0.0 2.0 1.0
2.8509387155098813 0.12755528760039414 1.998636782263437 0.9798572208468387
2.7175361372758493 0.24070462529903613 1.9988665735630926 0.9593741890117928
2.597454684726313 0.3415507403118677 2.0005115513803133 0.9384938736578673
2.4887322192333285 0.43186105881012404 2.0034063316093556 0.917186946434097
2.3897211709442616 0.5131219811627882 2.007397375496862 0.8954457766101879
2.299037536094468 0.5865843600527896 2.0123423903584623 0.8732795057887957
2.215518145791263 0.6533016064344891 2.018109731973423 0.8507100202523289
2.1381848652579523 0.7141616183919514 2.0245778148140325 0.8277686688358709
2.0662146005199458 0.7699135336974288 2.031634534856743 0.8044935993606319
];
inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=2.*ones(1,7);
inputs.PEsol.global_theta_min=-1.*ones(1,7);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=2.*ones(1,4);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=-1*ones(1,4)
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
%  inputs.nlpsol.eSS.local.nl2sol.maxiter = 1500;       % Parameters for local solver
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
