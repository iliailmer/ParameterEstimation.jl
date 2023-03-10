addpath(genpath('../src'))
addpath(genpath("./"))
addpath(genpath("../../matlab_tol/samples/crauste/"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='CrausteModel'; % Folder to keep results
inputs.pathd.short_name='Crauste';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=5;                                  % Number of states
inputs.model.n_par=13;                                 % Number of model parameters
%inputs.model.n_stimulus=0;                            % Number of inputs, stimuli or control variables
inputs.model.st_names=char('n', 'e', 's', 'm', 'p');    %x1=V, x2=R        % Names of the states
inputs.model.par_names=char('muN', 'muEE', 'muLE', 'muLL', 'muM', 'muP', 'muPE', 'muPL', 'deltaNE', 'deltaEL', 'deltaLM', 'rhoE', 'rhoP');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
% Equations describing system dynamics.
inputs.model.eqns=char('dn = -1 * n * muN - n * p * deltaNE;', 'de = n * p * deltaNE - e * e * muEE - e * deltaEL + e * p * rhoE;', 'ds = s * deltaEL - s * deltaLM - s * s * muLL - e * s * muLE;', 'dm = s * deltaLM - muM * m;','dp = p * p * rhoP - p * muP - e * p * muPE - s * p * muPL;');
inputs.model.par = [1 1.0 1.1 1.2 1.1 1 0.5 1.0 1.0 1.0 1.0 0.9 1.2]; %[1 1.3 1.1 1.2 1.1 1 0.5 1.0 1.0 1.0 1.0 0.9 1.2];         % Nominal value for the parameters
% inputs.model.AMIGOsensrhs = 1;                       % Generate the sensitivity equations for exact
%                                                      % Jacobian computation
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[1.0 1.0 1.0 1.0 1.0];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
% Names of the observables
inputs.exps.obs_names{1}=char('Y1', 'Y2', 'Y3', 'Y4');
inputs.exps.obs{1}=char('Y1=n', 'Y2=e', 'Y3=s+m', 'Y4=p');
inputs.exps.t_con{1}=[0 1];                 % Input swithching times: Initial and final time
inputs.exps.n_s{1}=10;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=readmatrix(sprintf('crauste-%i.csv', inputs.exps.n_s{1})); % read sample data
% inputs.exps.t_s{1} = inputs.exps.exp_data{1}(:, 1);
inputs.exps.exp_data{1}(:, 1) = [];
display(inputs.exps.exp_data{1});
% display(inputs.exps.exp_data{1});
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
% inputs.PEsol.global_theta_y0_guess=[1.0 1.0 1.0 1.0 1.0];              % [] Initial guess
inputs.PEsol.global_theta_y0_max=1.1 * ones(1,5);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=-1.1 * ones(1,5);
inputs.PEsol.id_global_theta='all';
% inputs.PEsol.global_theta_guess= [1 1.3 1.1 1.2 1.1 1 0.5 1.0 1.0 1.0 1.0 0.9 1.2];
inputs.PEsol.global_theta_max=2.*ones(1,13);
inputs.PEsol.global_theta_min=-2.*ones(1,13);
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
