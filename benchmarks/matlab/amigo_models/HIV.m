addpath(genpath('~/parameter-estimation/matlab/src'))
addpath(genpath("./"))
% lm   : 1.0001e-01  +-  5.4267e-05 (  0.0543%);
% d    : 1.0007e-01  +-  8.5117e-05 (  0.0851%);
% beta : 9.9941e-02  +-  5.8156e-05 (  0.0582%);
% a    : 9.9953e-02  +-  9.6916e-05 (   0.097%);
% k    : 1.0012e-01  +-  1.0840e-04 (   0.108%);
% uu   : 1.0011e-01  +-  1.5170e-04 (   0.152%);
% c    : 9.9992e-02  +-  2.8461e-05 (  0.0285%);
% q    : 9.9996e-02  +-  3.2260e-04 (   0.323%);
% b    : 9.9990e-02  +-  5.3582e-05 (  0.0536%);
% h    : 1.0000e-01  +-  3.4715e-05 (  0.0347%);
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
inputs.model.par = [0 0 0 0 0 0 0 0 0 0];%[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];         % Nominal value for the parameters
% inputs.model.AMIGOsensrhs = 1;                       % Generate the sensitivity equations for exact
%                                                      % Jacobian computation
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[1.0 1.0 1.0 1.0 1.0];        % Initial conditions
inputs.exps.t_f{1}=10;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('Y1', 'Y2', 'Y3', 'Y4'); % Names of the observables
inputs.exps.obs{1}=char('Y1=w', 'Y2=z', 'Y3=x', 'Y4=yy+vv');
inputs.exps.t_con{1}=[0 10];                 % Input swithching times including:
inputs.exps.n_s{1}=10;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
1.0 1.0 1.0 2.0
0.983119630378905 0.9052545409667376 0.9003738172091765 1.994260678960563
0.9561794490969373 0.8201348065908242 0.8206563764914304 1.9786096329620322
0.9214740875560065 0.7434687927995153 0.7569838596651506 1.9550801756278324
0.8811869741045105 0.6742647669255988 0.7062778662629082 1.925360062027259
0.837233886743163 0.6116812961168062 0.6660663034830178 1.890860719660369
0.7912061006833921 0.5550026435405538 0.6343674652448461 1.8527600966628055
0.7443722393031965 0.5036132711192001 0.6095771764596721 1.812043059644642
0.697706198982586 0.4569811526346315 0.5904022851383369 1.7695263845297464
0.6519335918497772 0.41464177285309733 0.5758023922503137 1.7258830104713447
];
inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=1.*ones(1,10);
inputs.PEsol.global_theta_min=0.0001.*ones(1,10);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=1.1*ones(1,5);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=0.9*ones(1,5);                % Minimum allowed values for the initial conditions
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

