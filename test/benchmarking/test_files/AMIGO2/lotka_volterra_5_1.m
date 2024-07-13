addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='lotka_volterra_5model'; % Folder to keep results
inputs.pathd.short_name='lotka_volterra_5';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=2;                                  % Number of states:\\\
inputs.model.n_par=3;                                 % Number of model parameters
inputs.model.st_names=char('r', 'w');    % Names of the states
inputs.model.par_names=char('k1', 'k2', 'k3');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char( 'dr = k1*r - k2*r*w;',  'dw = k2*r*w - k3*w;');               % Equations describing system dynamics.
inputs.model.par = [0.612, 0.215, 0.856];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.517, 0.432];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=1;                       % Number of observables
inputs.exps.obs_names{1}=char('y1'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = r');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.5170000000000000
0.5306401871657208
0.5447306711444647
0.5592844100990331
0.5743147953825504
0.5898356626780709
0.6058613036036767
0.6224064777867171
0.6394864254117314
0.6571168802473871
0.6753140831584821
0.6940947961097576
0.7134763166689762
0.7334764930173271
0.7541137394756514
0.7754070525557345
0.7973760275466911
0.8200408756459140
0.8434224416462268
0.8675422221894342
0.8924223845986755
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=1.0*ones(1,3);
inputs.PEsol.global_theta_min=0.0*ones(1,3);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=1.0*ones(1,2);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=0.0*ones(1,2);
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
inputs.nlpsol.eSS.log_var=1:(2+3); 
inputs.nlpsol.eSS.local.solver = 'nl2sol';
inputs.nlpsol.eSS.local.finish = 'nl2sol';
inputs.nlpsol.eSS.maxeval = 200000;                  % Maximum number of cost function evaluations
inputs.nlpsol.eSS.maxtime = 600;                    % Maximum time spent for optimization
inputs.nlpsol.eSS.local.nl2sol.maxiter             =      100000;
inputs.nlpsol.eSS.local.nl2sol.maxfeval            =      100000;
inputs.nlpsol.eSS.local.nl2sol.tolrfun             =     1e-13;
inputs.nlpsol.eSS.local.nl2sol.tolafun             =     1e-13;
inputs.nlpsol.eSS.local.nl2sol.objrtol			 =     1e-13;
% inputs.exps.u_interp{1}='sustained';          % Stimuli definition for experiment 1
                                              % Initial and final time
%inputs.exps.u{1}=1;                           % Values of the inputs for exp 1
AMIGO_Prep(inputs);
[PEresults] = AMIGO_PE(inputs);
PEresults.fit.global_theta_estimated
PEresults.fit.global_theta_y0_estimated
