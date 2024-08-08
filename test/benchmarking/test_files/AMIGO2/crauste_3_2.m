addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='crauste_3model'; % Folder to keep results
inputs.pathd.short_name='crauste_3';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=5;                                  % Number of states:\\\
inputs.model.n_par=13;                                 % Number of model parameters
inputs.model.st_names=char('n', 'e', 's', 'm', 'p');    % Names of the states
inputs.model.par_names=char('muN', 'muEE', 'muLE', 'muLL', 'muM', 'muP', 'muPE', 'muPL', 'deltaNE', 'deltaEL', 'deltaLM', 'rhoE', 'rhoP');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char( 'dn = -1 * n * muN - n * p * deltaNE;',  'de = n * p * deltaNE - e * e * muEE - e * deltaEL + e * p * rhoE;',  'ds = s * deltaEL - s * deltaLM - s * s * muLL - e * s * muLE;',  'dm = s * deltaLM - muM * m;',  'dp = p * p * rhoP - p * muP - e * p * muPE - s * p * muPL;');               % Equations describing system dynamics.
inputs.model.par = [0.267, 0.229, 0.622, 0.303, 0.473, 0.296, 0.227, 0.188, 0.625, 0.211, 0.257, 0.395, 0.757];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.178, 0.77, 0.177, 0.881, 0.475];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = n', 'y2 = e', 'y3 = s+m', 'y4 = p');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.1780000000000000 0.7700000000000000 1.0580000000000001 0.4750000000000000
0.1730608647128911 0.7649288899307114 1.0345974519719745 0.4715748776804923
0.1682767568120863 0.7598115420236969 1.0117434563530006 0.4681626648793320
0.1636423171937301 0.7546523885221187 0.9894241349267848 0.4647631191409535
0.1591523884994386 0.7494556356262904 0.9676260101214122 0.4613760111173148
0.1548020068595841 0.7442252742104696 0.9463359910876854 0.4580011242041332
0.1505863940048995 0.7389650900635069 0.9255413603538716 0.4546382541946068
0.1465009497282893 0.7336786736715822 0.9052297610306612 0.4512872089490162
0.1425412446799433 0.7283694295606022 0.8853891845416060 0.4479478080787268
0.1387030134796315 0.7230405852156179 0.8660079588551873 0.4446198826432332
0.1349821481308859 0.7176951995942376 0.8470747371956759 0.4413032748589974
0.1313746917226492 0.7123361712504898 0.8285784872110512 0.4379978378189337
0.1278768324047074 0.7069662460851578 0.8105084805772049 0.4347034352214946
0.1244848976239349 0.7015880247381471 0.7928542830185844 0.4314199411083891
0.1211953486090567 0.6962039696379613 0.7756057447263367 0.4281472396100563
0.1180047750922625 0.6908164117228844 0.7587529911558899 0.4248852246980832
0.1149098902566080 0.6854275568479781 0.7422864141867400 0.4216337999438329
0.1119075258987046 0.6800394918915155 0.7261966636280202 0.4183928782826048
0.1089946277967769 0.6746541905739335 0.7104746390542849 0.4151623817827210
0.1061682512745113 0.6692735190021039 0.6951114819563902 0.4119422414189695
0.1034255569518738 0.6638992409509012 0.6800985681935625 0.4087323968499093
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=2.0*ones(1,13);
inputs.PEsol.global_theta_min=0.0*ones(1,13);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=2.0*ones(1,5);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=0.0*ones(1,5);
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
inputs.nlpsol.eSS.log_var=1:(5+13); 
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
