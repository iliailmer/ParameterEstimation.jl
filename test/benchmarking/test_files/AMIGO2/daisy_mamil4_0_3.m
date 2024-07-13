addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='daisy_mamil4_0model'; % Folder to keep results
inputs.pathd.short_name='daisy_mamil4_0';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=4;                                  % Number of states:\\\
inputs.model.n_par=7;                                 % Number of model parameters
inputs.model.st_names=char('x1', 'x2', 'x3', 'x4');    % Names of the states
inputs.model.par_names=char('k01', 'k12', 'k13', 'k14', 'k21', 'k31', 'k41');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char( 'dx1 = -k01 * x1 + k12 * x2 + k13 * x3 + k14 * x4 - k21 * x1 - k31 * x1 - k41 * x1;',  'dx2 = -k12 * x2 + k21 * x1;',  'dx3 = -k13 * x3 + k31 * x1;',  'dx4 = -k14 * x4 + k41 * x1;');               % Equations describing system dynamics.
inputs.model.par = [0.539, 0.672, 0.582, 0.536, 0.439, 0.617, 0.45];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.813, 0.871, 0.407, 0.733];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=3;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = x1', 'y2 = x2', 'y3 = x3 + x4');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.8129999999999999 0.8710000000000000 1.1399999999999999
0.7917061406553852 0.8595345629227156 1.1511406408365945
0.7724068864637377 0.8480104182539394 1.1608918253003155
0.7548715730005824 0.8364700046344019 1.1693911976480795
0.7388966913459911 0.8249496803403117 1.1767611698264355
0.7243026836420433 0.8134804748481188 1.1831106865922434
0.7109311169762184 0.8020887505886896 1.1885367832731260
0.6986421909181093 0.7907967855298956 1.1931259606198483
0.6873125393129681 0.7796232859689999 1.1969553983158909
0.6768332915867556 0.7685838378076034 1.2000940261644597
0.6671083629205073 0.7576913036039562 1.2026034697280523
0.6580529462704464 0.7469561718349872 1.2045388852155914
0.6495921824004303 0.7363868640401109 1.2059496966665089
0.6416599869086006 0.7259900048479465 1.2068802469406139
0.6341980157110683 0.7157706592958535 1.2073703726649281
0.6271547526349489 0.7057325413304676 1.2074559120905315
0.6204847047031583 0.6958781969175099 1.2071691537562170
0.6141476923959014 0.6862091647835118 1.2065392329240201
0.6081082236749156 0.6767261174544492 1.2055924819301631
0.6023349418812995 0.6674289849407027 1.2043527398699136
0.5968001387845273 0.6583170631398092 1.2028416263962567
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=3.0*ones(1,7);
inputs.PEsol.global_theta_min=0.0*ones(1,7);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=3.0*ones(1,4);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=0.0*ones(1,4);
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
inputs.nlpsol.eSS.log_var=1:(4+7); 
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
