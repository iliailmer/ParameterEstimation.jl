% TITLE: LV model

addpath(genpath('~/parameter-estimation/matlab'))
addpath(genpath("./"))



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

inputs.model.st_names=char('x1','x2');    %x1=V, x2=R        % Names of the states

inputs.model.par_names=char('N', 'E', 'S', 'M', 'P');             % Names of the parameters

%inputs.model.stimulus_names=char('light');  % Names of the stimuli


% Equations describing system dynamics.
inputs.model.eqns=char('dN = -1 * N * mu_N - N * P * delta_NE;', 'dE = N * P * delta_NE - E^2 * mu_EE - E * delta_EL + E * P * rho_E;', 'dS = S * delta_EL - S * delta_LM - S^2 * mu_LL - E * S * mu_LE;', 'dM = S * delta_LM - mu_M * M;', 'dP = P^2 * rho_P - P * mu_P - E * P * mu_PE - S * P * mu_PL;');


inputs.model.par = [1 1 1 1 1 1 0 1 1 1 1 0 1];         % Nominal value for the parameters



% inputs.model.AMIGOsensrhs = 1;                       % Generate the sensitivity equations for exact

%                                                      % Jacobian computation

%==================================

% EXPERIMENTAL SCHEME RELATED DATA

%==================================

% EXPERIMENT DESIGN



inputs.exps.n_exp=1;                          % Number of experiments



% EXPERIMENT 1



inputs.exps.exp_y0{1}=[1.0 1.0 1.0 1.0 1.0];        % Initial conditions

inputs.exps.t_f{1}=5;                       % Experiments duration

inputs.exps.n_obs{1}=4;                       % Number of observables

% Names of the observables
inputs.exps.obs_names{1}=char('Y1', 'Y2', 'Y3', 'Y4');

inputs.exps.obs{1}=char('Y1=N', 'Y2 = E', 'Y3 = S + M', 'Y4 = P');

inputs.exps.t_con{1}=[0 5];                 % Input swithching times including:



inputs.exps.n_s{1}=10;

inputs.exps.data_type='real';

inputs.exps.exp_data{1}=[ 1.0 1.0 2.0 1.0
0.6262864365228055 0.382248281083486 1.1899703139289581 0.5264752527323943
0.35094171626208187 0.17593999852593248 0.7972736508483186 0.2904691769256627
0.19641206134424308 0.08939843309484297 0.5658887302539644 0.15953408533717922
0.1106228295084407 0.04799460408906465 0.4237780442067244 0.086898602754068
0.06265525487541022 0.026560233929535878 0.3344258308327127 0.047131290426414334
0.03563734358156775 0.01494357325016453 0.2767040722680609 0.025543666796321363
0.020330872307974358 0.008483388626826653 0.2380950449053304 0.013862016170187133
0.011622383867554731 0.004839385966218741 0.21116003159716396 0.007539702141978074
0.00665306373325196 0.0027679089960299676 0.19147479918617955 0.0041115921049275994
];



inputs.PEsol.id_global_theta=char( 'mu_N', 'mu_EE', 'mu_LE', 'mu_LL', 'mu_M', 'mu_P', 'mu_PE', 'mu_PL', 'delta_NE', 'delta_EL', 'delta_LM', 'rho_E', 'rho_P',);

inputs.PEsol.global_theta_max=2.*ones(1,13);

inputs.PEsol.global_theta_min=0.0001.*ones(1,13);

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



% % JRBanga: use global solver eSS to avoid converge to local solutions

 inputs.nlpsol.nlpsolver='eSS';                      % Solver used for optimization

 inputs.nlpsol.eSS.log_var = 1:3;                    % Index of parameters to be considered in log scale

 inputs.nlpsol.eSS.maxeval = 20000;                  % Maximum number of cost function evaluations

 inputs.nlpsol.eSS.maxtime = 600;                    % Maximum time spent for optimization

 inputs.nlpsol.eSS.local.solver = 'nl2sol';

 inputs.nlpsol.eSS.local.finish = 'nl2sol';

%  inputs.nlpsol.eSS.local.nl2sol.maxiter = 150;       % Parameters for local solver

%  inputs.nlpsol.eSS.local.nl2sol.maxfeval = 200;

  inputs.nlpsol.eSS.local.nl2sol.display = 1;

%  inputs.nlpsol.eSS.local.nl2sol.objrtol = 1e-6;

%  inputs.nlpsol.eSS.local.nl2sol.tolrfun = 1e-5;

% % JRBanga





% inputs.exps.u_interp{1}='sustained';          % Stimuli definition for experiment 1

                                              % Initial and final time

%inputs.exps.u{1}=1;                           % Values of the inputs for exp 1

AMIGO_Prep(inputs);

AMIGO_PE(inputs);



% % Long experiment; multiple data; low experimental error

% tf=5;    % Final time for plotting

% texp=5;  % Experiment duration; steady state of S=N1+N2

% ns=5;     % Number of sampling times

%

% initial_conditions=[100 100];

% inputs.exps.exp_y0=initial_conditions;

% inputs.exps.n_s=ns;

% inputs.exps.t_f=tf;

% inputs.exps.t_s=linspace(0,texp,ns);

%

% results=AMIGO_SModel(inputs);

%

% states{icase}=results.sim.states{1};

% tsim{icase}=results.sim.tsim{1};

% inputs.exps.std_dev{1}=percentage_error_data/100;

% inputs.exps.t_f{1}=tf;

% results=AMIGO_SData(inputs);

%

%  inputs.exps.n_exp=1;                                % Number of experiments

%  for iexp=1:inputs.exps.n_exp

%  inputs.exps.exp_y0{iexp}=[0.1 0.2 2.5];             % Initial conditions for each experiment

%  inputs.exps.t_f{iexp}=240;                          % Experiments duration

%  inputs.exps.t_con{iexp}=[0 240];                    % Input swithching times: Initial and final time

%  inputs.exps.ts_0{iexp}=1;                           % First sampling time

%  inputs.exps.n_s{iexp}=10;                           % Number of sampling times by default

%  inputs.exps.t_s{iexp}=linspace(1,240,10);           % Sampling times

%

% % OBSEVABLES DEFINITION

%  inputs.exps.n_obs{iexp}=2;                          % Number of observed quantities per experiment

%  inputs.exps.obs_names{iexp}=char('Xobs','Zobs');    % Name of the observed quantities per experiment

%  inputs.exps.obs{iexp}=char('Xobs=x1' ,'Zobs=x3');   % Observation function

%  end

%

%

%

% %==================================

% % EXPERIMENTAL DATA RELATED INFO

% %==================================

%

%  inputs.exps.data_type='pseudo';                     % Type of experimental data:

%                                                      % 'real'|'pseudo'|'pseudo_pos'(>=0)

% inputs.exps.noise_type = 'hetero_lin';

% inputs.exps.std_deva{1} = [0.003 0.1 ];

% inputs.exps.std_devb{1} = [0.1 0.1];                  % Type of experimental noise:

%                                                      % noise scales with the signal, with a low threshold

%

% %==================================

% % UNKNOWNS RELATED DATA

% %==================================

%

% % GLOBAL UNKNOWNS (SAME VALUE FOR ALL EXPERIMENTS)

% inputs.PEsol.id_global_theta='all';                 %  'all'|User selected

% inputs.PEsol.global_theta_max=1e3*ones(1,7);        % Maximum allowed values for the paramters

% inputs.PEsol.global_theta_max(8) = 12;              % 8th parameter is an exponent

% inputs.PEsol.global_theta_min=[1e-3*ones(1,7) 1];   % Minimum allowed values for the paramters

% inputs.PEsol.global_theta_guess= ...                % Initial guess

% AMIGO_logUniformIntialGuess(inputs.PEsol.global_theta_min,inputs.PEsol.global_theta_max);

%

%

%

% %==================================

% % COST FUNCTION RELATED DATA

% %==================================

% inputs.PEsol.PEcost_type='llk';                     %  'llk' (log likelihood)

% inputs.PEsol.llk_type='hetero';

% inputs.PEsol.PEcostJac_type = 'llk';                % AMIGO2 computes high quality Jacobian based on

%                                                     % exact sensitivities

%

%

%

% %==================================

% % NUMERICAL METHODS

% %==================================

%

%

% %

% % OPTIMIZATION

% %

% inputs.nlpsol.nlpsolver='eSS';                      % Solver used for optimization

% inputs.nlpsol.eSS.log_var = 1:7;                    % Index of parameters to be considered in log scale

% inputs.nlpsol.eSS.maxeval = 20000;                  % Maximum number of cost function evaluations

% inputs.nlpsol.eSS.maxtime = 600;                    % Maximum time spent for optimization

% inputs.nlpsol.eSS.local.solver = 'nl2sol';

% inputs.nlpsol.eSS.local.finish = 'nl2sol';

% inputs.nlpsol.eSS.local.nl2sol.maxiter = 150;       % Parameters for local solver

% inputs.nlpsol.eSS.local.nl2sol.maxfeval = 200;

% inputs.nlpsol.eSS.local.nl2sol.display = 2;

% inputs.nlpsol.eSS.local.nl2sol.objrtol = 1e-6;

% inputs.nlpsol.eSS.local.nl2sol.tolrfun = 1e-5;

%

% %======================

% % DISPLAY OF RESULTS

% %======================

%

% inputs.plotd.plotlevel='noplot'; % supress plots, since Regularization will call AMIGO_PE several times

%

%

%

%

% %======================

% % REGULARIZATION

% %======================

% inputs.nlpsol.regularization.ison = 1;

% inputs.nlpsol.regularization.method = 'tikhonov';     % Regularization approach

% inputs.nlpsol.regularization.tikhonov.gW = eye(8);    % Weighting matrix for global parameters,

%                                                       % use lW{iexp} for

%                                                       % experiment-wise local parameters if defined

% inputs.nlpsol.regularization.tikhonov.gx0 = ...

%                      inputs.PEsol.global_theta_min;   % Reference parameters for global unknowns

%                                                       % (prior knowledge). Default: lowest bound

%                                                       % lx0{iexp} for experiment-wise local parameters

% alpha1 = 5e-3;

% inputs.nlpsol.regularization.alphaSet = ...

%     logspace(log10(alpha1*100), log10(alpha1/20),7);  % Set of regularization parameters to be evaluated

%                                                       % by generalized cross validation

%

%

% % Save a copy of inputs for later, for crossvalidation and for regularization-free

% % estimation

% copy_inputs = inputs;