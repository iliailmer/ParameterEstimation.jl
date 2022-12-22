% TITLE: LV model

addpath(genpath('~/parameter-estimation/matlab'))
addpath(genpath("./"))



%======================

% PATHS RELATED DATA

%======================

inputs.pathd.results_folder='DaisyEx3Model'; % Folder to keep results
inputs.pathd.short_name='Dex3';                 % To identify figures and reports

%======================

% MODEL RELATED DATA

%======================

clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=4;                                  % Number of states
inputs.model.n_par=5;                                 % Number of model parameters
%inputs.model.n_stimulus=0;                            % Number of inputs, stimuli or control variables
inputs.model.st_names=char('x1','x2', 'x3', 'u0');    %x1=V, x2=R        % Names of the states
inputs.model.par_names=char('p1', 'p3', 'p4', 'p6', 'p7');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char('dx1 = -1 * p1 * x1 + x2 + u0;', 'dx2 = p3 * x1 - p4 * x2 + x3;', 'dx3 = p6 * x1 - p7 * x3;', 'du0 = 1;'); % Equations describing system dynamics.
inputs.model.par = [0.2, 0.4, 0.3, 0.2, -0.2];         % Nominal value for the parameters
% inputs.model.AMIGOsensrhs = 1;                       % Generate the sensitivity equations for exact
%==================================

% EXPERIMENTAL SCHEME RELATED DATA

%==================================

% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[.9, 2.0, 2.0, 1];        % Initial conditions
inputs.exps.t_f{1}=10;                       % Experiments duration
inputs.exps.n_obs{1}=2;                       % Number of observables
inputs.exps.obs_names{1}=char('y1','y2'); % Names of the observables
inputs.exps.obs{1}=char('y1=x2','y2=x1+x3');
inputs.exps.t_con{1}=[0 1];                 % Input swithching times including:
inputs.exps.n_s{1}=10;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[2.0 2.0
2.0428356360309623 2.4164313539808657
2.105926682978807 2.867853717076878
2.1912287529021253 3.3577094330748594
2.3008181222902726 3.8897933929380573
2.43691107752443 4.468279643158513
2.6018843914523857 5.097751084304832
2.7982970663480975 5.78323250927625
3.0289134912729128 6.530227253032823
3.2967281771635824 7.344757752723998
];
inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=2.*ones(1,5);
inputs.PEsol.global_theta_min=-1.*ones(1,5);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=2.*ones(1,4);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=-1.*ones(1,4);

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