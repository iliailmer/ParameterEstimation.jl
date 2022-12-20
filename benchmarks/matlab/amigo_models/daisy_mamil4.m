% TITLE: LV model

addpath(genpath('~/parameter-estimation/matlab'))
addpath(genpath("./"))



%======================

% PATHS RELATED DATA

%======================

inputs.pathd.results_folder='FHNModel'; % Folder to keep results

inputs.pathd.short_name='FHN';                 % To identify figures and reports



%======================

% MODEL RELATED DATA

%======================

clear



inputs.model.input_model_type='charmodelC';           % Model type- C



inputs.model.n_st=2;                                  % Number of states

inputs.model.n_par=3;                                 % Number of model parameters

%inputs.model.n_stimulus=0;                            % Number of inputs, stimuli or control variables

inputs.model.st_names=char('x1','x2');    %x1=V, x2=R        % Names of the states

inputs.model.par_names=char('a','b','g');             % Names of the parameters

%inputs.model.stimulus_names=char('light');  % Names of the stimuli



inputs.model.eqns=char('dx1 = g * (x1 - x1^3/3 + x2);', 'dx2 = - 1/g * (x1 - a + b * x2);');                                 % Equations describing system dynamics.

                            %Time derivatives are regarded 'd'st_name''





inputs.model.par = [0.2 0.2 2];         % Nominal value for the parameters



% inputs.model.AMIGOsensrhs = 1;                       % Generate the sensitivity equations for exact

%                                                      % Jacobian computation

%==================================

% EXPERIMENTAL SCHEME RELATED DATA

%==================================

% EXPERIMENT DESIGN



inputs.exps.n_exp=1;                          % Number of experiments



% EXPERIMENT 1



inputs.exps.exp_y0{1}=[1 -1];        % Initial conditions

inputs.exps.t_f{1}=5;                       % Experiments duration



inputs.exps.n_obs{1}=1;                       % Number of observables

inputs.exps.obs_names{1}=char('Y'); % Names of the observables

inputs.exps.obs{1}=char('Y=x1');

inputs.exps.t_con{1}=[0 5];                 % Input swithching times including:



inputs.exps.n_s{1}=50;

inputs.exps.data_type='real';

inputs.exps.exp_data{1}=[ 1.0

0.9347004839815006

0.8734370360176208

0.8143346316470316

0.7556774113342264

0.6957586392702363

0.6326920688469742

0.5642544040070117

0.4878097834390997

0.39981783663941944

0.29566680712598137

0.1693528801321386

0.01303853127332005

-0.18250153495546084

-0.4268597226140922

-0.7234527632926807

-1.0602927803013233

-1.4019854874805533

-1.700666684179913

-1.9245447982432382

-2.0716028661983974

-2.166840808939588

-2.231643000262245

-2.2799841389406557

-2.320811552193883

-2.3575047396475557

-2.3925100770746317

-2.4258790955413474

-2.4593497682737318

-2.491449731593281

-2.5238470970160867

-2.555815775576268

-2.5862994927517304

-2.6185068247553063

-2.6489487673251686

-2.6784212728524386

-2.7108562794521878

-2.738899350988676

-2.7693365580748965

-2.8005380383472023

-2.825879263710476

-2.8599427903551993

-2.88492714886106

-2.916265913279319

-2.9441919264887457

-2.97188507713319

-3.0017944253376596

-3.0281950036978156

-3.0579007463380976

-3.083802610662054];



inputs.PEsol.id_global_theta=char('a', 'b', 'g');

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