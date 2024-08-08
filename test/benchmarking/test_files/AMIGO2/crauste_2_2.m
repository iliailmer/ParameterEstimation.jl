addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='crauste_2model'; % Folder to keep results
inputs.pathd.short_name='crauste_2';                 % To identify figures and reports
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
inputs.model.par = [0.59, 0.594, 0.855, 0.645, 0.388, 0.45, 0.658, 0.148, 0.633, 0.637, 0.268, 0.203, 0.352];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.391, 0.556, 0.451, 0.891, 0.182];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = n', 'y2 = e', 'y3 = s+m', 'y4 = p');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.3910000000000000 0.5560000000000000 1.3420000000000001 0.1820000000000000
0.3774971009869568 0.5330116859081938 1.3222195792086584 0.1747627857266967
0.3645416633014960 0.5112105455433987 1.3032394168395791 0.1679264531142495
0.3521049336792028 0.4905163170263881 1.2850046655338296 0.1614604992192587
0.3401602336095659 0.4708557092473020 1.2674665707395745 0.1553373820123845
0.3286827595225796 0.4521616619157810 1.2505816255522355 0.1495321744864892
0.3176494069097969 0.4343726976324226 1.2343108600377812 0.1440222657028017
0.3070386150098897 0.4174323529518057 1.2186192410994492 0.1387871015499346
0.2968302292286252 0.4012886774697531 1.2034751636467571 0.1338079592314968
0.2870053789049702 0.3858937916642149 1.1888500175134182 0.1290677505016980
0.2775463684009816 0.3712034956303723 1.1747178174915658 0.1245508494882241
0.2684365797972992 0.3571769220277923 1.1610548861754688 0.1202429416143236
0.2596603857279172 0.3437762275355742 1.1478395811618405 0.1161308906827107
0.2512030710994738 0.3309663179351899 1.1350520596479017 0.1122026216405593
0.2430507626183371 0.3187146026349600 1.1226740746750192 0.1084470169238507
0.2351903651970345 0.3069907750300986 1.1106887982400968 0.1048538245924387
0.2276095044386755 0.2957666155893435 1.0990806672956950 0.1014135767312966
0.2202964745047471 0.2850158149773919 1.0878352493102237 0.0981175168131774
0.2132401907628608 0.2747138148797622 1.0769391245948507 0.0949575349033212
0.2064301466886522 0.2648376645010979 1.0663797830446695 0.0919261097429773
0.1998563745626053 0.2553658909688989 1.0561455333068372 0.0890162568808087
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
