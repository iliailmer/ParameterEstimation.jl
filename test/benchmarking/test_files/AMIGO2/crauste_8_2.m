addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='crauste_8model'; % Folder to keep results
inputs.pathd.short_name='crauste_8';                 % To identify figures and reports
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
inputs.model.par = [0.745, 0.663, 0.18, 0.836, 0.671, 0.899, 0.22, 0.795, 0.23, 0.592, 0.199, 0.778, 0.746];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.555, 0.426, 0.155, 0.658, 0.463];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = n', 'y2 = e', 'y3 = s+m', 'y4 = p');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.5550000000000000 0.4260000000000000 0.8130000000000001 0.4630000000000000
0.5319212341862344 0.4179003606495041 0.7942690702040981 0.4453609120710308
0.5099043972260947 0.4096225691188739 0.7762340994981501 0.4281289824395003
0.4888945631348754 0.4012026328419896 0.7588731049750057 0.4113161576998267
0.4688398951721718 0.3926738814822740 0.7421647968829785 0.3949323552304742
0.4496914717897017 0.3840670619059854 0.7260885556806275 0.3789855435551530
0.4314031212063168 0.3754104398039503 0.7106244099278797 0.3634818307053924
0.4139312643350815 0.3667299063705726 0.6957530149733264 0.3484255591132835
0.3972347657798423 0.3580490886016781 0.6814556324000232 0.3338194056079219
0.3812747926173609 0.3493894619165148 0.6677141101950799 0.3196644851497505
0.3660146806743400 0.3407704639520525 0.6545108636100452 0.3059604570177407
0.3514198080113653 0.3322096085075848 0.6418288566819913 0.2927056322528606
0.3374574753213682 0.3237225987474070 0.6296515843864777 0.2798970812644051
0.3240967929554527 0.3153234388843843 0.6179630553961561 0.2675307406093574
0.3113085742872805 0.3070245436820870 0.6067477754198991 0.2556015180685539
0.2990652351333334 0.2988368452137293 0.5959907310993084 0.2441033952515209
0.2873406989491949 0.2907698964132944 0.5856773744408375 0.2330295270739377
0.2761103075267558 0.2828319710422213 0.5757936077628967 0.2223723375578108
0.2653507369247000 0.2750301597734298 0.5663257691389344 0.2121236115055683
0.2550399183691300 0.2673704621691273 0.5572606183182396 0.2022745816973528
0.2451569638696437 0.2598578743922521 0.5485853231074175 0.1928160113480900
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
