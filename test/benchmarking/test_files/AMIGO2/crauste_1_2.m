addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='crauste_1model'; % Folder to keep results
inputs.pathd.short_name='crauste_1';                 % To identify figures and reports
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
inputs.model.par = [0.723, 0.796, 0.883, 0.739, 0.469, 0.724, 0.195, 0.612, 0.215, 0.856, 0.517, 0.432, 0.312];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.719, 0.465, 0.555, 0.115, 0.594];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = n', 'y2 = e', 'y3 = s+m', 'y4 = p');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.7190000000000000 0.4650000000000000 0.6700000000000000 0.5940000000000000
0.6891630869853785 0.4473001347914089 0.6684820282663858 0.5659435456145435
0.6607586121162276 0.4301105355163440 0.6673335428525521 0.5392782452529776
0.6337019647336311 0.4134407757141161 0.6665121204077256 0.5139261473327328
0.6079147933522117 0.3972957282447690 0.6659804656207795 0.4898148057230319
0.5833244413038073 0.3816763813951085 0.6657057594354236 0.4668767580574586
0.5598634423512669 0.3665805290932351 0.6656590947586090 0.4450490669541278
0.5374690689179602 0.3520033542647910 0.6658149872518803 0.4242729153382373
0.5160829266042802 0.3379379215304371 0.6661509506405542 0.4044932484149577
0.4956505895284845 0.3243755930094418 0.6666471275377742 0.3856584559704136
0.4761212717625039 0.3113063789132614 0.6672859681020809 0.3677200896242419
0.4574475307548486 0.2987192328402057 0.6680519499659446 0.3506326104495089
0.4395849991641689 0.2866023001730943 0.6689313338229583 0.3343531630438646
0.4224921419805031 0.2749431267018193 0.6699119498670975 0.3188413726984533
0.4061300362002437 0.2637288335074751 0.6709830109622130 0.3040591627867700
0.3904621706556057 0.2529462632254000 0.6721349490025293 0.2899705898987294
0.3754542638881441 0.2425821020262885 0.6733592714210709 0.2765416945873932
0.3610740982056311 0.2326229809960159 0.6746484352259471 0.2637403658868588
0.3472913682781973 0.2230555600375698 0.6759957363055127 0.2515362180078634
0.3340775428179760 0.2138665969468723 0.6773952120520630 0.2399004778294750
0.3214057380508464 0.2050430039149594 0.6788415556179935 0.2288058819866124
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
