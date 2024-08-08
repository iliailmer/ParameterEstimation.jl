addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='crauste_6model'; % Folder to keep results
inputs.pathd.short_name='crauste_6';                 % To identify figures and reports
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
inputs.model.par = [0.278, 0.862, 0.458, 0.777, 0.66, 0.338, 0.751, 0.417, 0.805, 0.565, 0.805, 0.654, 0.68];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.501, 0.865, 0.615, 0.439, 0.585];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = n', 'y2 = e', 'y3 = s+m', 'y4 = p');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.5010000000000000 0.8650000000000000 1.0540000000000000 0.5850000000000000
0.4828217483244787 0.8370331110745929 1.0306808219874504 0.5610221732002226
0.4657384178286545 0.8098435426053429 1.0085796131749365 0.5385112017257521
0.4496543971528523 0.7834621101351719 0.9875364367926516 0.5173497680109930
0.4344847334618718 0.7579056697010640 0.9674186164662896 0.4974316597041598
0.4201537212349509 0.7331801644968544 0.9481156593183564 0.4786605931494619
0.4065937042166720 0.7092830729135189 0.9295352089839608 0.4609491727344572
0.3937440546818666 0.6862053725870602 0.9115998000179076 0.4442179686933785
0.3815503007777465 0.6639331135070091 0.8942442404940472 0.4283946987055774
0.3699633779913837 0.6424486755689737 0.8774134903922264 0.4134135007515562
0.3589389850399481 0.6217317715514182 0.8610609337488284 0.3992142863903652
0.3484370279103093 0.6017602448035042 0.8451469653680925 0.3857421650213667
0.3384211385587396 0.5825107014543022 0.8296378301901919 0.3729469308742988
0.3288582570476364 0.5639590092930289 0.8145046666185694 0.3607826054804251
0.3197182677501104 0.5460806892809917 0.7997227152705598 0.3492070292511245
0.3109736817749883 0.5288512206604182 0.7852706624856798 0.3381814965514712
0.3025993590180833 0.5122462765943891 0.7711300940607377 0.3276704293224509
0.2945722642817575 0.4962419040191398 0.7572850394920052 0.3176410848903766
0.2868712527642426 0.4808146587644379 0.7437215907996124 0.3080632941164724
0.2794768809348125 0.4659417048781373 0.7304275830167207 0.2989092264916741
0.2723712394076036 0.4516008853784238 0.7173923258236700 0.2901531791799983
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=3.0*ones(1,13);
inputs.PEsol.global_theta_min=0.0*ones(1,13);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=3.0*ones(1,5);                % Maximum allowed values for the initial conditions
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
