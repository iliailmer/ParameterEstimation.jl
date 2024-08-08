addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='hiv_3model'; % Folder to keep results
inputs.pathd.short_name='hiv_3';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=5;                                  % Number of states:\\\
inputs.model.n_par=10;                                 % Number of model parameters
inputs.model.st_names=char('x', 'yy', 'vv', 'w', 'z');    % Names of the states
inputs.model.par_names=char('lm', 'd', 'beta', 'a', 'k', 'uu', 'c', 'q', 'b', 'h');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char( 'dx = lm - d * x - beta * x * vv;',  'dyy = beta * x * vv - a * yy;',  'dvv = k * yy - uu * vv;',  'dw = c * x * yy * w - c * q * yy * w - b * w;',  'dz = c * q * yy * w - h * z;');               % Equations describing system dynamics.
inputs.model.par = [0.637, 0.268, 0.203, 0.352, 0.391, 0.556, 0.451, 0.891, 0.182, 0.267];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.229, 0.622, 0.303, 0.473, 0.296];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = w', 'y2 = z', 'y3 = x', 'y4 = yy+vv');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.4730000000000000 0.2960000000000000 0.2290000000000000 0.9250000000000000
0.4645083636159486 0.2978458897456130 0.2568425517041134 0.9184808358328423
0.4564080932630584 0.2994735585665779 0.2842202728619913 0.9119414446371784
0.4486739869411567 0.3008961308522079 0.3111397377816203 0.9053898534135923
0.4412827645255092 0.3021258799590117 0.3376076472831905 0.8988335050090914
0.4342129020510644 0.3031742925881939 0.3636308059798169 0.8922792917079708
0.4274444818708946 0.3040521275381734 0.3892161012329683 0.8857335871925116
0.4209590570205362 0.3047694693828958 0.4143704836980774 0.8792022769293975
0.4147395283102162 0.3053357775676373 0.4391009493775670 0.8726907870381919
0.4087700328331140 0.3057599313620906 0.4634145231004161 0.8662041116985455
0.4030358427244309 0.3060502710643366 0.4873182433495578 0.8597468391526808
0.3975232731346059 0.3062146358085497 0.5108191483606829 0.8533231763594786
0.3922195984944616 0.3062603982926233 0.5339242634186181 0.8469369723556850
0.3871129762497028 0.3061944967098723 0.5566405892799255 0.8405917403791342
0.3821923773311812 0.3060234641400529 0.5789750916531028 0.8342906788078079
0.3774475227063148 0.3057534556290595 0.6009346916705672 0.8280366909672097
0.3728688254256083 0.3053902731641486 0.6225262572892715 0.8218324038575318
0.3684473376401680 0.3049393887309461 0.6437565955596189 0.8156801858505827
0.3641747021203257 0.3044059656204010 0.6646324457050307 0.8095821634050639
0.3600431078550051 0.3037948781370615 0.6851604729574254 0.8035402368469130
0.3560452493531545 0.3031107298460609 0.7053472630964327 0.7975560952603791
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=1.0*ones(1,10);
inputs.PEsol.global_theta_min=0.0*ones(1,10);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=1.0*ones(1,5);                % Maximum allowed values for the initial conditions
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
inputs.nlpsol.eSS.log_var=1:(5+10); 
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
