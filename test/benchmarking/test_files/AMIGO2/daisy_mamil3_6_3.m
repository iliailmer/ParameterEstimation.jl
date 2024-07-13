addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='daisy_mamil3_6model'; % Folder to keep results
inputs.pathd.short_name='daisy_mamil3_6';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=3;                                  % Number of states:\\\
inputs.model.n_par=5;                                 % Number of model parameters
inputs.model.st_names=char('x1', 'x2', 'x3');    % Names of the states
inputs.model.par_names=char('a12', 'a13', 'a21', 'a31', 'a01');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char( 'dx1 = -(a21 + a31 + a01) * x1 + a12 * x2 + a13 * x3;',  'dx2 = a21 * x1 - a12 * x2;',  'dx3 = a31 * x1 - a13 * x3;');               % Equations describing system dynamics.
inputs.model.par = [0.352, 0.391, 0.556, 0.451, 0.891];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.182, 0.267, 0.229];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=2;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = x1', 'y2 = x2');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.1820000000000000 0.2670000000000000
0.1742749748082725 0.2672488236743937
0.1672416895015331 0.2672901545912315
0.1608314861734258 0.2671457094033790
0.1549825811523953 0.2668350305630272
0.1496393762377908 0.2663757039998363
0.1447518389510582 0.2657835549917670
0.1402749448865577 0.2650728244136493
0.1361681759394488 0.2642563273296087
0.1323950688117512 0.2633455956984082
0.1289228087585384 0.2623510067835641
0.1257218640410355 0.2612818987005845
0.1227656570079575 0.2601466743900563
0.1200302681346320 0.2589528951763284
0.1174941697175354 0.2577073649552355
0.1151379862527949 0.2564162059497467
0.1129442788247968 0.2550849268783992
0.1108973510990409 0.2537184842966951
0.1089830747545035 0.2523213377954568
0.1071887324076559 0.2508974996716056
0.1055028762752577 0.2494505796252266
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=3.0*ones(1,5);
inputs.PEsol.global_theta_min=0.0*ones(1,5);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=3.0*ones(1,3);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=0.0*ones(1,3);
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
inputs.nlpsol.eSS.log_var=1:(3+5); 
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
