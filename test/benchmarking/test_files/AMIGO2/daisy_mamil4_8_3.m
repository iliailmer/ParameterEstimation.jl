addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='daisy_mamil4_8model'; % Folder to keep results
inputs.pathd.short_name='daisy_mamil4_8';                 % To identify figures and reports
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
inputs.model.par = [0.561, 0.843, 0.355, 0.634, 0.205, 0.673, 0.332];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.247, 0.569, 0.116, 0.763];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=3;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = x1', 'y2 = x2', 'y3 = x3 + x4');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.2470000000000000 0.5690000000000000 0.8790000000000000
0.2734599760014158 0.5481310430403750 0.8661008259904810
0.2964008098333949 0.5283708887847428 0.8549200868850205
0.3162040896635891 0.5096404144902540 0.8452485791815995
0.3332115771876493 0.4918673882207092 0.8368999002080272
0.3477293215934701 0.4749858058621626 0.8297080555315661
0.3600313498952317 0.4589352946150209 0.8235253141743994
0.3703629772164363 0.4436605761910641 0.8182202860879771
0.3789437761192086 0.4291109836377194 0.8136761989687695
0.3859702400568005 0.4152400263344496 0.8097893538565784
0.3916181724184978 0.4020049982648651 0.8064677410730191
0.3960448294005493 0.3893666251693469 0.8036297999558787
0.3993908420323159 0.3772887466329582 0.8012033075487683
0.4017819400826466 0.3657380295669571 0.7991243829330653
0.4033304982337602 0.3546837099045669 0.7973365952603005
0.4041369228134268 0.3440973596566714 0.7957901647727837
0.4042908954950873 0.3339526767648115 0.7944412472035156
0.4038724886872760 0.3242252954507376 0.7932512929365169
0.4029531658203978 0.3148926149965523 0.7921864731961357
0.4015966783795869 0.3059336451004673 0.7912171663309465
0.3998598703139260 0.2973288661423555 0.7903174979722731
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
