addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='daisy_mamil4_6model'; % Folder to keep results
inputs.pathd.short_name='daisy_mamil4_6';                 % To identify figures and reports
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
inputs.model.par = [0.757, 0.178, 0.77, 0.177, 0.881, 0.475, 0.881];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.584, 0.691, 0.131, 0.326];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=3;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = x1', 'y2 = x2', 'y3 = x3 + x4');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.5840000000000000 0.6909999999999999 0.4570000000000000
0.5161988411968281 0.7089600991054379 0.4860538446947205
0.4583691672362094 0.7240125336450635 0.5104163497368699
0.4090226643189218 0.7365867978064274 0.5307978409306379
0.3668937069773684 0.7470483330035340 0.5478020129655941
0.3309060213181989 0.7557081094965798 0.5619419037456101
0.3001443393695156 0.7628307737056630 0.5736534767065645
0.2738302972585524 0.7686415759570091 0.5833071692694486
0.2513019418288050 0.7733322612495022 0.5912177119528186
0.2319963054372512 0.7770660782922291 0.5976524770744066
0.2154345895615649 0.7799820388193615 0.6028385772035996
0.2012095666247325 0.7821985394247573 0.6069689005632413
0.1889748679243619 0.7838164413538375 0.6102072425513571
0.1784358752792720 0.7849216894011849 0.6126926687213960
0.1693419762846488 0.7855875389130684 0.6145432242965432
0.1614799790173072 0.7858764495633814 0.6158590880636752
0.1546685125988133 0.7858416957879348 0.6167252538429943
0.1487532660153333 0.7855287362931970 0.6172138102725533
0.1436029396917117 0.7849763787052635 0.6173858790556395
0.1391058031079796 0.7842177700251097 0.6172932628131412
0.1351667677231959 0.7832812389650663 0.6169798460258576
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=2.0*ones(1,7);
inputs.PEsol.global_theta_min=0.0*ones(1,7);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=2.0*ones(1,4);                % Maximum allowed values for the initial conditions
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
