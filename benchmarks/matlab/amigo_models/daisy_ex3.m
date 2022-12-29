addpath(genpath('../src'))
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
inputs.exps.t_con{1}=[0 10];                 % Input swithching times including:
inputs.exps.n_s{1}=20;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
  2.0 2.0
2.398397183165376 4.311297896215865
3.437496304116641 7.756032933976405
5.432983428532577 12.964346593394282
8.855355771213679 20.89723022202182
14.420341755103157 33.0140832468806
23.22217324005539 51.531459044211715
36.93369584591418 79.81933096444803
58.10949948167877 123.00588449956251
90.64684806270999 188.8993648362735
140.48757650341386 289.39259564654105
216.68744401444195 442.60268197507503
333.04550022023477 676.1307380511039
510.58673341437856 1032.0280414343965
781.344754503689 1574.3620704312536
1194.1251429194422 2400.7436780239086
1823.2864201976292 3659.889359963373
2782.1185344259616 5578.378418935048
4243.22591449712 8501.419177957636
6469.582410768947 12954.958906505457
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
% %
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
% %
% inputs.exps.u_interp{1}='sustained';          % Stimuli definition for experiment 1
                                              % Initial and final time
%inputs.exps.u{1}=1;                           % Values of the inputs for exp 1
AMIGO_Prep(inputs);
AMIGO_PE(inputs);
