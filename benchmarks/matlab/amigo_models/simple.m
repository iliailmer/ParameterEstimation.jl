addpath(genpath('../src'))
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='SimpleModel'; % Folder to keep results
inputs.pathd.short_name='Simple';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=2;                                  % Number of states
inputs.model.n_par=2;                                 % Number of model parameters
%inputs.model.n_stimulus=0;                            % Number of inputs, stimuli or control variables
inputs.model.st_names=char('x1','x2');    %x1=V, x2=R        % Names of the states
inputs.model.par_names=char('a','b');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char('dx1 = -a * x2;','dx2 = 1 / b * (x1);');                                 % Equations describing system dynamics.
inputs.model.par = [0.2 0.2];         % Nominal value for the parameters
% inputs.model.AMIGOsensrhs = 1;                       % Generate the sensitivity equations for exact
%                                                      % Jacobian computation
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[1 1];        % Initial conditions
inputs.exps.t_f{1}=2.0 * 3.14159265 * sqrt(1.3 / 9.8);                       % Experiments duration
inputs.exps.n_obs{1}=1;                       % Number of observables
inputs.exps.obs_names{1}=char('Y1', 'Y2'); % Names of the observable
inputs.exps.obs{1}=char('Y1=x2', 'Y2=x1');
inputs.exps.t_con{1}=[0 2.0 * 3.14159265 * sqrt(1.3 / 9.8)];                 % Input swithching times including:
inputs.exps.n_s{1}=20;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[1.0 1.0
1.0367869631032143 -0.2131370089757947
0.9612219713463254 -1.4031773158605787
0.7814936640945697 -2.441161588047151
0.5170783922159042 -3.2146081236261184
0.19662965324200823 -3.6397019892197804
-0.14512695968147638 -3.670377668485929
-0.4711568146481816 -3.3033109755914225
-0.7461295179881701 -2.578279282359504
-0.9402475106698658 -1.5738510227428189
-1.0324750961248383 -0.3988715840136487
-1.0128179844141494 0.8193317799828911
-0.8834063286038646 1.9487478323598202
-0.658263889624748 2.866986819178007
-0.36178834419170564 3.4745432982439945
-0.02610741794112031 3.70557909785514
0.31240265214351715 3.5350579042280135
0.6170590474382538 2.9814583346094343
0.854847520288421 2.1047714923332395
0.9999999999859759 0.9999999999744436
];
inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=10.*ones(1,2);
inputs.PEsol.global_theta_min=0.0001.*ones(1,2);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=1.5.*ones(1,2);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=0.*ones(1,2);
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
