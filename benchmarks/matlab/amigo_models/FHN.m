addpath(genpath('../src'))
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='FHNModel'; % Folder to keep results
inputs.pathd.short_name='FHN';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=2;                                  % Number of states
inputs.model.n_par=3;                                 % Number of model parameters
%inputs.model.n_stimulus=0;                            % Number of inputs, stimuli or control variables
inputs.model.st_names=char('x1','x2');    %x1=V, x2=R        % Names of the states
inputs.model.par_names=char('a','b','g');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char('dx1 = g * (x1 - x1^3/3 + x2);', 'dx2 = - 1/g * (x1 - a + b * x2);');                                 % Equations describing system dynamics.
inputs.model.par = [0.2 0.2 2];         % Nominal value for the parameters
% inputs.model.AMIGOsensrhs = 1;                       % Generate the sensitivity equations for exact
%                                                      % Jacobian computation
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[1 -1];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=1;                       % Number of observables
inputs.exps.obs_names{1}=char('Y'); % Names of the observables
inputs.exps.obs{1}=char('Y=x1');
inputs.exps.t_con{1}=[0 1];                 % Input swithching times including:
inputs.exps.n_s{1}=50;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[1.0
0.9865161738969255
0.9732626036107086
0.9602205595438085
0.9473722004384822
0.9347004680289507
0.9221889913641396
0.9098219994143794
0.8975842408051563
0.8854609097506593
0.8734375772154442
0.8615001267279976
0.8496346936983219
0.8378276081883363
0.8260653404882949
0.8143344485649917
0.8026215271698885
0.7909131589058271
0.7791958654565811
0.7674560596070894
0.7556799982697886
0.7438537339289144
0.7319630671786492
0.7199934979748167
0.7079301752128767
0.695757847083321
0.6834608068484177
0.6710228402258328
0.6584271675490089
0.6456563858105527
0.6326924069682979
0.6195163929532004
0.606108689372379
0.5924487520996189
0.5785150751481465
0.5642851085191068
0.5497351777452164
0.5348403938714783
0.5195745618912201
0.5039100836774605
0.48781785413809176
0.47126715568168104
0.454225544859093
0.4366587358017351
0.4185304800580947
0.39980244160012274
0.38043407010207164
0.3603824732191721
0.3396022892593786
0.31804556259642747
];
inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=2.*ones(1,3);
inputs.PEsol.global_theta_min=-1.*ones(1,3);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=2.*ones(1,2);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=-1.*ones(1,2);
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
