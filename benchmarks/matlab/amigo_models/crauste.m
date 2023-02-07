addpath(genpath('../src'))
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='CrausteModel'; % Folder to keep results
inputs.pathd.short_name='Crauste';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=5;                                  % Number of states
inputs.model.n_par=13;                                 % Number of model parameters
%inputs.model.n_stimulus=0;                            % Number of inputs, stimuli or control variables
inputs.model.st_names=char('n', 'e', 's', 'm', 'p');    %x1=V, x2=R        % Names of the states
inputs.model.par_names=char('muN', 'muEE', 'muLE', 'muLL', 'muM', 'muP', 'muPE', 'muPL', 'deltaNE', 'deltaEL', 'deltaLM', 'rhoE', 'rhoP');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
% Equations describing system dynamics.
inputs.model.eqns=char( 'dn = -1 * n * muN - n * p * deltaNE;', 'de = n * p * deltaNE - e * e * muEE - e * deltaEL + e * e * rhoE;', 'ds = s * deltaEL - s * deltaLM - s * s * muLL - e * s * muLE;', 'dm = s * deltaLM - muM * m;', 'dp = p * p * rhoP - p * muP - e * p * muPE - s * p * muPL;');
inputs.model.par = [1 1.3 1.1 1.2 1.1 1 0.5 1.0 1.0 1.0 1.0 0.9 1.2];         % Nominal value for the parameters
% inputs.model.AMIGOsensrhs = 1;                       % Generate the sensitivity equations for exact
%                                                      % Jacobian computation
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[1.0 1.0 1.0 1.0 1.0];        % Initial conditions
inputs.exps.t_f{1}=5;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
% Names of the observables
inputs.exps.obs_names{1}=char('Y1', 'Y2', 'Y3', 'Y4');
inputs.exps.obs{1}=char('Y1=e', 'Y2=n', 'Y3=s+m', 'Y4=p');
inputs.exps.t_con{1}=[0 1];                 % Input swithching times including:
inputs.exps.n_s{1}=20;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
1.0 1.0 2.0 1.0
0.9748203865016476 0.901654713757981 1.8817413878377107 0.9350058134120561
0.943347445426423 0.8156323245950008 1.777306727235583 0.8761075966493853
0.9079222171927342 0.7400045430471617 1.6839562031649669 0.8223637276867406
0.8702449097965126 0.6732108218409861 1.5996505976905997 0.773035031165633
0.8315374183857062 0.6139748727925719 1.5228475501525338 0.7275322306009067
0.7926678354998998 0.561243198372053 1.4523632196845826 0.6853792958497149
0.7542442024752796 0.5141390488711327 1.3872763059618005 0.6461871652736274
0.7166843142465519 0.4719274120521719 1.3268602697313623 0.6096344732979303
0.6802673260817139 0.4339880355696268 1.2705347510812985 0.5754531481934801
0.6451717428101934 0.3997943895931915 1.2178302971067458 0.5434174800175874
0.6115033121731885 0.3688970858299796 1.1683624628765426 0.5133357145963535
0.579315478604056 0.3409106807336174 1.1218125990678502 0.4850435186966817
0.5486243706785097 0.31550307902818403 1.0779134667535657 0.4583988530605419
0.5194197777836486 0.2923869562314639 1.0364383730341487 0.43327791861873644
0.49167318493966367 0.271312763683446 0.9971928986212648 0.409571929907168
0.4653436483809577 0.25206298522799536 0.960008550346444 0.3871845324184962
0.440382084946922 0.2344473918875798 0.9247378545211429 0.3660297253711122
0.4167343950949036 0.21829909798985445 0.8912505365733951 0.3460301838872221
0.39434372673007195 0.2034712662961829 0.8594305260723616 0.3271158990871384
];
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=1.1 * ones(1,5);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=-1.1 * ones(1,5);
inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=2.*ones(1,13);
inputs.PEsol.global_theta_min=-2.*ones(1,13);
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
