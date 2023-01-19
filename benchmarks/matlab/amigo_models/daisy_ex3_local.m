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
inputs.exps.exp_y0{1}=[.9, -1.2, 0.8, -1.];        % Initial conditions
inputs.exps.t_f{1}=6;                       % Experiments duration
inputs.exps.n_obs{1}=2;                       % Number of observables
inputs.exps.obs_names{1}=char('y1','y2'); % Names of the observables
inputs.exps.obs{1}=char('y1=u0','y2=x1');
inputs.exps.t_con{1}=[0 6];                 % Input swithching times including:
inputs.exps.n_s{1}=50;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
  -1.0 1.0
  -0.8775510204081634 0.6838478755214819
  -0.7551020408163268 0.4548589875457535
  -0.63265306122449 0.2939552967454753
  -0.5102040816326533 0.1862685099119232
  -0.38775510204081653 0.12029042312170225
  -0.26530612244897983 0.08718426058875367
  -0.142857142857143 0.08022883173797922
  -0.02040816326530633 0.09437176935290122
  0.10204081632653049 0.1258719674087034
  0.22448979591836732 0.17201465510832742
  0.3469387755102039 0.23088537740238854
  0.469387755102041 0.3011915543919238
  0.5918367346938777 0.38212231652799455
  0.7142857142857145 0.4732390087388776
  0.8367346938775517 0.5743901692322241
  0.9591836734693885 0.6856459611281706
  1.081632653061225 0.8072480013725957
  1.2040816326530612 0.9395713262849236
  1.3265306122448979 1.0830958821328212
  1.4489795918367343 1.2383854590420975
  1.5714285714285714 1.4060724151063664
  1.6938775510204074 1.5868468847857902
  1.816326530612244 1.7814494443315287
  1.9387755102040805 1.9906664309001039
  2.0612244897959173 2.2153272899738465
  2.183673469387754 2.4563034679785365
  2.306122448979591 2.714508478811951
  2.4285714285714275 2.9908988616884273
  2.5510204081632644 3.2864758169372705
  2.673469387755101 3.602287360612938
  2.795918367346938 3.9394308807963117
  2.9183673469387745 4.299056011098294
  3.0408163265306105 4.68236776169366
  3.163265306122448 5.090629867755042
  3.2857142857142847 5.525168329297093
  3.4081632653061216 5.987375127341126
  3.5306122448979584 6.47871211012
  3.653061224489795 7.000715049289326
  3.775510204081632 7.554997869244719
  3.8979591836734677 8.143257058037435
  4.020408163265304 8.767276270423418
  4.1428571428571415 9.428931133574759
  4.265306122448978 10.13019427023488
  4.387755102040814 10.873140554151385
  4.510204081632652 11.659952611389738
  4.632653061224488 12.492926586024346
  4.7551020408163245 13.374478186117535
  4.877551020408162 14.307149026383103
  4.999999999999998 15.293613288882552

];
inputs.PEsol.id_global_theta=char('p1', 'p3');
inputs.PEsol.global_theta_max=2.*ones(1,2);
inputs.PEsol.global_theta_min=-1.*ones(1,2);
inputs.PEsol.id_global_theta_y0=char('x1', 'x2', 'u0');               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=2.*ones(1,3);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=-1.*ones(1,3);
inputs.PEsol.id_local_theta{1}=char('p4','p6','p7');                % [] 'all'|User selected| 'none' (default)
inputs.PEsol.local_theta_max{1}=2.*ones(1,3);              % Maximum allowed values for the paramters
inputs.PEsol.local_theta_min{1}=rand(1, 3, 1);              % Minimum allowed values for the parameters
inputs.PEsol.local_theta_guess{1}=[];            % [] Initial guess
inputs.PEsol.id_local_theta_y0{1}=char('x3');             % [] 'all'|User selected| 'none' (default)
inputs.PEsol.local_theta_y0_max{1}=2.*ones(1,1);           % Maximum allowed values for the initial conditions
inputs.PEsol.local_theta_y0_min{1}=-1.*ones(1,1);           % Minimum allowed values for the initial conditions
inputs.PEsol.local_theta_y0_guess{1}=rand(1, 1, 1);         % [] Initial guess


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
