addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='biohydrogenation_6model'; % Folder to keep results
inputs.pathd.short_name='biohydrogenation_6';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=4;                                  % Number of states:\\\
inputs.model.n_par=6;                                 % Number of model parameters
inputs.model.st_names=char('x4', 'x5', 'x6', 'x7');    % Names of the states
inputs.model.par_names=char('k5', 'k6', 'k7', 'k8', 'k9', 'k10');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char( 'dx4 = - k5 * x4 / (k6 + x4);',  'dx5 = k5 * x4 / (k6 + x4) - k7 * x5/(k8 + x5 + x6);',  'dx6 = k7 * x5 / (k8 + x5 + x6) - k9 * x6 * (k10 - x6) / k10;',  'dx7 = k9 * x6 * (k10 - x6) / k10;');               % Equations describing system dynamics.
inputs.model.par = [0.227, 0.188, 0.625, 0.211, 0.257, 0.395];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.757, 0.178, 0.77, 0.177];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=2;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = x4', 'y2 = x5');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.7570000000000000 0.1780000000000000
0.7479189130374432 0.1822629193475097
0.7388599215421852 0.1864694463277067
0.7298234041965543 0.1906239586446644
0.7208097483125036 0.1947306880629483
0.7118193500570725 0.1987937362394648
0.7028526146835158 0.2028170898478370
0.6939099567681575 0.2068046351011285
0.6849918004530260 0.2107601717730417
0.6760985796943200 0.2146874268146045
0.6672307385167282 0.2185900676591407
0.6583887312736288 0.2224717153073259
0.6495730229131712 0.2263359572829252
0.6407840892502248 0.2301863605499277
0.6320224172441624 0.2340264844832340
0.6232885052824257 0.2378598939874635
0.6145828634697925 0.2416901728620959
0.6059060139232411 0.2455209375161184
0.5972584910722769 0.2493558511417572
0.5886408419645544 0.2531986384649854
0.5800536265765862 0.2570531012000476
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta=char('k5', 'k6', 'k7');
inputs.PEsol.global_theta_max=2.0*ones(1,3);
inputs.PEsol.global_theta_min=0.0*ones(1,3);
inputs.PEsol.id_global_theta_y0=char('x4', 'x5', 'x7');               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=2.0*ones(1,3);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=0.0*ones(1,3)

inputs.PEsol.id_local_theta{1}=char('k8', 'k9', 'k10');                % [] 'all'|User selected| 'none' (default)
inputs.PEsol.local_theta_max{1}=2.0*ones(1,3);              % Maximum allowed values for the paramters
inputs.PEsol.local_theta_min{1}=0.0*ones(1,3)              % Minimum allowed values for the parameters
inputs.PEsol.local_theta_guess{1}=[];            % [] Initial guess
inputs.PEsol.id_local_theta_y0{1}=char('x6');             % [] 'all'|User selected| 'none' (default)
inputs.PEsol.local_theta_y0_max{1}=2.0*ones(1,1);           % Maximum allowed values for the initial conditions
inputs.PEsol.local_theta_y0_min{1}=0.0*ones(1,1);           % Minimum allowed values for the initial conditions
inputs.PEsol.local_theta_y0_guess{1}=[];         % [] Initial guess

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
inputs.nlpsol.eSS.log_var=1:(4+6); 
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

%val names
inputs.model.par_names
inputs.model.st_names
%true vals:
inputs.model.par
inputs.exps.exp_y0{1}



