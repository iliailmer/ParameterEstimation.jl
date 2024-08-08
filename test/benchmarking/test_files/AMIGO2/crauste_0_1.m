addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='crauste_0model'; % Folder to keep results
inputs.pathd.short_name='crauste_0';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=5;                                  % Number of states:\\\
inputs.model.n_par=13;                                 % Number of model parameters
inputs.model.st_names=char('n', 'e', 's', 'm', 'p');    % Names of the states
inputs.model.par_names=char('muN', 'muEE', 'muLE', 'muLL', 'muM', 'muP', 'muPE', 'muPL', 'deltaNE', 'deltaEL', 'deltaLM', 'rhoE', 'rhoP');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char( 'dn = -1 * n * muN - n * p * deltaNE;',  'de = n * p * deltaNE - e * e * muEE - e * deltaEL + e * p * rhoE;',  'ds = s * deltaEL - s * deltaLM - s * s * muLL - e * s * muLE;',  'dm = s * deltaLM - muM * m;',  'dp = p * p * rhoP - p * muP - e * p * muPE - s * p * muPL;');               % Equations describing system dynamics.
inputs.model.par = [0.539, 0.672, 0.582, 0.536, 0.439, 0.617, 0.45, 0.813, 0.871, 0.407, 0.733, 0.523, 0.554];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.84, 0.157, 0.17, 0.116, 0.766];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = n', 'y2 = e', 'y3 = s+m', 'y4 = p');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.8400000000000000 0.1570000000000000 0.2860000000000000 0.7660000000000000
0.7911043265948752 0.1828602063795358 0.2852440397323129 0.7504589906423921
0.7455640203672603 0.2062429454270441 0.2842671225606205 0.7346359189762957
0.7031331065436138 0.2272833287676913 0.2830907845300768 0.7185955642618805
0.6635827800238200 0.2461170537186840 0.2817349638706559 0.7023982726299909
0.6267005219041951 0.2628789752907657 0.2802180428539484 0.6860999193221017
0.5922891969905865 0.2777018997904639 0.2785569071491489 0.6697519353753321
0.5601661474526684 0.2907155867004264 0.2767670178624668 0.6534013876015446
0.5301622950404167 0.3020459426364223 0.2748624921435397 0.6370911015391310
0.5021212618127441 0.3118143894430454 0.2728561889190836 0.6208598180251251
0.4758985171329837 0.3201373876975292 0.2707597969452596 0.6047423750748173
0.4513605567715172 0.3271260968502869 0.2685839229416269 0.5887699078105221
0.4283841183115333 0.3328861537640510 0.2663381780740698 0.5729700602102753
0.4068554356666754 0.3375175523639809 0.2640312614893429 0.5573672034206539
0.3866695343656166 0.3411146083432906 0.2616710399718756 0.5419826562769893
0.3677295683168197 0.3437659942719209 0.2592646230982371 0.5268349044870719
0.3499461980105130 0.3455548319403413 0.2568184335122185 0.5119398156566587
0.3332370095195960 0.3465588302683555 0.2543382721404178 0.4973108479666928
0.3175259732083933 0.3468504585675268 0.2518293783208395 0.4829592508546304
0.3027429407167999 0.3464971463316698 0.2492964849330011 0.4688942565146510
0.2888231785488668 0.3455615020168123 0.2467438687023843 0.4551232614168310
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=1.0*ones(1,13);
inputs.PEsol.global_theta_min=0.0*ones(1,13);
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
inputs.nlpsol.eSS.log_var=1:(5+13); 
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
