addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='crauste_7model'; % Folder to keep results
inputs.pathd.short_name='crauste_7';                 % To identify figures and reports
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
inputs.model.par = [0.115, 0.341, 0.628, 0.332, 0.594, 0.443, 0.208, 0.339, 0.556, 0.573, 0.559, 0.623, 0.622];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.445, 0.817, 0.394, 0.449, 0.814];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = n', 'y2 = e', 'y3 = s+m', 'y4 = p');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.4450000000000000 0.8169999999999999 0.8430000000000000 0.8139999999999999
0.4326074004883782 0.8127391866912219 0.8283902968723780 0.8042614938009350
0.4206733553382767 0.8079853861070453 0.8140251635994996 0.7945973001161656
0.4091780242419233 0.8027747122910718 0.7999049577203210 0.7850067087124086
0.3981025383187725 0.7971414392274664 0.7860294837806192 0.7754890025750710
0.3874289476700515 0.7911180627206000 0.7723980650044528 0.7660434654928384
0.3771401719363901 0.7847353644329854 0.7590096068071515 0.7566693886112744
0.3672199536878708 0.7780224773096683 0.7458626530099156 0.7473660760615443
0.3576528144846342 0.7710069517202059 0.7329554355362345 0.7381328497597406
0.3484240134532260 0.7637148217458720 0.7202859182980093 0.7289690534629303
0.3395195082324892 0.7561706711245545 0.7078518359098220 0.7198740561593656
0.3309259181501848 0.7483976984434267 0.6956507278066194 0.7108472548625627
0.3226304894990981 0.7404177812384140 0.6836799682816149 0.7018880768719608
0.3146210627885889 0.7322515387208843 0.6719367929077066 0.6929959815565830
0.3068860418545960 0.7239183929060379 0.6604183217567718 0.6841704617124519
0.2994143647177908 0.7154366279650836 0.6491215797868056 0.6754110445394348
0.2921954760859607 0.7068234476649083 0.6380435147266370 0.6667172922786383
0.2852193014028065 0.6980950307950635 0.6271810127516080 0.6580888025473803
0.2784762223511402 0.6892665845130820 0.6165309122108729 0.6495252084051171
0.2719570537239688 0.6803523955658891 0.6060900156375584 0.6410261781804182
0.2656530215821637 0.6713658793678420 0.5958551002466739 0.6325914150861581
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=2.0*ones(1,13);
inputs.PEsol.global_theta_min=0.0*ones(1,13);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=2.0*ones(1,5);                % Maximum allowed values for the initial conditions
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
