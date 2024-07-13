addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='crauste_9model'; % Folder to keep results
inputs.pathd.short_name='crauste_9';                 % To identify figures and reports
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
inputs.model.par = [0.678, 0.793, 0.88, 0.785, 0.109, 0.388, 0.684, 0.237, 0.517, 0.143, 0.26, 0.115, 0.735];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.279, 0.376, 0.842, 0.664, 0.125];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = n', 'y2 = e', 'y3 = s+m', 'y4 = p');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.2790000000000000 0.3760000000000000 1.5060000000000000 0.1250000000000000
0.2688464160820143 0.3689674395353237 1.4684333470204547 0.1204186739803621
0.2590922334161361 0.3621114303432829 1.4341317589131486 0.1160723769357315
0.2497192934675315 0.3554278044596481 1.4026836041911119 0.1119431485627929
0.2407104884952256 0.3489122936276555 1.3737415363696446 0.1080150161318365
0.2320496773441354 0.3425605739965818 1.3470105872003786 0.1042737013431940
0.2237216107022548 0.3363683020179464 1.3222387958402291 0.1007063804910865
0.2157118644076275 0.3303311432118669 1.2992097673265119 0.0973014866038163
0.2080067796531320 0.3244447951355259 1.2777367138214155 0.0940485449506242
0.2005934091424331 0.3187050056215257 1.2576576461945597 0.0909380353041357
0.1934594684129006 0.3131075871470049 1.2388314658961108 0.0879612758334181
0.1865932916702846 0.3076484280313091 1.2211347671775448 0.0851103246177683
0.1799837915840796 0.3023235010303564 1.2044592040875981 0.0823778956179657
0.1736204225766610 0.2971288697922550 1.1887093097090151 0.0797572865891473
0.1674931472083003 0.2920606935553809 1.1738006799620866 0.0772423169198496
0.1615924053167099 0.2871152304027825 1.1596584531436194 0.0748272737707476
0.1559090856171824 0.2822888393320082 1.1462160308043994 0.0725068651926142
0.1504344995082259 0.2775779813548297 1.1334139966637378 0.0702761791437157
0.1451603568610602 0.2729792198047089 1.1211991989021948 0.0681306475194598
0.1400787435987464 0.2684892199997770 1.1095239679037463 0.0660660144603140
0.1351821008947789 0.2641047483842524 1.0983454468362002 0.0640783083283742
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=3.0*ones(1,13);
inputs.PEsol.global_theta_min=0.0*ones(1,13);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=3.0*ones(1,5);                % Maximum allowed values for the initial conditions
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
