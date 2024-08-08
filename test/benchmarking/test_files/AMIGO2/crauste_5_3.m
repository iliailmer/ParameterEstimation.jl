addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='crauste_5model'; % Folder to keep results
inputs.pathd.short_name='crauste_5';                 % To identify figures and reports
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
inputs.model.par = [0.355, 0.634, 0.205, 0.673, 0.332, 0.247, 0.569, 0.116, 0.763, 0.104, 0.642, 0.316, 0.688];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.87, 0.299, 0.561, 0.574, 0.558];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = n', 'y2 = e', 'y3 = s+m', 'y4 = p');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.8700000000000000 0.2990000000000000 1.1350000000000000 0.5580000000000001
0.8367364851338506 0.3152300016204651 1.1164355795379897 0.5551508060476220
0.8048354303980548 0.3304036601439784 1.0985359205656875 0.5520958213858890
0.7742437549273141 0.3445508582995369 1.0812361528128744 0.5488461613652542
0.7449098411638229 0.3577033469626814 1.0644800434688495 0.5454130928674710
0.7167835572358224 0.3698943808064739 1.0482187001705712 0.5418079337696080
0.6898162736421802 0.3811583761230886 1.0324094903955705 0.5380419650603893
0.6639608741089289 0.3915305930045497 1.0170151374576928 0.5341263547633416
0.6391717606517086 0.4010468432777046 1.0020029613068284 0.5300720928170968
0.6154048530098120 0.4097432248993022 0.9873442385909926 0.5258899360626541
0.5926175827179356 0.4176558829293890 0.9730136613569260 0.5215903624928794
0.5707688821557922 0.4248207967134997 0.9589888776585164 0.5171835339325225
0.5498191689670312 0.4312735925081365 0.9452501004359564 0.5126792663375368
0.5297303262716200 0.4370493804748406 0.9317797735051150 0.5080870069303527
0.5104656791129922 0.4421826147348494 0.9185622854845690 0.5034158174217055
0.4919899675857975 0.4467069750112700 0.9055837240931524 0.4986743626088158
0.4742693170846682 0.4506552682794266 0.8928316645526927 0.4938709036828003
0.4572712061011688 0.4540593487894849 0.8802949868902292 0.4890132956239760
0.4409644319769675 0.4569500548103554 0.8679637178001431 0.4841089881110603
0.4253190749979257 0.4593571604623980 0.8558288934372240 0.4791650294181281
0.4103064611875256 0.4613093410515528 0.8438824400966145 0.4741880728207105
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
