addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='crauste_4model'; % Folder to keep results
inputs.pathd.short_name='crauste_4';                 % To identify figures and reports
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
inputs.model.par = [0.881, 0.584, 0.691, 0.131, 0.326, 0.196, 0.337, 0.195, 0.354, 0.431, 0.151, 0.654, 0.553];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.312, 0.519, 0.175, 0.561, 0.843];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = n', 'y2 = e', 'y3 = s+m', 'y4 = p');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.3120000000000000 0.5190000000000000 0.7360000000000000 0.8430000000000000
0.2941262353473193 0.5188086661915069 0.7273516679354821 0.8456166698884126
0.2772633771728470 0.5184233756335711 0.7188294580758141 0.8483152237692477
0.2613546177416955 0.5178634723969073 0.7104325299842267 0.8511009596565763
0.2463463211253271 0.5171471327970272 0.7021599472115644 0.8539791075694231
0.2321878489044934 0.5162914309729127 0.6940106858874944 0.8569548481527162
0.2188313950187223 0.5153124027308480 0.6859836427683667 0.8600333311404929
0.2062318293442074 0.5142251074700145 0.6780776427506692 0.8632196936957077
0.1943465495913449 0.5130436880586157 0.6702914458632694 0.8665190786689045
0.1831353411263809 0.5117814285729761 0.6626237537547455 0.8699366528249451
0.1725602443331800 0.5104508098501562 0.6550732156945014 0.8734776250937036
0.1625854291444174 0.5090635628364836 0.6476384341079668 0.8771472649068252
0.1531770763850286 0.5076307197412716 0.6403179696672275 0.8809509206886771
0.1443032655841519 0.5061626630276242 0.6331103459589937 0.8848940385755978
0.1359338689253894 0.5046691722910270 0.6260140537519758 0.8889821814435808
0.1280404510189724 0.5031594690918934 0.6190275548855732 0.8932210483307336
0.1205961741929109 0.5016422598210192 0.6121492858013642 0.8976164943473767
0.1135757090134741 0.5001257766873783 0.6053776607382778 0.9021745511736500
0.1069551497585099 0.4986178169261834 0.5987110746115470 0.9069014482520340
0.1007119345798680 0.4971257803320988 0.5921479055946618 0.9118036347904740
0.0948247701036751 0.4956567052281818 0.5856865174225676 0.9168878027009364
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
