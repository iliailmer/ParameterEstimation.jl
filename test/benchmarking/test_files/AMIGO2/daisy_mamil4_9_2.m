addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='daisy_mamil4_9model'; % Folder to keep results
inputs.pathd.short_name='daisy_mamil4_9';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=4;                                  % Number of states:\\\
inputs.model.n_par=7;                                 % Number of model parameters
inputs.model.st_names=char('x1', 'x2', 'x3', 'x4');    % Names of the states
inputs.model.par_names=char('k01', 'k12', 'k13', 'k14', 'k21', 'k31', 'k41');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char( 'dx1 = -k01 * x1 + k12 * x2 + k13 * x3 + k14 * x4 - k21 * x1 - k31 * x1 - k41 * x1;',  'dx2 = -k12 * x2 + k21 * x1;',  'dx3 = -k13 * x3 + k31 * x1;',  'dx4 = -k14 * x4 + k41 * x1;');               % Equations describing system dynamics.
inputs.model.par = [0.104, 0.642, 0.316, 0.688, 0.87, 0.299, 0.561];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.574, 0.558, 0.278, 0.862];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=3;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = x1', 'y2 = x2', 'y3 = x3 + x4');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.5740000000000000 0.5580000000000001 1.1399999999999999
0.5732681978252925 0.5649297554857344 1.1308191045496845
0.5724402083309379 0.5716071495002856 1.1219908195328476
0.5715294480005325 0.5780363312168398 1.1134980438547528
0.5705477487316511 0.5842218582969054 1.1053247869562326
0.5695055397093809 0.5901686201435092 1.0974560714443888
0.5684120085352597 0.5958817709092308 1.0898778458223406
0.5672752439736710 0.6013666711149435 1.0825769062016866
0.5661023624097523 0.6066288368661137 1.0755408260068211
0.5648996198740046 0.6116738957701180 1.0687578927922363
0.5636725112776586 0.6165075487611107 1.0622170513929057
0.5624258583157050 0.6211355371303002 1.0559078527156305
0.5611638873286436 0.6255636141403950 1.0498204075570370
0.5598902982670878 0.6297975206746237 1.0439453449028653
0.5586083257726365 0.6338429644343702 1.0382737742246584
0.5573207932739782 0.6377056022553484 1.0327972513436023
0.5560301608935800 0.6413910251624840 1.0275077474801737
0.5547385678711101 0.6449047458272892 1.0223976211500081
0.5534478701280056 0.6482521881310138 1.0174595926051189
0.5521596735275928 0.6514386785710746 1.0126867205523715
0.5508753633213546 0.6544694392790940 1.0080723809113663
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=2.0*ones(1,7);
inputs.PEsol.global_theta_min=0.0*ones(1,7);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=2.0*ones(1,4);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=0.0*ones(1,4);
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
inputs.nlpsol.eSS.log_var=1:(4+7); 
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
