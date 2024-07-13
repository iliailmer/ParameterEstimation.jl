addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='hiv_0model'; % Folder to keep results
inputs.pathd.short_name='hiv_0';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=5;                                  % Number of states:\\\
inputs.model.n_par=10;                                 % Number of model parameters
inputs.model.st_names=char('x', 'yy', 'vv', 'w', 'z');    % Names of the states
inputs.model.par_names=char('lm', 'd', 'beta', 'a', 'k', 'uu', 'c', 'q', 'b', 'h');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char( 'dx = lm - d * x - beta * x * vv;',  'dyy = beta * x * vv - a * yy;',  'dvv = k * yy - uu * vv;',  'dw = c * x * yy * w - c * q * yy * w - b * w;',  'dz = c * q * yy * w - h * z;');               % Equations describing system dynamics.
inputs.model.par = [0.539, 0.672, 0.582, 0.536, 0.439, 0.617, 0.45, 0.813, 0.871, 0.407];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.733, 0.523, 0.554, 0.84, 0.157];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = w', 'y2 = z', 'y3 = x', 'y4 = yy+vv');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.8400000000000000 0.1570000000000000 0.7330000000000000 1.0770000000000000
0.8034042124831940 0.1616004828939251 0.7237943588094109 1.0691470764604003
0.7683261449353638 0.1657355479399720 0.7151466703585995 1.0612170889443830
0.7347120244300112 0.1694291351674156 0.7070250232618626 1.0532211553668609
0.7025089825526932 0.1727045734740796 0.6993995100901538 1.0451695635468199
0.6716651530646888 0.1755845251696008 0.6922420852749919 1.0370718368874261
0.6421297556123075 0.1780909433973471 0.6855264345907458 1.0289367939664031
0.6138531661878021 0.1802450407105925 0.6792278551216586 1.0207726026922432
0.5867869751657481 0.1820672672709434 0.6733231447385032 1.0125868296029985
0.5608840338097465 0.1835772973104546 0.6677905002137920 1.0043864848159181
0.5360984901790620 0.1847940226550048 0.6626094231961304 0.9961780630764514
0.5123858153719059 0.1857355522468269 0.6577606333451262 0.9879675813030686
0.4897028210284430 0.1864192167300112 0.6532259879996730 0.9797606129788036
0.4680076689877773 0.1868615772756182 0.6489884078156422 0.9715623197006016
0.4472598739535573 0.1870784379239747 0.6450318078649266 0.9633774801626918
0.4274202999759320 0.1870848608118418 0.6413410337376294 0.9552105168195424
0.4084511515061971 0.1868951837325424 0.6379018022332252 0.9470655204471414
0.3903159597266331 0.1865230395486411 0.6347006462660535 0.9389462727976324
0.3729795648034957 0.1859813770403239 0.6317248636456116 0.9308562675215122
0.3564080946568654 0.1852824828290845 0.6289624694233682 0.9227987295132662
0.3405689407884726 0.1844380040662867 0.6264021515258908 0.9147766328201244
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=1.0*ones(1,10);
inputs.PEsol.global_theta_min=0.0*ones(1,10);
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
inputs.nlpsol.eSS.log_var=1:(5+10); 
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
