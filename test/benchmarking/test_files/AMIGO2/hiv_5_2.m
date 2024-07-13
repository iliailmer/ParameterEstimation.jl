addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='hiv_5model'; % Folder to keep results
inputs.pathd.short_name='hiv_5';                 % To identify figures and reports
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
inputs.model.par = [0.131, 0.326, 0.196, 0.337, 0.195, 0.354, 0.431, 0.151, 0.654, 0.553];         % Nominal value for the parameters
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
inputs.exps.obs{1}=char( 'y1 = w', 'y2 = z', 'y3 = x', 'y4 = yy+vv');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.5610000000000001 0.8430000000000000 0.3120000000000000 0.6940000000000000
0.5439256388873058 0.8209232253054490 0.3129180685251345 0.6877683584551655
0.5273613703351790 0.7994068554586351 0.3138141078751096 0.6815684847876849
0.5112922274257974 0.7784375189730779 0.3146888145902325 0.6754009711479371
0.4957036833113344 0.7580021207923006 0.3155428666190820 0.6692663781958366
0.4805816371757404 0.7380878388570397 0.3163769236211448 0.6631652362032368
0.4659124007711342 0.7186821205720915 0.3171916272761637 0.6570980461160842
0.4516826854955219 0.6997726791849466 0.3179876015991671 0.6510652805779948
0.4378795899808206 0.6813474900875750 0.3187654532602192 0.6450673849168594
0.4244905881629641 0.6633947870521122 0.3195257719079803 0.6391047780959738
0.4115035178068842 0.6459030584100693 0.3202691304962614 0.6331778536311866
0.3989065694620324 0.6288610431845139 0.3209960856127779 0.6272869804754384
0.3866882758253549 0.6122577271836492 0.3217071778093918 0.6214325038720184
0.3748375014899796 0.5960823390636930 0.3224029319331794 0.6156147461778432
0.3633434330603409 0.5803243463686841 0.3230838574576908 0.6098340076579221
0.3521955696143513 0.5649734515537113 0.3237504488138590 0.6040905672522330
0.3413837134962212 0.5500195879983701 0.3244031857200034 0.5983846833160535
0.3308979614229607 0.5354529160158835 0.3250425335104779 0.5927165943348703
0.3207286958903160 0.5212638188638045 0.3256689434624924 0.5870865196148211
0.3108665768633982 0.5074428987609074 0.3262828531207325 0.5814946599496881
0.3013025337395620 0.4939809729153583 0.3268846866193815 0.5759411982653229
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=2.0*ones(1,10);
inputs.PEsol.global_theta_min=0.0*ones(1,10);
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
