addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='hiv_2model'; % Folder to keep results
inputs.pathd.short_name='hiv_2';                 % To identify figures and reports
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
inputs.model.par = [0.312, 0.719, 0.465, 0.555, 0.115, 0.594, 0.59, 0.594, 0.855, 0.645];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.388, 0.45, 0.658, 0.148, 0.633];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = w', 'y2 = z', 'y3 = x', 'y4 = yy+vv');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.1480000000000000 0.6330000000000000 0.3880000000000000 1.1080000000000001
0.1414180226845196 0.6140257998836990 0.3838974825524164 1.0847078216845851
0.1351273320798774 0.5955885050075047 0.3801419804768162 1.0618445397354366
0.1291158653257322 0.5776746299808505 0.3767097165163073 1.0394090343107505
0.1233719173127756 0.5602709631740170 0.3735786575396553 1.0173993987633820
0.1178841513467685 0.5433645586533127 0.3707283710775440 0.9958130358257992
0.1126416059522960 0.5269427295277814 0.3681398949589438 0.9746467437585824
0.1076336983647894 0.5109930424347815 0.3657956187410042 0.9538967935223894
0.1028502252013754 0.4955033129341880 0.3636791757651083 0.9335589979162976
0.0982813607472349 0.4804616016179516 0.3617753447988956 0.9136287735200443
0.0939176532438502 0.4658562107732250 0.3600699603320695 0.8941011961872560
0.0897500195195348 0.4516756814645106 0.3585498306933350 0.8749710507545119
0.0857697382609154 0.4379087909235846 0.3572026632422719 0.8562328755596342
0.0819684421860486 0.4245445501557802 0.3560169959664377 0.8378810022992944
0.0783381093458861 0.4115722016881676 0.3549821348830591 0.8199095916993806
0.0748710537503917 0.3989812173995515 0.3540880967056247 0.8023126654215925
0.0715599154884869 0.3867612963843697 0.3533255562897719 0.7850841345854767
0.0683976504869991 0.3749023628128560 0.3526857984215537 0.7682178252454868
0.0653775200325345 0.3633945637584437 0.3521606735542847 0.7517075011275982
0.0624930801615110 0.3522282669706101 0.3517425571388790 0.7355468838986718
0.0597381710071838 0.3413940585773748 0.3514243122271774 0.7197296712138568
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
