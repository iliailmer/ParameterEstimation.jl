addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='hiv_6model'; % Folder to keep results
inputs.pathd.short_name='hiv_6';                 % To identify figures and reports
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
inputs.model.par = [0.355, 0.634, 0.205, 0.673, 0.332, 0.247, 0.569, 0.116, 0.763, 0.104];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.642, 0.316, 0.688, 0.87, 0.299];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = w', 'y2 = z', 'y3 = x', 'y4 = yy+vv');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.8700000000000000 0.2990000000000000 0.6420000000000000 1.0040000000000000
0.8413396996610073 0.2983306856321507 0.6350181245739387 0.9946791851328636
0.8135027191404470 0.2976194731743931 0.6283225478860796 0.9854315661075570
0.7864751528077271 0.2968688428990563 0.6219023128148252 0.9762576793717302
0.7602424259512479 0.2960811600277993 0.6157468771421125 0.9671579919355434
0.7347894037545489 0.2952586786953025 0.6098460984384384 0.9581329043898064
0.7101004905747826 0.2944035459273149 0.6041902194322916 0.9491827538648556
0.6861597202002220 0.2935178056109532 0.5987698538543004 0.9403078169258384
0.6629508377313206 0.2926034024381516 0.5935759727460214 0.9315083124007788
0.6404573736966787 0.2916621858058533 0.5885998912230525 0.9227844041384260
0.6186627109822056 0.2906959136589702 0.5838332556819507 0.9141362036934568
0.5975501451188975 0.2897062562643207 0.5792680314402799 0.9055637729371344
0.5771029384423596 0.2886947999057191 0.5748964907989997 0.8970671265919802
0.5573043686056856 0.2876630504921313 0.5707112015163452 0.8886462346894555
0.5381377718967675 0.2866124370723785 0.5667050156823091 0.8803010249500252
0.5195865817816682 0.2855443152512466 0.5628710589828434 0.8720313850853193
0.5016343630674325 0.2844599705030822 0.5592027203429342 0.8638371650224190
0.4842648420507847 0.2833606213800385 0.5556936419377493 0.8557181790505630
0.4674619329934467 0.2822474226130750 0.5523377095611510 0.8476742078908208
0.4512097612404603 0.2811214681046380 0.5491290433409841 0.8397050006894833
0.4354926832749689 0.2799837938126742 0.5460619887906535 0.8318102769361192
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=3.0*ones(1,10);
inputs.PEsol.global_theta_min=0.0*ones(1,10);
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
