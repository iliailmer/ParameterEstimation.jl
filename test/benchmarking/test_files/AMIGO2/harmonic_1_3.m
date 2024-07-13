addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='harmonic_1model'; % Folder to keep results
inputs.pathd.short_name='harmonic_1';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=2;                                  % Number of states:\\\
inputs.model.n_par=2;                                 % Number of model parameters
inputs.model.st_names=char('x1', 'x2');    % Names of the states
inputs.model.par_names=char('a', 'b');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char( 'dx1 = -a*x2;',  'dx2 = +x1/b;');               % Equations describing system dynamics.
inputs.model.par = [0.439, 0.617];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.45, 0.813];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=2;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = x1', 'y2 = x2');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
 0.4500000000000000 0.8129999999999999
 0.4317597764557221 0.8487330025755986
 0.4127516661711303 0.8829565296180878
 0.3930094751574514 0.9156097144278276
 0.3725683149918278 0.9466344831660748
 0.3514645403712587 0.9759756581394966
 0.3297356844555895 1.0035810559339819
 0.3074203921146054 1.0294015802230949
 0.2845583511979588 1.0533913090861071
 0.2611902219501484 1.0755075766803419
 0.2373575646960889 1.0957110491225694
 0.2131027659258841 1.1139657944444985
 0.1884689629102591 1.1302393464979481
 0.1634999669807216 1.1445027626960431
 0.1382401856109214 1.1567306754877316
 0.1127345434377420 1.1669013374740971
 0.0870284023626246 1.1749966600862076
 0.0611674808752236 1.1810022457557150
 0.0351977727428746 1.1849074135209934
 0.0091654652104627 1.1867052180232820
-0.0168831431438014 1.1863924618590369
];


inputs.ivpsol.rtol=1.0e-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-13;

inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=3.0*ones(1,2);
inputs.PEsol.global_theta_min=0.0*ones(1,2);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=3.0*ones(1,2);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=0.0*ones(1,2);
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
inputs.nlpsol.eSS.log_var=1:(2+2); 
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
