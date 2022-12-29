addpath(genpath('../src'))
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='CRNModel'; % Folder to keep results
inputs.pathd.short_name='CRN';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=4;                                  % Number of states
inputs.model.n_par=7;                                 % Number of model parameters
%inputs.model.n_stimulus=0;                            % Number of inputs, stimuli or control variables
inputs.model.st_names=char('x1','x2','x3', 'x4');    %x1=V, x2=R        % Names of the states
inputs.model.par_names=char('k01', 'k12', 'k13', 'k14', 'k21', 'k31', 'k41');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char('dx1 = -k01 * x1 + k12 * x2 + k13 * x3 + k14 * x4 - k21 * x1 - k31 * x1 - k41 * x1;', 'dx2 = -k12 * x2 + k21 * x1;', 'dx3 = -k13 * x3 + k31 * x1;', 'dx4 = -k14 * x4 + k41 * x1;');                                 % Equations describing system dynamics.
inputs.model.par = [0.18, 0.21, 0.58, 0.64, -0.12, 0.9, 0.02];         % Nominal value for the parameters
% inputs.model.AMIGOsensrhs = 1;                       % Generate the sensitivity equations for exact
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[1.0, 2.2, 0.8, -1.0];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1','y2','y3','y4'); % Names of the observables
inputs.exps.obs{1}=char('y1=x1+x2', 'y2=x3+x4', 'y3=x1+x3', 'y4=x2+x4');
inputs.exps.t_con{1}=[0 10];                 % Input swithching times including:
inputs.exps.n_s{1}=20;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[3.0 0.0 2.0 1.0
2.4149170746088684 0.49254090103636394 2.006248995341815 0.9012089803034171
2.0337914416363048 0.794725460186959 2.0351526247164964 0.7933642771067672
1.7550483486985677 0.9987654116159914 2.0742447200384033 0.6795690402761561
1.533444211750992 1.1463992909759586 2.115188150637521 0.5646553520894297
1.34814793239913 1.2573828483741372 2.1526597706671575 0.45287101010610964
1.1889544250391981 1.3418125154808664 2.1834787928443147 0.3472881476757495
1.0503353659135324 1.4054578785164402 2.2059529985374002 0.24984024589257242
0.928873639088659 1.4520732591548455 2.2194003948777303 0.16154650336577425
0.8221562259419125 1.4844051838181562 2.223805229525575 0.08275618023449369
0.7282941690873559 1.5046374827000906 2.219574080848242 0.01335757093920445
0.6457112433820796 1.5145940606597004 2.207365233351038 -0.047059929309257904
0.5730472960249151 1.515837148588883 2.1879710988767043 -0.09908665426290585
0.5091103413587413 1.5097201993444456 2.1622386830882734 -0.14340814238508634
0.4528494018681479 1.4974207248544344 2.1310171111491343 -0.1807469844265518
0.40333631980596907 1.4799638372505186 2.095124230584699 -0.21182407352821156
0.3597516658777004 1.4582410535178252 2.055326527901702 -0.23733380850617605
0.3213727708685027 1.433026298320708 2.0123282197932344 -0.2579291506040238
0.28756310045472855 1.4049899430076913 1.9667665626950295 -0.2742135192326098
0.2577626685404456 1.3747112700665123 1.9192112817985731 -0.2867373431916155
];
inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=2.*ones(1,7);
inputs.PEsol.global_theta_min=-2.*ones(1,7);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=2.*ones(1,4);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=-2*ones(1,4)
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
%  inputs.nlpsol.eSS.log_var = 1:3;                    % Index of parameters to be considered in log scale
 inputs.nlpsol.eSS.maxeval = 20000;                  % Maximum number of cost function evaluations
 inputs.nlpsol.eSS.maxtime = 600;                    % Maximum time spent for optimization
 inputs.nlpsol.eSS.local.solver = 'nl2sol';
 inputs.nlpsol.eSS.local.finish = 'nl2sol';
%  inputs.nlpsol.eSS.local.nl2sol.maxiter = 1500;       % Parameters for local solver
%  inputs.nlpsol.eSS.local.nl2sol.maxfeval = 200;
  inputs.nlpsol.eSS.local.nl2sol.display = 1;
%  inputs.nlpsol.eSS.local.nl2sol.objrtol = 1e-6;
%  inputs.nlpsol.eSS.local.nl2sol.tolrfun = 1e-5;
% %
% inputs.exps.u_interp{1}='sustained';          % Stimuli definition for experiment 1
                                              % Initial and final time
%inputs.exps.u{1}=1;                           % Values of the inputs for exp 1
AMIGO_Prep(inputs);
AMIGO_PE(inputs);
