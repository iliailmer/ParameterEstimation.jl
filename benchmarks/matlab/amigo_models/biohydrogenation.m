addpath(genpath('../src'))
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='DaisyEx3Model'; % Folder to keep results
inputs.pathd.short_name='Dex3';                 % To identify figures and reports
%======================
% MODEL RELATED DATA
%======================
clear
inputs.model.input_model_type='charmodelC';           % Model type- C
inputs.model.n_st=3;                                  % Number of states
inputs.model.n_par=6;                                 % Number of model parameters
%inputs.model.n_stimulus=0;                            % Number of inputs, stimuli or control variables
inputs.model.st_names=char('x4','x5', 'x6');    %x1=V, x2=R        % Names of the states
inputs.model.par_names=char('k5', 'k6', 'k7', 'k8', 'k9', 'k10');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli

inputs.model.eqns=char('dx4 = -k5 * x4 / (k6 + x4);', 'dx5 = k5 * x4 / (k6 + x4) - k7 * x5 / (k8 + x5 + x6);', 'dx6 = k7 * x5 / (k8 + x5 + x6) - k9 * x6 * (k10 - x6) / k10;'); % Equations describing system dynamics.
inputs.model.par = [1.0, 1.1, 1.1, 1.1, 1.1, 1.1];         % Nominal value for the parameters
% inputs.model.AMIGOsensrhs = 1;                       % Generate the sensitivity equations for exact
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[1.0, 1.0, 1.1];        % Initial conditions
inputs.exps.t_f{1}=20;                       % Experiments duration
inputs.exps.n_obs{1}=2;                       % Number of observables
inputs.exps.obs_names{1}=char('y1','y2'); % Names of the observables
inputs.exps.obs{1}=char('y1=x4','y2=x5');
inputs.exps.t_con{1}=[0 1];                 % Input swithching times including:
inputs.exps.n_s{1}=20;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
1.0 1.0
0.9772650316941635 1.0046655895507275
0.9548280505482635 1.0090842409763083
0.9326909927599267 1.0132626083421328
0.9108557054936655 1.0172077323325575
0.889323942038991 1.0209270950290972
0.8680973569743798 1.0244286816551698
0.8471775013519242 1.027721050694502
0.8265658179206239 1.03081341408483
0.8062636364076868 1.0337157295847201
0.7862721688773661 1.0364388079362472
0.76659250518763 1.038994438112768
0.7472256085655986 1.0413955348090624
0.7281723113229969 1.0436563134831578
0.7094333107341264 1.0457924997170307
0.6910091650972194 1.0478215817811696
0.6729002900017038 1.0497631179708886
0.6551069548243583 1.0516391138719476
0.6376292794752342 1.0534744901662192
0.6204672314132962 1.0552976688737516

];
inputs.PEsol.id_global_theta=char('k6', 'k5', 'k7');
inputs.PEsol.global_theta_max=2.*ones(1,3);
inputs.PEsol.global_theta_min=-1.*ones(1,3);
inputs.PEsol.id_global_theta_y0=char('x4','x5');               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=2.*ones(1,2);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=-1.*ones(1,2);
inputs.PEsol.id_local_theta{1}='none';                % [] 'all'|User selected| 'none' (default)
% inputs.PEsol.local_theta_max{iexp}=[];              % Maximum allowed values for the paramters
% inputs.PEsol.local_theta_min{iexp}=[];              % Minimum allowed values for the parameters
% inputs.PEsol.local_theta_guess{iexp}=[];            % [] Initial guess
inputs.PEsol.id_local_theta_y0{1}='none';             % [] 'all'|User selected| 'none' (default)
% inputs.PEsol.local_theta_y0_max{iexp}=[];           % Maximum allowed values for the initial conditions
% inputs.PEsol.local_theta_y0_min{iexp}=[];           % Minimum allowed values for the initial conditions
% inputs.PEsol.local_theta_y0_guess{iexp}=[];         % [] Initial guess


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
%  inputs.nlpsol.eSS.local.nl2sol.maxiter = 150;       % Parameters for local solver
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
