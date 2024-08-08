addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='hiv_9model'; % Folder to keep results
inputs.pathd.short_name='hiv_9';                 % To identify figures and reports
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
inputs.model.par = [0.573, 0.559, 0.623, 0.622, 0.445, 0.817, 0.394, 0.449, 0.814, 0.745];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.663, 0.18, 0.836, 0.671, 0.899];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = w', 'y2 = z', 'y3 = x', 'y4 = yy+vv');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.6710000000000000 0.8990000000000000 0.6630000000000000 1.0160000000000000
0.6447348719617880 0.8671883879086489 0.6563419127913270 0.9976847009174008
0.6195101636276522 0.8365563941106822 0.6505996340365398 0.9797058182307550
0.5952823748768049 0.8070525944366616 0.6456892902247593 0.9620760207007752
0.5720104421349226 0.7786287185103390 0.6415355003600179 0.9448040489794886
0.5496554539934596 0.7512393637136527 0.6380703754003849 0.9278953699233899
0.5281804234738571 0.7248417428286794 0.6352326518508881 0.9113527292056860
0.5075501041551191 0.6993954607211819 0.6329669395753755 0.8951766184261655
0.4877308403246590 0.6748623161539374 0.6312230671230076 0.8793656701925662
0.4686904435642591 0.6512061254140133 0.6299555105303484 0.8639169924108938
0.4503980899118380 0.6283925649317160 0.6291228937654179 0.8488264511810988
0.4328242330704141 0.6063890304800537 0.6286875508114502 0.8340889101739395
0.4159405301618551 0.5851645108873241 0.6286151409123983 0.8196984331060224
0.3997197773162204 0.5646894744850119 0.6288743097771026 0.8056484548836469
0.3841358530007367 0.5449357667568270 0.6294363906044143 0.7919319261159119
0.3691636674679459 0.5258765178620087 0.6302751396888532 0.7785414349701851
0.3547791170706640 0.5074860588817851 0.6313665021202782 0.7654693097347183
0.3409590424771421 0.4897398457883067 0.6326884037289000 0.7527077049428255
0.3276811900409392 0.4726143902637032 0.6342205659655034 0.7402486734840701
0.3149241757517723 0.4560871966077058 0.6359443408653755 0.7280842267661269
0.3026674513261964 0.4401367040671222 0.6378425636329833 0.7162063846859945
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
