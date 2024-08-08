addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='hiv_4model'; % Folder to keep results
inputs.pathd.short_name='hiv_4';                 % To identify figures and reports
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
inputs.model.par = [0.227, 0.188, 0.625, 0.211, 0.257, 0.395, 0.757, 0.178, 0.77, 0.177];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.881, 0.475, 0.881, 0.584, 0.691];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = w', 'y2 = z', 'y3 = x', 'y4 = yy+vv');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.5840000000000000 0.6909999999999999 0.8810000000000000 1.3560000000000001
0.5691244410894108 0.6867848702654364 0.8603491194504859 1.3636418550516129
0.5546798507518967 0.6826277984500015 0.8407263337234340 1.3707007481725102
0.5406348914401946 0.6785224751509052 0.8220706697011689 1.3772173511040495
0.5269627770748309 0.6744631866609221 0.8043255434467297 1.3832290485498593
0.5136405724070653 0.6704447600915415 0.7874383859136094 1.3887702406032210
0.5006486011037127 0.6664625149564314 0.7713603053531594 1.3938726141188846
0.4879699445780013 0.6625122201919685 0.7560457823776332 1.3985653865501400
0.4755900168668813 0.6585900557840009 0.7414523941291431 1.4028755253346781
0.4634962034708077 0.6546925783231647 0.7275405644313748 1.4068279455347359
0.4516775541752229 0.6508166899337454 0.7142733371731784 1.4104456881079153
0.4401245215641635 0.6469596101187141 0.7016161704917304 1.4137500809036334
0.4288287383213747 0.6431188501433718 0.6895367496098870 1.4167608842281374
0.4177828275309124 0.6392921896429699 0.6780048164216929 1.4194964226101732
0.4069802411180408 0.6354776551920650 0.6669920141402464 1.4219737042069189
0.3964151223273143 0.6316735006148722 0.6564717455071521 1.4242085291277926
0.3860821887659650 0.6278781888502138 0.6464190432291230 1.4262155878086735
0.3759766330645298 0.6240903752126559 0.6368104514530379 1.4280085504423621
0.3660940386439493 0.6203088919143038 0.6276239172186013 1.4296001483601262
0.3564303084452145 0.6165327337305406 0.6188386909403417 1.4310022481616280
0.3469816047870143 0.6127610447086085 0.6104352350709010 1.4322259193040292
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
