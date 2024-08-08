addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='hiv_7model'; % Folder to keep results
inputs.pathd.short_name='hiv_7';                 % To identify figures and reports
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
inputs.model.par = [0.561, 0.574, 0.558, 0.278, 0.862, 0.458, 0.777, 0.66, 0.338, 0.751];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.417, 0.805, 0.565, 0.805, 0.654];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = w', 'y2 = z', 'y3 = x', 'y4 = yy+vv');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.8050000000000000 0.6540000000000000 0.4170000000000000 1.3700000000000001
0.7856476130267719 0.6459623437718031 0.4261752042823289 1.3870295173633682
0.7670005318135612 0.6377560053817057 0.4347011757914956 1.4038501880171808
0.7490107642020628 0.6294132885414367 0.4426046072210708 1.4204813522959141
0.7316346431797566 0.6209633103061959 0.4499120698577295 1.4369405110416595
0.7148323985703742 0.6124322777393073 0.4566498882965644 1.4532434774389751
0.6985677735681285 0.6038437409547961 0.4628440321582009 1.4694045185023434
0.6828076813290800 0.5952188246736497 0.4685200234552696 1.4854364865581275
0.6675218973321211 0.5865764402177560 0.4737028582814072 1.5013509411233406
0.6526827836716436 0.5779334796768997 0.4784169415373856 1.5171582616268664
0.6382650418507663 0.5693049938141133 0.4826860334612344 1.5328677514466702
0.6242454910086076 0.5607043551236918 0.4865332067895291 1.5484877337528835
0.6106028688431334 0.5521434073210393 0.4899808134432766 1.5640256396528285
0.5973176527863104 0.5436326024220971 0.4930504597015360 1.5794880891318437
0.5843718992521658 0.5351811264621237 0.4957629888971879 1.5948809652755702
0.5717490990153465 0.5267970148064735 0.4981384707409325 1.6102094822459130
0.5594340469895994 0.5184872579191878 0.5001961964503423 1.6254782474658445
0.5474127248651476 0.5102578983771905 0.5019546789297151 1.6406913184483922
0.5356721952332522 0.5021141198477862 0.5034316573128973 1.6558522546834964
0.5242005059773790 0.4940603286841107 0.5046441052446472 1.6709641649736136
0.5129866038451676 0.4861002287363612 0.5056082423361494 1.6860297505855957
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
