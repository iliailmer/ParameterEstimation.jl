addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='hiv_8model'; % Folder to keep results
inputs.pathd.short_name='hiv_8';                 % To identify figures and reports
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
inputs.model.par = [0.68, 0.501, 0.865, 0.615, 0.439, 0.585, 0.115, 0.341, 0.628, 0.332];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.594, 0.443, 0.208, 0.339, 0.556];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = w', 'y2 = z', 'y3 = x', 'y4 = yy+vv');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.3390000000000000 0.5560000000000000 0.5940000000000000 0.6510000000000000
0.3287362482296451 0.5471315259136722 0.6074990843866807 0.6464524796591979
0.3187899778536173 0.5383955063794711 0.6204575560659987 0.6420847002036607
0.3091506258454814 0.5297905929433178 0.6328968814124905 0.6378922481069481
0.2998080704281611 0.5213153998201919 0.6448379749311434 0.6338704275700541
0.2907526024567475 0.5129685088950844 0.6563011659969396 0.6300143294061712
0.2819748995131980 0.5047484743419504 0.6673061737546606 0.6263188910632036
0.2734660024238559 0.4966538268829588 0.6778720891822232 0.6227789486798616
0.2652172939395529 0.4886830777097890 0.6880173634085761 0.6193892820023363
0.2572204793446791 0.4808347220880681 0.6977598014609991 0.6161446529216293
0.2494675687860451 0.4731072426653033 0.7071165606958529 0.6130398383273841
0.2419508611346641 0.4654991125018580 0.7161041532411533 0.6100696579129499
0.2346629292138428 0.4580087978436726 0.7247384518486237 0.6072289975087306
0.2275966062453189 0.4506347606545590 0.7330346986171012 0.6045128284667692
0.2207449733817327 0.4433754609249966 0.7410075161083098 0.6019162235690838
0.2141013482085620 0.4362293587734728 0.7486709204304159 0.5994343698855100
0.2076592741119558 0.4291949163555182 0.7560383359143575 0.5970625789636589
0.2014125104208443 0.4222705995947043 0.7631226110529513 0.5947962946939178
0.1953550232422900 0.4154548797490157 0.7699360354138530 0.5926310991560315
0.1894809769183974 0.4087462348251914 0.7764903572743905 0.5905627167207592
0.1837847260416448 0.4021431508527920 0.7827968017591855 0.5885870166498621
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
