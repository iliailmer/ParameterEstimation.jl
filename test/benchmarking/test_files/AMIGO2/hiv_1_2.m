addpath(genpath('AMIGO2_PATH')) % Replace the string with path to the package
addpath(genpath("./"))
%======================
% PATHS RELATED DATA
%======================
inputs.pathd.results_folder='hiv_1model'; % Folder to keep results
inputs.pathd.short_name='hiv_1';                 % To identify figures and reports
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
inputs.model.par = [0.17, 0.116, 0.766, 0.723, 0.796, 0.883, 0.739, 0.469, 0.724, 0.195];         % Nominal value for the parameters
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[0.612, 0.215, 0.856, 0.517, 0.432];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=4;                       % Number of observables
inputs.exps.obs_names{1}=char('y1', 'y2', 'y3', 'y4'); % Names of the observables
inputs.exps.obs{1}=char( 'y1 = w', 'y2 = z', 'y3 = x', 'y4 = yy+vv');
inputs.exps.t_con{1}=[-0.5, 0.5];                 % Input swithching times including:
inputs.exps.n_s{1}=21;
inputs.exps.data_type='real';
inputs.exps.exp_data{1}=[
0.5170000000000000 0.4320000000000000 0.6120000000000000 1.0710000000000000
0.4991713767505533 0.4297434056347244 0.5975001946979785 1.0541391919025769
0.4819262223401090 0.4275329993554783 0.5841528147007989 1.0374616576742266
0.4652434126526222 0.4253528972409336 0.5718539978717804 1.0209900365053317
0.4491039947536127 0.4231898765827715 0.5605115140577890 1.0047416850944051
0.4334906639258374 0.4210329998529662 0.5500432343518911 0.9887296365493580
0.4183873704498324 0.4188732920135101 0.5403758297839700 0.9729633887091892
0.4037790252349675 0.4167034628585465 0.5314436610907848 0.9574495531912814
0.3896512807762681 0.4145176675805815 0.5231878282537663 0.9421923904149156
0.3759903694655096 0.4123112999355854 0.5155553541309056 0.9271942510398544
0.3627829854914449 0.4100808133226448 0.5084984810418246 0.9124559404214836
0.3500161997667546 0.4078235658512461 0.5019740628321213 0.8979770196081037
0.3376773997670075 0.4055376860833338 0.4959430379164657 0.8837560539339885
0.3257542480449689 0.4032219566404513 0.4903699712245889 0.8697908182659752
0.3142346546286184 0.4008757132815206 0.4852226549578581 0.8560784663453026
0.3031067596258614 0.3984987574020626 0.4804717596931137 0.8426156703534164
0.2923589232203201 0.3960912801946615 0.4760905287132596 0.8293987357603654
0.2819797209091841 0.3936537969539116 0.4720545095551792 0.8164236956396826
0.2719579423503023 0.3911870902151564 0.4683413176876685 0.8036863879166647
0.2622825925855304 0.3886921605915333 0.4649304279994933 0.7911825184279537
0.2529428947168799 0.3861701843238815 0.4618029904203329 0.7789077121847061
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
