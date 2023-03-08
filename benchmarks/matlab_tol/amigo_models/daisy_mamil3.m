addpath(genpath('../src'))
addpath(genpath("./"))
addpath(genpath("../samples/daisy_mamil3/"))
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
inputs.model.n_st=3;                                  % Number of states
inputs.model.n_par=5;                                 % Number of model parameters
%inputs.model.n_stimulus=0;                            % Number of inputs, stimuli or control variables
inputs.model.st_names=char('x1','x2','x3');    %x1=V, x2=R        % Names of the states
inputs.model.par_names=char('a12', 'a13', 'a21', 'a31', 'a01');             % Names of the parameters
%inputs.model.stimulus_names=char('light');  % Names of the stimuli
inputs.model.eqns=char('dx1 = -(a21 + a31 + a01) * x1 + a12 * x2 + a13 * x3;', 'dx2 = a21 * x1 - a12 * x2;', 'dx3 = a31 * x1 - a13 * x3;');                                 % Equations describing system dynamics.
inputs.model.par = [0.03, 0.02, 0.05, 0.03, 0.02];         % Nominal value for the parameters
% inputs.model.AMIGOsensrhs = 1;                       % Generate the sensitivity equations for exact
%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
% EXPERIMENT DESIGN
inputs.exps.n_exp=1;                          % Number of experiments
% EXPERIMENT 1
inputs.exps.exp_y0{1}=[1.0, 1.0, 1.0];        % Initial conditions
inputs.exps.t_f{1}=1;                       % Experiments duration
inputs.exps.n_obs{1}=2;                       % Number of observables
inputs.exps.obs_names{1}=char('y1','y2'); % Names of the observables
inputs.exps.obs{1}=char('y1=x1','y2=x2');
inputs.exps.t_con{1}=[0 1];                 % Input swithching times including:
% inputs.exps.n_s{1}=20;
inputs.exps.data_type='real';
% inputs.exps.exp_data{1}=[2.0 1.0
% 2.03952942072742 0.9276311719129872
% 2.063221883532521 0.8994348378831487
% 2.080799851429902 0.8947677729545105
% 2.0972158733694894 0.9026760235677073
% 2.114913303689913 0.917355779057597
% 2.135025223998584 0.9357399522024451
% 2.1580091372075536 0.9562150772012912
% 2.1839819715155784 0.9779392274173088
% 2.212896022034377 1.000479501924943
% 2.244630573775909 1.0236194211694714
% 2.279038919356856 1.0472566533738128
% 2.3159718618589618 1.0713487571870446
% 2.3552888930281184 1.0958844426164835
% 2.396862976081291 1.1208683880658927
% 2.440582067035102 1.1463132544093857
% 2.486349024393318 1.1722355161506313
% 2.5340807699641132 1.1986533132215766
% 2.583707146614094 1.2255853696080141
% 2.635169698746112 1.2530504723319549
% ];
inputs.PEsol.id_global_theta='all';
inputs.PEsol.global_theta_max=2.*ones(1,5);
inputs.PEsol.global_theta_min=-1.*ones(1,5);
inputs.PEsol.id_global_theta_y0='all';               % [] 'all'|User selected| 'none' (default)
inputs.PEsol.global_theta_y0_max=2.*ones(1,3);                % Maximum allowed values for the initial conditions
inputs.PEsol.global_theta_y0_min=-1*ones(1,3)
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
inputs.ivpsol.rtol=1.0e-10;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0e-10;

for i = 3:21
  tic
  start = cputime;
  inputs.exps.exp_data{1}= readmatrix(sprintf('daisy_mamil3-%d.csv', i));
  inputs.exps.exp_data{1}(:, 1) = [];

  inputs.exps.n_s{1}=size(inputs.exps.exp_data{1}, 1);
  AMIGO_Prep(inputs);
  AMIGO_PE(inputs);
  endTime = cputime - start;
  endtoc = toc;
  cputime_s = sprintf('CPU Time taken for %d is %d', i, endTime);
  elapsed_s = sprintf('Elapsed Time taken for %d is %d', i, endtoc);
  disp(cputime_s);
  disp(elapsed_s);
end

