%% Below you find a template script, allowing to run parameter estimation
% using command-line arguments. It is constructed in such a way that you can
% use the "cell-mode" of the MATLAB editor, increasing the ease of use.
% The cell-model can be enabled by choosing in the menu "Cell->Enable Cell Mode".
% Once enabled, you can execute the yellow cells in this script (assuming you
% opened it in the MATLAB editor) by selecting them and pressing "Ctrl-Enter".
clc; clear all;close all
run ../../src/'IQMtools V1.2.2.2'/IQMlite/installIQMlite.m;
run ../../src/'IQMtools V1.2.2.2'/IQMpro/installIQMpro.m;
run ../../src/'IQMtools V1.2.2.2'/installIQMtools.m;

%% LOAD THE PROJECT
lv = IQMprojectSB('project');

%% DISPLAY INFORMATION ABOUT THE PROJECT
IQMinfo(lv);

%% KEEP THE ORIGINAL PROJECT UNCHANGED
lvopt = lv;

%% COMPARE MEASUREMENTS WITH MODEL
IQMcomparemeasurements(lv);

%% SELECT PARAMETERS/STATES TO ESTIMATE AND CORRESPONDING BOUNDS
% Global parameters
% Names         Lower bounds  Upper bounds
paramdata = {
'k1'             0.0001        .1
'k2'             0.0001        .1
'k3'             0.0001        .1
'k4'             0.0001        .1
'k5'             0.0001        .1
'k6'             0.0001        .1
% 'kRGact'       0.114811      11.4811
% 'k1Gact'       999.422       99942.2
% 'k2Gact'       0.419216      41.9216
% 'kGactPDEact'  0.0365684     3.65684
% 'kRArr1'       0.0101633     1.01633
% 'kRArr2'       0.0405202     4.05202
% 'kGr1'         0.00175855    0.175855
% 'kGr2'         0.23023       23.023
% 'kG'           0.238471      23.8471
% 'magStim'      0.2           20
% 'durStim'      0.01          1
};

% Local (experiment dependend) parameters
% Names         Lower bounds  Upper bounds
paramdatalocal = {
};

% Initial conditions (always experiment dependend)
% Names         Lower bounds  Upper bounds
icdata = {
'x1'             0.5             1.5
'x2'             0.5             1.5
'x3'             0.5             1.5
'x4'             0.5             1.5
'x5'             0.5             1.5
'x6'             0.5             1.5
% 'Arr'          0.5           50
% 'G'            300           30000
% 'Gact'         0             100
% 'GactPDEact'   0             100
% 'Gr'           0             100
% 'PDE'          10            1000
% 'R'            50            5000
% 'Ract'         0             100
% 'RactArr'      0             100
% 'RactG'        0             100
};


%% DEFINE THE ESTIMATION INFORMATION (STRUCTURE)
estimation = [];

% Model and experiment settings
estimation.modelindex = 1;
estimation.experiments.indices = [1];%, 2, 3, 4];
estimation.experiments.weight = [1];%, 1, 1, 1];

% Optimization settings
estimation.optimization.method = 'simplexIQM';
estimation.optimization.options.maxfunevals = 2000;

% Integrator settings
estimation.integrator.options.abstol = 1e-006;
estimation.integrator.options.reltol = 1e-006;
estimation.integrator.options.minstep = 0;
estimation.integrator.options.maxstep = Inf;
estimation.integrator.options.maxnumsteps = 1000;

% Flags
estimation.displayFlag = 2; % show iterations and final message
estimation.scalingFlag = 2; % scale by mean values
estimation.timescalingFlag = 0; % do not apply time-scaling
estimation.initialconditionsFlag = 1; % do use initial conditions from measurement data (if available)

% Always needed (do not change ... unless you know what you do)
estimation.parameters = paramdata;
estimation.parameterslocal = paramdatalocal;
estimation.initialconditions = icdata;

% Run estimation
output = IQMparameterestimation(lvopt,estimation)
% Get optimized project
lvopt = output.projectopt;

%% COMPARE OPTIMIZED PROJECT WITH MEASUREMENTS
IQMcomparemeasurements(lvopt,estimation.modelindex);
display(lvopt);
% %% ANALYSIS OF RESIDUALS
% IQManalyzeresiduals(lvopt,estimation);

% %% RUN A-POSTERIORI IDENTIFIABILITY ANALYSIS (only considering global variables)
% IQMidentifiability(lvopt,paramdata(:,1));

%% RUN SOME FIT ANALYSIS
% (after completion click in lower figure to remove outliers, corresponding
%  to local minima. Finish with "Enter")
% output = IQMparameterfitanalysis(lvopt,estimation);

%% FITANALYSIS EVALUATION
% IQMfaboxplot(output);
% IQMfahist(output);
% IQMfacorr(output);
% IQMfaclustering(output);
% IQMfadetcorr(output);
% IQMfasigncorr(output);
