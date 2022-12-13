% Copyright 2018 Jake Alan Pitt and Julio R. Banga

% This file is part of GEARS.

% GEARS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% GEARS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with GEARS.  If not, see <http://www.gnu.org/licenses/>.

% File author: Jake Alan Pitt (jp00191.su@gmail.com)

% addpath(genpath("./"))



%% Model information

Model_information.Model_name = 'FHN'; % The name of the model. (string)

Model_information.Diff_eqns = {'dV = g*(V - V^3/3 + R)'; 'dR = (-1/g)*(V - a + b*R)'}; % The model equations. (cell array)

Model_information.Var_names = {'V', 'R'}; % The variable names. (cell array)

Model_information.Fitting_params = {'a', 'b', 'g'}; % The parameters that should be estimated. (cell array)

% Model_information.Fixed_params = []; % The parameters that have known values. (cell array) [Optional]
% All parameters are fit here.

Model_information.Param_bounds_lower = 10^-5*ones(1, length(Model_information.Fitting_params)); % The lower bounds for the parameter estimation. (vector)

Model_information.Param_bounds_upper = 10^5*ones(1, length(Model_information.Fitting_params)); % The upper bounds for the parameter estimation. (vector)


%% Fitting data

Fitting_exp_name = 'exp1'; % The name of the fitting experiment. (string)

Fitting_data.(char(Fitting_exp_name)).Observables = {'V'}; % The observables, should be either a subsection of the variable names or {'Custom_function'} for observed functions. (cell array)

Fitting_data.(char(Fitting_exp_name)).Time_points = [0.; 0.07500000000; 0.1500000000; 0.225000000; 0.3000000000; 0.3750000000; 0.4500000000; 0.5250000000; 0.6000000000; 0.6750000000; 0.7500000000; 0.8250000000; 0.9000000000; 0.9750000000; 1.050000000; 1.125000000; 1.200000000]; %[1.0000;  4.8000; 8.6000; 12.4000; 16.2000; 20.0000]; % The time points for the experiment. (vector)

Fitting_data.(char(Fitting_exp_name)).Measurements = [-1.; -0.947108404127384; -0.888048828846400; -0.822073987730191; -0.748170856386973; -0.665003889913865; -0.570840882888612; -0.463472658537472; -0.340144058248838; -0.197542464040998; -0.0319585250114381; 0.160173667732012; 0.381023346358266; 0.629086206140822; 0.896167264554509; 1.16547182841767; 1.41410548672615]; %[1.4480; 1.0070; -1.3020; 1.6260; -1.3930; 1.8890]; % The measurements, should be of the size # time points vs # observables. (matrix)

Fitting_data.(char(Fitting_exp_name)).Standard_deviation = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; % The standard deviation for the measurements, should be same size as the measurements. (matrix)

Fitting_data.(char(Fitting_exp_name)).Initial_conditions = [-1 1]; % The initial conditions for the states. Use a cell string for initial conditions that depend on parameters. (vector or cell array)


%% Cross-validation data

% Initial_exp
Validation_exp_name = 'Initial_exp'; % The name of the fitting experiment. (string)

Validation_data.(char(Validation_exp_name)).Observables = {'V'}; % The observables, should be either a subsection of the variable names or {'Custom_function'} for observed functions. (cell array)

Validation_data.(char(Validation_exp_name)).Time_points = [1.0000; 7.8000; 14.6000; 21.4000; 28.2000; 35.0000]; % The time points for the experiment. (vector)

Validation_data.(char(Validation_exp_name)).Measurements = [-1.0670; -2.0690; 1.1710; 1.7200; -1.1620; -2.0670]; % The measurements, should be of the size # time points vs # observables. (matrix)

Validation_data.(char(Validation_exp_name)).Standard_deviation = [0.1596; 0.2194; 0.1629; 0.2137; 0.1451; 0.2149]; % The standard deviation for the measurements, should be same size as the measurements. (matrix)

Validation_data.(char(Validation_exp_name)).Initial_conditions = [-7.4850; 0.1792]; % The initial conditions for the states. Use a cell string for initial conditions that depend on parameters. (vector or cell array)

% Follow_up_exp
Validation_exp_name = 'Follow_up_exp'; % The name of the fitting experiment. (string)

Validation_data.(char(Validation_exp_name)).Observables = {'V'}; % The observables, should be either a subsection of the variable names or {'Custom_function'} for observed functions. (cell array)

Validation_data.(char(Validation_exp_name)).Time_points = [1.0000; 7.8000; 14.6000; 21.4000; 28.2000; 35.0000]; % The time point for the experiment. (vector)

Validation_data.(char(Validation_exp_name)).Measurements = [1.9900; -1.7880; 1.3180; 1.9030; -1.0800; -2.2870]; % The measurements, should be of the size # time points vs # observables. (matrix)

Validation_data.(char(Validation_exp_name)).Standard_deviation = [0.2599; 0.2179; 0.1667; 0.2162; 0.1514; 0.2175]; % The standard deviation for the measurements, should be same size as the measurements. (matrix)

Validation_data.(char(Validation_exp_name)).Initial_conditions = [-8.7420; 2.3690]; % The initial conditions for the states. Use a cell string for initial conditions that depend on parameters. (vector or cell array)


%% Results folder

Results_folder = 'FHN_example_results'; % The results folder. (string)


%% Run GEARS

[Results, Processed_fitting_data, Processed_validation_data] = Run_GEARS(Results_folder, Model_information, Fitting_data, Validation_data, 'FHN_example_options'); % Run GEARS

