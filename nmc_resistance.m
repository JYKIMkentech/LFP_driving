clc; clear; close all;

%% Interface

data_file = 'HNE_Formation_1_13.txt'; % Replace with the actual path if necessary
I_1C = 55.6; %[A] Given total capacity
n_hd = 1; % headline number used in 'readtable' option
sample_plot = 1;

%% Engine

% Load the data from the text file
data_now = readtable(data_file, 'FileType', 'text', 'Delimiter', '\t', 'NumHeaderLines', 1);

% Extract the relevant columns using Var variables
data1.I = data_now.Var9; % Current (A)
data1.V = data_now.Var8; % Voltage (V)
data1.t2 = data_now.Var7; % Experiment time
data1.t1 = data_now.Var6; % Step time
data1.cycle = data_now.Var14; % Cycle index
data1.T = data_now.Var21; % Temperature (°C)

% Check the format of t2 (Total Time), assuming it is in seconds
% If it is already in numeric format, no conversion is needed
if isnumeric(data1.t2)
    data1.t = data1.t2; % Use directly if t2 is already in seconds
else
    % If t2 is in string or other format, convert it properly
    % Example: data1.t2 might be a duration string 'HH:MM:SS'
    try
        data1.t = seconds(duration(data1.t2, 'InputFormat', 'hh:mm:ss'));
    catch
        error('t2 format is not recognized. Please check the time format.');
    end
end

% Absolute current
data1.I_abs = abs(data1.I);

% Define the type of operation: Charge (C), Rest (R), Discharge (D)
data1.type = char(zeros([length(data1.t), 1]));
data1.type(data1.I > 0) = 'C';
data1.type(data1.I == 0) = 'R';
data1.type(data1.I < 0) = 'D';

% Assign step numbers
data1_length = length(data1.t);
data1.step = zeros(data1_length, 1);
m = 1;
data1.step(1) = m;
for j = 2:data1_length
    if data1.type(j) ~= data1.type(j - 1)
        m = m + 1;
    end
    data1.step(j) = m;
end

% Check for errors, if any step has more than one type
vec_step = unique(data1.step);
num_step = length(vec_step);
for i_step = 1:num_step
    type_in_step = unique(data1.type(data1.step == vec_step(i_step)));
    if size(type_in_step, 1) ~= 1 || size(type_in_step, 2) ~= 1
        disp('ERROR: step assignment is not unique for a step');
        return;
    end
end

% Plot for selected samples
if sample_plot
    figure
    title('Voltage and Current vs Time');
    hold on
    plot(data1.t / 3600, data1.V, '-')
    xlabel('Time (hours)')
    ylabel('Voltage (V)')
    
    yyaxis right
    plot(data1.t / 3600, data1.I / I_1C, '-')
    ylabel('Current (C)')
    grid on;
end

% Save parsed data into a .mat file
save_path = fullfile(pwd, 'parsed_data.mat');
save(save_path, 'data1');

