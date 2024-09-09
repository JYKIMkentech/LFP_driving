clc; clear; close all;

%% Interface

data_folder = 'C:\Users\USER\Documents\GitHub\LFP_driving';
file_name = 'LFP_k1_0_05C_25degC.xlsx';
file_path = [data_folder filesep file_name];

save_path = data_folder;
sample_plot = 1;
max_capacity_Ah = 2.57; % Example maximum capacity from your previous calculation (Ah)

%% Engine

% Load the data from the Excel file
data_now = readtable(file_path);

% Extract relevant columns
data1.I = data_now.Current_A_;
data1.V = data_now.Voltage_V_;
data1.t = data_now.Test_Time_s_; % experiment time in seconds
data1.T = data_now.Surface_Temp_degC_; % surface temperature

% absolute current
data1.I_abs = abs(data1.I);

% Type assignment based on current values
data1.type = char(zeros([length(data1.t),1]));
data1.type(data1.I > 0) = 'C';
data1.type(data1.I == 0) = 'R';
data1.type(data1.I < 0) = 'D';

% Step assignment
data1_length = length(data1.t);
data1.step = zeros(data1_length,1);
m  = 1;
data1.step(1) = m;
for j = 2:data1_length
    if data1.type(j) ~= data1.type(j-1)
        m = m + 1;
    end
    data1.step(j) = m;
end

% Check for errors: if any step has more than one type
vec_step = unique(data1.step);
num_step = length(vec_step);
for i_step = 1:num_step
    type_in_step = unique(data1.type(data1.step == vec_step(i_step)));
    
    if size(type_in_step,1) ~= 1 || size(type_in_step,2) ~= 1
        disp('ERROR: step assignment is not unique for a step')
        return
    end
end

% Calculate cumulative charge (Q) using cumtrapz for numerical integration
data1.Q = zeros(data1_length, 1); % Initialize Q
data1.SOC = zeros(data1_length, 1); % Initialize SOC

% Identify the discharge step index
discharge_step = find(strcmp(cellstr(data1.type), 'D'), 1);

if ~isempty(discharge_step)
    % Set SOC to 1 at the start of discharge
    data1.SOC(discharge_step) = 1;

    % Calculate Q using cumulative trapezoidal integration
    Q_discharge = cumtrapz(data1.t(discharge_step:end)/3600, data1.I(discharge_step:end));
    data1.Q(discharge_step:end) = Q_discharge;
    
    % Calculate SOC based on Q
    data1.SOC(discharge_step:end) = 1 + data1.Q(discharge_step:end) / max_capacity_Ah;
end

% Plot for selected samples
if sample_plot
    figure
    title(file_name(1:end-5))
    hold on
    plot(data1.t/3600, data1.V, '-')
    xlabel('Time (hours)')
    ylabel('Voltage (V)')
    
    yyaxis right
    plot(data1.t/3600, data1.I, '-')
    ylabel('Current (C)')
end

% Create the output struct (data_line format)
data_line = struct('V', zeros(1,1), 'I', zeros(1,1), 't', zeros(1,1), 'indx', zeros(1,1), ...
    'type', char('R'), 'steptime', zeros(1,1), 'T', zeros(1,1), 'cycle', 0, 'Q', zeros(1,1), 'SOC', zeros(1,1));
data = repmat(data_line, num_step, 1);

% Fill in the struct
n = 1; 
for i_step = 1:num_step
    range = find(data1.step == vec_step(i_step));
    data(i_step).V = data1.V(range);
    data(i_step).I = data1.I(range);
    data(i_step).t = data1.t(range);
    data(i_step).indx = range;
    data(i_step).type = data1.type(range(1));
    data(i_step).steptime = data1.t(range) - data1.t(range(1));
    data(i_step).T = data1.T(range);
    data(i_step).cycle = 0; % No cycle information in the provided data
    data(i_step).Q = data1.Q(range);
    data(i_step).SOC = data1.SOC(range);
    
    % Display progress
    if i_step > num_step/10*n
        fprintf('%6.1f%%\n', round(i_step/num_step*100));
        n = n + 1;
    end
end

% Save output data
if ~isfolder(save_path)
    mkdir(save_path)
end
save_fullpath = [save_path filesep file_name(1:end-5) '.mat'];
save(save_fullpath, 'data')

