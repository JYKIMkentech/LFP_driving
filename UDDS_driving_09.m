clc; clear; close all;

%% Parameters and File Paths
% File paths
file_paths = struct( ...
    'udds', 'C:\Users\deu04\OneDrive\바탕 화면\ECM_github\LFP_driving\uddscol.txt', ...
    'hwycol', 'C:\Users\deu04\OneDrive\바탕 화면\ECM_github\LFP_driving\hwycol.txt');

% Physical constants for Power model (Tesla model 3)
a = 34.98 * 4.44822; % lbf to Newton
b = 0.08650 * 4.44822 / 0.44704; % lbf/mph to N/(m/s)
c = 0.014800 * 4.44822 / 0.44704^2; % lbf/mph^2 to N/(m/s)^2
m_vehicle = 1927.768; % vehicle mass in [kg]
epsilon = 1.05;

% Battery pack configuration
m_series = 106; % Number of cells in series
n_parallel = 1; % Number of parallel strings

% Battery parameters for a single cell
OCV_cell = 3.2; % [V] 
R_cell = 0.0009; % Resistance [ohm]
nominal_capacity_Ah = 161; % [Ah]
Scaling_nominal_capacity_Ah = 16; % [Ah]

%% File Selection
disp('Select the file to analyze:');
disp('1. UDDS (uddscol.txt)');
disp('2. HWYCOL (hwycol.txt)');
file_choice = input('Enter the number of the file you want to process (1 or 2): ');

if file_choice == 1
    file_path = file_paths.udds;
    disp('You have selected UDDS.');
elseif file_choice == 2
    file_path = file_paths.hwycol;
    disp('You have selected HWYCOL.');
else
    error('Invalid selection. Please enter 1 or 2.');
end

%% Data Loading and Preprocessing
% Read the file
data_unit = readtable(file_path, 'Delimiter', '\t');
data_unit.Properties.VariableNames{1} = 'time';
data_unit.Properties.VariableNames{2} = 'speed_mph';

% Separate the first column as time and the second as speed
time = data_unit.time;
speed_mph = data_unit.speed_mph;

% Convert speed from mph to m/s (1 mph = 0.44704 m/s)
speed_ms = speed_mph * 0.44704;

%% Acceleration, Distance Calculation
acceleration = zeros(size(speed_ms));

% Central difference for interior points
for i = 2:length(time)-1
    acceleration(i) = (speed_ms(i+1) - speed_ms(i-1)) / (time(i+1) - time(i-1));
end

% Forward and Backward difference for the first and last points
acceleration(1) = (speed_ms(2) - speed_ms(1)) / (time(2) - time(1));
acceleration(end) = (speed_ms(end) - speed_ms(end-1)) / (time(end) - time(end-1));

% Add the calculated acceleration to the data_unit table
data_unit.acceleration = acceleration;

% Total Distance
total_distance_km = sum((speed_ms(1:end-1) .* diff(time))) / 1000;
fprintf('Total Distance: %.2f km\n', total_distance_km);

%% Power Model  
% Calculate pack power
pack_power = a * speed_ms + b * speed_ms.^2 + c * speed_ms.^3 + (1 + epsilon) * m_vehicle * speed_ms .* acceleration;
data_unit.pack_power = pack_power;

% Convert to cell power by dividing by the total number of cells (m * n)
cell_power = pack_power / (m_series * n_parallel);
data_unit.cell_power = cell_power;

%% Current Calculation

% -I^2 * R + OCV * I - P = 0 , P = I * V, V = OCV - I * R 
current = zeros(size(cell_power));

for i = 1:length(cell_power)
    P_cell = cell_power(i);
    discriminant = OCV_cell^2 - 4 * (-R_cell) * (-P_cell);

    if discriminant >= 0
        root1 = (-OCV_cell + sqrt(discriminant)) / (-2 * R_cell);
        root2 = (-OCV_cell - sqrt(discriminant)) / (-2 * R_cell);
        current(i) = min(root1, root2);
    else
        current(i) = NaN;
    end
end

data_unit.current = current;

%% C-rate and Scaled Current Calculation
C_rate = current / nominal_capacity_Ah;
scaled_current = C_rate * Scaling_nominal_capacity_Ah;

data_unit.C_rate = C_rate;
data_unit.scaled_current = scaled_current;

% Total charge and energy calculations
positive_current = current(current > 0);
positive_time = time(current > 0);

total_charge_As = trapz(positive_time, positive_current);
total_charge_Ah = total_charge_As / 3600;
total_energy_Wh = total_charge_Ah * OCV_cell * m_series;
total_used_Ah = trapz(time,scaled_current) / 3600; 

fprintf('Total Charge: %.2f Ah\n', total_charge_Ah);
fprintf('Total Energy: %.2f Wh\n', total_energy_Wh);
fprintf('Total used Cap: %.2f Ah\n', total_used_Ah);


%% Plot Speed, Acceleration, and Distance
figure;
subplot(3,1,1);
plot(time, speed_ms);
xlabel('Time (seconds)');
ylabel('Speed (m/s)');
title('Speed vs Time');
grid on;

subplot(3,1,2);
plot(time, acceleration);
xlabel('Time (seconds)');
ylabel('Acceleration (m/s^2)');
title('Acceleration vs Time');
grid on;

subplot(3,1,3);
plot(time, [0; cumsum(speed_ms(1:end-1) .* diff(time))] / 1000);
xlabel('Time (seconds)');
ylabel('Distance (km)');
title('Distance vs Time');
grid on;

%% Plot Power and Current
% Plot Power and Current on the same figure
figure;
subplot(2,1,1);
plot(time, pack_power);
xlabel('Time (seconds)');
ylabel('Power (W)');
title('Pack Power vs Time');
grid on;

subplot(2,1,2);
plot(time, current);
xlabel('Time (seconds)');
ylabel('Current (A)');
title('Cell Current vs Time');
grid on;

%% Plot C-rate and Scaled Current
figure;
subplot(2,1,1);
plot(time, C_rate);
xlabel('Time (seconds)');
ylabel('C-rate');
title('Cell C-rate vs Time');
grid on;

subplot(2,1,2);
plot(time, scaled_current);
xlabel('Time (seconds)');
ylabel('Current (A)');
title('Scaled Cell Current vs Time');
grid on;

%% Save Results to Excel
output_table = table(time, scaled_current);

if file_choice == 1
    output_file_path = 'udds_unit_time_scaled_current.xlsx';
else
    output_file_path = 'hwycol_unit_time_scaled_current.xlsx';
end

writetable(output_table, output_file_path);
fprintf('Excel file created successfully: %s\n', output_file_path);
