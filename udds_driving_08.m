clc; clear; close all;

%% File Selection
disp('Select the file to analyze:');
disp('1. UDDS (uddscol.txt)');
disp('2. HWYCOL (hwycol.txt)');
file_choice = input('Enter the number of the file you want to process (1 or 2): ');

if file_choice == 1
    file_path = 'C:\Users\USER\Documents\GitHub\LFP_driving\uddscol.txt';
    disp('You have selected UDDS.');
elseif file_choice == 2
    file_path = 'C:\Users\USER\Documents\GitHub\LFP_driving\hwycol.txt';
    disp('You have selected HWYCOL.');
else
    error('Invalid selection. Please enter 1 or 2.');
end

% Physical constants
a = 34.98 * 4.44822; % lbf to Newton
b = 0.08650 * 4.44822 / 0.44704; % lbf/mph to N/(m/s)
c = 0.014800 * 4.44822 / 0.44704^2; % lbf/mph^2 to N/(m/s)^2
m = 1927.768; % vehicle mass in kg
epsilon = 1.05;

% Read the file
data_unit = readtable(file_path, 'Delimiter', '\t');
data_unit.Properties.VariableNames{1} = 'time';
data_unit.Properties.VariableNames{2} = 'speed(mph)';

% Separate the first column as time and the second as speed

time = data_unit{:, 1};
speed_mph = data_unit{:, 2};

% Convert speed from mph to m/s (1 mph = 0.44704 m/s)
speed_ms = speed_mph * 0.44704;

% Calculate acceleration using central difference method
acceleration = zeros(size(speed_ms));

% Central difference for interior points
for i = 2:length(time)-1
    acceleration(i) = (speed_ms(i+1) - speed_ms(i-1)) / (time(i+1) - time(i-1));
end

% Forward difference for the first point
acceleration(1) = (speed_ms(2) - speed_ms(1)) / (time(2) - time(1));

% Backward difference for the last point
acceleration(end) = (speed_ms(end) - speed_ms(end-1)) / (time(end) - time(end-1));

% Add the calculated acceleration to the data_unit table
data_unit.acceleration = acceleration;

% Calculate power
power = a * speed_ms + b * speed_ms.^2 + c * speed_ms.^3 + (1 + epsilon) * m * speed_ms .* acceleration;
data_unit.power = power;

% Output results
fprintf('Total Distance: %.2f km\n', sum((speed_ms(1:end-1) .* diff(time))) / 1000);

% Plot graphs
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

figure(2);
plot(time, power);
xlabel('Time (seconds)');
ylabel('Power (W)');
title('Pack Power');
grid on;

%% Calculate current

% Parameters
OCV = 3.2 * 106; % 342V
R = 0.0009 * 106; % 16Ah to 23Ah

% Initialize current 
current = zeros(size(data_unit.power));

% 2차 방정식 해 구하기
for i = 1:length(data_unit.power)
    P = data_unit.power(i);
    a = -R;  % 106 series Resistance
    b = OCV; % 106 series OCV
    c = -P; % pack power

    % 판별식 계산
    discriminant = b^2 - 4 * a * c;

    % 판별식 > 0 만 계산
    if discriminant >= 0
        % Compute the roots
        root1 = (-b + sqrt(discriminant)) / (2 * a);
        root2 = (-b - sqrt(discriminant)) / (2 * a);
        
        % 전류의 두 해 중, min 값 가져오기
        current(i) = min(root1, root2); % pack current = cell current
    else
        % 허근일 경우 NaN 할당
        current(i) = NaN;
    end
end

data_unit.current = current;

positive_current = current(current > 0);
positive_time = time(current > 0);

total_charge_As = trapz(positive_time, positive_current); 

total_charge_Ah = total_charge_As / 3600;

% Calculate the total energy in watt-hours (Wh)
total_energy_Wh = total_charge_Ah * OCV;

% Display the results
fprintf('Total Charge: %.2f Ah\n', total_charge_Ah);
fprintf('Total Energy: %.2f Wh\n', total_energy_Wh);

% Plot current over time
figure(3);
plot(time, current);
xlabel('Time (seconds)');
ylabel('Current (A)');
title('Pack Current');
grid on;

%% Calculate C-rate and Scaled Current

nominal_capacity_Ah = 161; % Capacity of one cell in Ah
scaling_factor_Ah = 16; % Scaling factor

% C-rate [h-1] = Pack current [A] / Cell capacity [Ah]
C_rate = current / nominal_capacity_Ah;

% scaled current [A] = C-rate [h-1] * Cell 1개 용량 [Ah]
scaled_current = C_rate * scaling_factor_Ah;

data_unit.C_rate = C_rate;
data_unit.scaled_current = scaled_current;

results = struct();
results.total_distance_km = sum((speed_ms(1:end-1) .* diff(time))) / 1000;
results.total_charge_Ah = total_charge_Ah;
results.total_energy_Wh = total_energy_Wh;
results.C_rate = C_rate;
results.scaled_current = scaled_current;

% Display results
fprintf('Total Distance: %.2f km\n', results.total_distance_km);
fprintf('Total Charge: %.2f Ah\n', results.total_charge_Ah);
fprintf('Total Energy: %.2f Wh\n', results.total_energy_Wh);

%% Plot

figure(4);
plot(time, C_rate);
xlabel('Time (seconds)');
ylabel('C-rate');
title('Cell C-rate');
grid on;


figure(5);
plot(time, scaled_current);
xlabel('Time (seconds)');
ylabel('Current (A)');
title('Scaled Cell Current');
grid on;

%% Excel

time = data_unit.time;
scaled_current = data_unit.scaled_current;

output_table = table(time, scaled_current);

if file_choice == 1
    output_file_path = 'udds_unit_time_scaled_current.xlsx';
else
    output_file_path = 'hwycol_unit_time_scaled_current.xlsx';
end

writetable(output_table, output_file_path);

fprintf('Excel file created successfully: %s\n', output_file_path);
