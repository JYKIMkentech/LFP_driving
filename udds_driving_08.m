clc; clear; close all;

%% UDDS convert to Power

% Physical constants
a = 34.98 * 4.44822; % lbf to Newton
b = 0.08650 * 4.44822 / 0.44704; % lbf/mph to N/(m/s)
c = 0.014800 * 4.44822 / 0.44704^2; % lbf/mph^2 to N/(m/s)^2
m = 1927.768; % vehicle mass in kg
epsilon = 1.05;

file_path = 'C:\Users\deu04\OneDrive\바탕 화면\ECM_github\LFP_driving\uddscol.txt';

% Read the file
UDDS_unit = readtable(file_path, 'Delimiter', '\t');
UDDS_unit.Properties.VariableNames{1} = 'time';
UDDS_unit.Properties.VariableNames{2} = 'speed(mph)';

% Separate the first column as time and the second as speed
time = UDDS_unit{:, 1};
speed_mph = UDDS_unit{:, 2};

% Convert speed from mph to m/s (1 mph = 0.44704 m/s)
speed_ms = speed_mph * 0.44704;

% Calculate time intervals (assuming uniform time intervals)
dt = diff(time); % time intervals

% Calculate acceleration
acceleration = diff(speed_ms) ./ dt;
acceleration = [acceleration; 0]; % set last acceleration to 0
UDDS_unit.acceleration = acceleration;

% Calculate power
power = a * speed_ms + b * speed_ms.^2 + c * speed_ms.^3 + (1 + epsilon) * m * speed_ms .* acceleration;
UDDS_unit.power = power;

% Output results
fprintf('Total Distance: %.2f km\n', sum(speed_ms(1:end-1) .* dt) / 1000);

% Plot graphs
figure;
subplot(4,1,1);
plot(time, speed_ms);
xlabel('Time (seconds)');
ylabel('Speed (m/s)');
title('UDDS Speed');
grid on;

subplot(4,1,2);
plot(time, acceleration);
xlabel('Time (seconds)');
ylabel('Acceleration (m/s^2)');
title('UDDS Acceleration');
grid on;

subplot(4,1,3);
plot(time, [0; cumsum(speed_ms(1:end-1) .* dt)] / 1000);
xlabel('Time (seconds)');
ylabel('Distance (km)');
title('UDDS Distance');
grid on;

figure(2);
plot(time, power);
xlabel('Time (seconds)');
ylabel('Power (W)');
title('Power vs Time');
grid on;

%% Calculate current

% Parameters
OCV = 3.0; % Open-circuit voltage in volts
R = 0.009 * (16/23); % Resistance in ohms

% Initialize current array
current = zeros(size(UDDS_unit.power));

% Solve quadratic equation for each time step
for i = 1:length(UDDS_unit.power)
    P = UDDS_unit.power(i);
    a = -R;
    b = OCV;
    c = -P;

    % Calculate discriminant
    discriminant = b^2 - 4 * a * c;

    % Check if discriminant is non-negative
    if discriminant >= 0
        % Compute the roots
        root1 = (-b + sqrt(discriminant)) / (2 * a);
        root2 = (-b - sqrt(discriminant)) / (2 * a);
        
        % Assuming we want the positive root since current cannot be negative
        current(i) = max(root1, root2);
    else
        % If no real solution exists, set current to NaN
        current(i) = NaN;
    end
end

% Add current to the UDDS_unit table
UDDS_unit.current = current;

% Integrate the positive current over time to get the total charge in Ampere-seconds
positive_current = current(current > 0);
positive_time = time(current > 0);

total_charge_As = trapz(positive_time, positive_current); 

% Convert charge from Ampere-seconds to Ampere-hours (1 Ah = 3600 As)
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
title('Current vs Time');
grid on;
