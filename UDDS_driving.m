clc; clear; close all;

%% UDDS convert to Power

% 물리 상수
a = 34.98 * 4.44822; % lbf to Newton
b = 0.08650 * 4.44822 / 0.44704; % lbf/mph to N/(m/s)
c = 0.014800 * 4.44822 / 0.44704^2; % lbf/mph^2 to N/(m/s)^2
m = 1927.768; % 차량 질량 kg
epsilon = 1.05;

file_path = 'C:\Users\USER\Documents\GitHub\LFP_driving\uddscol.txt';

% 파일을 읽어오기
UDDS_unit = readtable(file_path, 'Delimiter', '\t');
UDDS_unit.Properties.VariableNames{1} = 'time';
UDDS_unit.Properties.VariableNames{2} = 'speed(mph)';

% 첫 번째 열을 시간, 두 번째 열을 속도로 분리
time = UDDS_unit{:, 1};
speed_mph = UDDS_unit{:, 2};

% 속도를 mph에서 m/s로 변환 (1 mph = 0.44704 m/s)
speed_ms = speed_mph * 0.44704;

% 시간 간격 계산 (가정: 균일한 시간 간격)
dt = diff(time); % 시간 간격

% 가속도 계산
acceleration = diff(speed_ms) ./ dt;
acceleration = [acceleration; 0]; % 마지막 가속도는 0으로 처리
UDDS_unit.acceleration = acceleration;

% 파워 계산
power = a * speed_ms + b * speed_ms.^2 + c * speed_ms.^3 + (1 + epsilon) * m * speed_ms .* acceleration;
UDDS_unit.power = power;

% power 변환

power_scaled = (power * 23) /(106 * 161);
UDDS_unit.power_scaled = power_scaled;

% 결과 출력
fprintf('총 주행거리: %.2f km\n', sum(speed_ms(1:end-1) .* dt) / 1000);

% 그래프 표시
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

figure(3);
plot(time, power_scaled);
xlabel('Time (seconds)');
ylabel('Power scaled (W)');
title('Power scaled vs Time');
grid on;

%% 전류 구하기

% Parameters
OCV = 3.0; % Open-circuit voltage in volts
R = 0.009; % Resistance in ohms

% Initialize current array
current = zeros(size(UDDS_unit.power_scaled));

% Solve quadratic equation for each time step
for i = 1:length(UDDS_unit.power_scaled)
    P = UDDS_unit.power_scaled(i);
    a = R;
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

% Plot current over time
figure(4);
plot(time, current);
xlabel('Time (seconds)');
ylabel('Current (A)');
title('Current vs Time');
grid on;














