clc; clear; close all

%% UDDS convert to Power

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

% 변환된 속도를 새로운 변수에 저장
UDDS_unit.('speed(m/s)') = speed_ms;

% 시간 간격 계산 (가정: 균일한 시간 간격)
dt = diff(time);

% 주행거리 계산 (속도와 시간의 적분)
distance = sum(speed_ms(1:end-1) .* dt); % m 단위

% km로 변환
distance_km = distance / 1000;

% 결과 출력
fprintf('총 주행거리: %.2f km\n', distance_km);

% 속도와 가속도에 대해 플롯 그리기
acceleration = diff(speed_ms) ./ dt;
acceleration = [acceleration; NaN];
UDDS_unit.('acceleration(m/s^2)') = acceleration;

figure;
subplot(3,1,1);
plot(time, speed_ms);
xlabel('Time (seconds)');
ylabel('Speed (m/s)');
title('UDDS Speed');
grid on;

subplot(3,1,2);
plot(time, acceleration);
xlabel('Time (seconds)');
ylabel('Acceleration (m/s^2)');
title('UDDS Acceleration');
grid on;

subplot(3,1,3);
plot(time, [0; cumsum(speed_ms(1:end-1) .* dt)] / 1000);
xlabel('Time (seconds)');
ylabel('Distance (km)');
title('UDDS Distance');
grid on;


