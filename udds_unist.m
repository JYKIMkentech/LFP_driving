clc; clear; close all

%% UDDS convert to Power

file_path = 'C:\Users\김준연\Documents\MATLAB\LFP_driving\uddscol.txt';

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

% 가속도 계산 (m/s^2)
acceleration = diff(speed_ms) ./ dt;

% 가속도 벡터 길이를 시간 벡터와 맞추기 위해 마지막 값을 NaN으로 채우기
acceleration = [acceleration; 0];

% 가속도를 테이블에 추가
UDDS_unit.('acceleration(m/s^2)') = acceleration;


% 속도와 가속도에 대해 플롯 그리기
figure;
subplot(2,1,1);
plot(time, speed_ms);
xlabel('Time (seconds)');
ylabel('Speed (m/s)');
title('UDDS Speed Over Time (in m/s)');
grid on;

subplot(2,1,2);
plot(time, acceleration);
xlabel('Time (seconds)');
ylabel('Acceleration (m/s^2)');
title('UDDS Acceleration Over Time');
grid on;


