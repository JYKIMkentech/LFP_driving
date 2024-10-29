clc; clear; close all;

% 데이터 파일 불러오기
data = readtable('G:\공유 드라이브\Battery Software Lab\Driving cycle\16Ah_BYD\Processed\BSL_CITY1_time_scaled_current.xlsx');

% 변수 설정
time = data.time;
current = data.scaled_current;

% 전하량(Q) 계산
Q = cumtrapz(time, current); % 적분하여 누적 전하량 계산
total_capacity = trapz(time, current); % 전체 전하량 (최종 용량)

% 그래프 생성
figure;
yyaxis left
plot(time, current, '-');
xlabel('Time (s)');
ylabel('Current (A)');
title('Current and Charge (Q) vs Time');
grid on;

yyaxis right
plot(time, Q, '-r');
ylabel('Charge (Q) (As)');

% 전체 용량 표시
text(time(end) * 0.7, max(Q) * 0.9, ['Total Capacity: ', num2str(total_capacity), ' As'], 'FontSize', 12, 'Color', 'r');
