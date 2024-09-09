clc; clear; close all;

%% Load the 1C Formation Data
formation_data = load('C:\Users\USER\Documents\GitHub\LFP_driving\4ch_byd_1c_3cycle.mat');
formation_voltage = formation_data.data(12).V; 
formation_current = formation_data.data(12).I;
formation_time = formation_data.data(12).t / 3600; % Convert time to hours
formation_capacity_Ah = formation_data.data(12).Q; % Use the final Q value as the maximum capacity

% Calculate SOC for the formation data using the provided Q (capacity)
formation_Q = cumtrapz(formation_time * 3600, formation_current) / 3600; % Q in Ah
formation_SOC = 1 + formation_Q / formation_capacity_Ah; % SOC starting from 1 (100%)

%% Load the OCV Data
ocv_data = load('C:\Users\USER\Documents\GitHub\LFP_driving\LFP_k1_0_05C_25degC.mat');
ocv_voltage = ocv_data.data(4).V; 
ocv_current = ocv_data.data(4).I;
ocv_time = ocv_data.data(4).t / 3600; % Convert time to hours

%% Plot SOC vs Voltage for 1C Formation Data
figure;
yyaxis left;
plot(formation_SOC, formation_voltage, 'b', 'LineWidth', 1.5);
ylabel('Voltage (V)');
xlabel('SOC');
title('SOC vs Voltage');
grid on;

%% Overlay OCV Data on the Same Plot
hold on;
plot(ocv_data.data(4).SOC  , ocv_data.data(4).V, 'R', 'LineWidth', 1.5);
hold off;

