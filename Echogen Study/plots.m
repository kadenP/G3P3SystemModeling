%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots used for Echogen study
% Kaden Plewe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all

%% load simulation data
t = load('time.mat'); t = t.t;
Ts_in = load('Ts_in.mat'); Ts_in = Ts_in.Ts_in;
Ts_out_1 = load('Ts_out_case1.mat'); Ts_out_1 = Ts_out_1.Ts_out;
Ts_out_2 = load('Ts_out_case2.mat'); Ts_out_2 = Ts_out_2.Ts_out;
Ts_out_3 = load('Ts_out_case3.mat'); Ts_out_3 = Ts_out_3.Ts_out;
Ts_out_4 = load('Ts_out_case4.mat'); Ts_out_4 = Ts_out_4.Ts_out;
qLoss_1 = load('qLoss_case1.mat'); qLoss_1 = qLoss_1.qLoss;
qLoss_2 = load('qLoss_case2.mat'); qLoss_2 = qLoss_2.qLoss;
qLoss_3 = load('qLoss_case3.mat'); qLoss_3 = qLoss_3.qLoss;
qLoss_4 = load('qLoss_case1.mat'); qLoss_4 = qLoss_4.qLoss;

%% plot temperatures

f = figure('Units', 'normalized', 'Position', [0 0 0.2 0.2], 'Visible', 'on');
plot(t/3600, Ts_in, '-k', 'LineWidth', 1.5); hold on;
plot(t/3600, Ts_out_1, '-', 'Color', [0.8500 0.7*0.3250 0.0980], 'LineWidth', 1.5);
plot(t/3600, Ts_out_2, '-', 'Color', [0.8500 0.8*0.3250 0.0980], 'LineWidth', 1.5);
plot(t/3600, Ts_out_3, '-', 'Color', [0.8500 0.9*0.3250 0.0980], 'LineWidth', 1.5);
plot(t/3600, Ts_out_4, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);

ylabel('$T$ ($^\circ C$)', 'interpreter', 'latex', 'FontSize', 14);
xlabel('$t$ ($h$)');

