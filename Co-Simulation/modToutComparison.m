%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing the comparison of the lumped modelica storage model
% G3P3 simulation to the modified MATLAB storage model simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all

%% load data
load('G3P3_3day.mat', 'x');
xMod = load('G3P3_3day_Tout_mod2.mat'); xMod = xMod.x;
xCte = load('G3P3_3day_Tout_cte.mat'); xCte = xCte.x;
load('G3P3_3day_HX.mat', 'HX_params');
HX_params_mod = load('G3P3_3day_Tout_mod2_HX.mat');
HX_params_mod = HX_params_mod.HX_params;
HX_params_cte = load('G3P3_3day_T_out_cte_HX.mat');
HX_params_cte = HX_params_cte.HX_params;
jcn_params = load('G3P3_3day_jcn.mat');
jcn_params = jcn_params.Junction_params;
jcn_params_mod = load('G3P3_3day_Tout_mod2_jcn.mat');
jcn_params_mod = jcn_params_mod.Junction_params;
%% hot storage bin outlet temperatures
figure('Units', 'normalized', 'Position', [0 0 0.4 0.3], 'Visible', 'on');
plot(x.hour(4:end), x.hotBin_T_s_out(4:end) - 273.15, '--k'); hold on;
plot(xMod.hour(4:end), xMod.hotBin_T_s_out(4:end) - 273.15, '-r');
plot(xCte.hour(4:end), xCte.hotBin_T_s_out(4:end), '.-b');
ylabel('$T_{out}$ ($^\circ$C)', 'interpreter', 'latex', 'FontSize', 14);
xlabel('hour from start-up', 'interpreter', 'latex', 'FontSize', 14);
xlim([0, 72]);
ylim([350, 800]);
set(gca, 'TickLabelInterpreter', 'latex')
set(gcf, 'Color', [1 1 1])
legend('Original Outlet Temp', 'Modified Outlet Temp', 'Constant Outlet Temp', ...
                    'interpreter', 'latex', 'FontSize', 14, 'Location', 'SouthEast');


%% heater input
figure('Units', 'normalized', 'Position', [0 0 0.4 0.3], 'Visible', 'on');
plot(HX_params.hour(2:end), HX_params.Qheater(2:end), '--k'); hold on;
plot(HX_params_mod.hour(2:end), HX_params_mod.Qheater(2:end), '-r');
plot(HX_params_cte.hour(2:end), HX_params_cte.Qheater(2:end), '.-b');
ylabel('$\dot{Q}_{heater}$ (W)', 'interpreter', 'latex', 'FontSize', 14);
xlabel('hour from start-up', 'interpreter', 'latex', 'FontSize', 14);
xlim([0, 72]);
% ylim([350, 800]);
set(gca, 'TickLabelInterpreter', 'latex')
set(gcf, 'Color', [1 1 1])
legend('Original Outlet Temp', 'Modified Outlet Temp', 'Constant Outlet Temp', ...
                    'interpreter', 'latex', 'FontSize', 14, 'Location', 'SouthEast');


%% HX CO2 output
% figure('Units', 'normalized', 'Position', [0 0 0.4 0.3], 'Visible', 'on');
% plot(HX_params.hour(2:end), HX_params.m_dot_CO2_out(2:end).*HX_params.hCO2_20(2:end), '--k'); hold on;
% plot(HX_params_mod.hour(2:end), HX_params_mod.m_dot_CO2_out(2:end).*HX_params_mod.hCO2_20(2:end), '-r');
% plot(HX_params_cte.hour(2:end), HX_params_cte.m_dot_CO2_out(2:end).*HX_params_cte.hCO2_20(2:end), '.-b');
% ylabel('$\dot{Q}_{heater}$ (W)', 'interpreter', 'latex', 'FontSize', 14);
% xlabel('hour from start-up', 'interpreter', 'latex', 'FontSize', 14);
% xlim([0, 72]);
% % ylim([350, 800]);
% set(gca, 'TickLabelInterpreter', 'latex')
% set(gcf, 'Color', [1 1 1])
% legend('Original Outlet Temp', 'Modified Outlet Temp', 'Constant Outlet Temp', ...
%                     'interpreter', 'latex', 'FontSize', 14, 'Location', 'SouthEast');
                
%% junction mass flow rates
figure('Units', 'normalized', 'Position', [0 0 0.4 0.3], 'Visible', 'on');
plot(jcn_params.hour(4:end), jcn_params.hotBinBypass_m_dot_s_out(4:end), '--k'); hold on;
plot(jcn_params_mod.hour(4:end), jcn_params_mod.hotBinBypass_m_dot_s_out(4:end), '-r');
ylabel('$\dot{m}$ (kg/s)', 'interpreter', 'latex', 'FontSize', 14);
xlabel('hour from start-up', 'interpreter', 'latex', 'FontSize', 14);
xlim([0, 72]);
% ylim([0, 0.2]);
set(gca, 'TickLabelInterpreter', 'latex')
set(gcf, 'Color', [1 1 1])
legend('Original Outlet Temp', 'Modified Outlet Temp', ...
                    'interpreter', 'latex', 'FontSize', 14, 'Location', 'SouthEast');
