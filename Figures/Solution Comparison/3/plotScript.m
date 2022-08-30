
clc, clear, close all

%% load data
baselineFile = 'baselineVariables.csv';
modifiedFile = 'modifiedTsoutVariables.csv';
baselineTable = readtable(baselineFile);
modifiedTable = readtable(modifiedFile);


%% plot a comparison of baseline and modified results

% hot bin outlet temperature
f1 = figure('Units', 'normalized', 'Position', [0 0 0.4 0.3], ...
    'Visible', 'on'); 
plot(baselineTable.time/3600, baselineTable.HotBinDischarge_ParticleInlet_T - 273.15, '-k', 'DisplayName', 'baseline'); 
hold on;
plot(modifiedTable.time/3600, modifiedTable.HotBinDischarge_ParticleInlet_T - 273.15, '-r', 'DisplayName', 'modified storage Tout');
legend('interpreter', 'latex', 'FontSize', 12);
ylim([350, 800]);
xlim([0, 72]);
ylabel('Hot Bin Particle Outlet Temp ($^\circ$C)', 'interpreter', 'latex', ...
        'FontSize', 14);
xlabel('$t$ (h)', 'interpreter', 'latex', ...
        'FontSize', 14);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
set(gcf, 'Color', [1 1 1])

% particle heat exchanger particle inlet temperature
f2 = figure('Units', 'normalized', 'Position', [0 0 0.4 0.3], ...
    'Visible', 'on'); 
plot(baselineTable.time/3600, baselineTable.ParticleHeatExchanger_ParticleInlet_T - 273.15, '-k', 'DisplayName', 'baseline'); 
hold on;
plot(modifiedTable.time/3600, modifiedTable.ParticleHeatExchanger_ParticleInlet_T - 273.15, '-r', 'DisplayName', 'modified storage Tout');
legend('interpreter', 'latex', 'FontSize', 12);
ylim([350, 800]);
xlim([0, 72]);
ylabel('Particle HX Inlet Temp ($^\circ$C)', 'interpreter', 'latex', ...
        'FontSize', 14);
xlabel('$t$ (h)', 'interpreter', 'latex', ...
        'FontSize', 14);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
set(gcf, 'Color', [1 1 1])

% particle heater heat input
f3 = figure('Units', 'normalized', 'Position', [0 0 0.4 0.3], ...
    'Visible', 'on'); 
plot(baselineTable.time/3600, baselineTable.ParticleHeater_Q_in, '-k', 'DisplayName', 'baseline'); 
hold on;
plot(modifiedTable.time/3600, modifiedTable.ParticleHeater_Q_in, '-r', 'DisplayName', 'modified storage Tout');
legend('interpreter', 'latex', 'FontSize', 12);
% ylim([350, 800]);
xlim([0, 72]);
ylabel('Particle Heater Input (W)', 'interpreter', 'latex', ...
        'FontSize', 14);
xlabel('$t$ (h)', 'interpreter', 'latex', ...
        'FontSize', 14);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
set(gcf, 'Color', [1 1 1])

% particle heat exchanger CO2 outlet temperature
f4 = figure('Units', 'normalized', 'Position', [0 0 0.4 0.3], ...
    'Visible', 'on'); 
plot(baselineTable.time/3600, baselineTable.CO2Outlet_T - 273.15, '-k', 'DisplayName', 'baseline'); 
hold on;
plot(modifiedTable.time/3600, modifiedTable.CO2Outlet_T - 273.15, '-r', 'DisplayName', 'modified storage Tout');
legend('interpreter', 'latex', 'FontSize', 12);
ylim([500, 800]);
xlim([0, 72]);
ylabel('CO2 HX Outlet Temp ($^\circ$C)', 'interpreter', 'latex', ...
        'FontSize', 14);
xlabel('$t$ (h)', 'interpreter', 'latex', ...
        'FontSize', 14);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
set(gcf, 'Color', [1 1 1])