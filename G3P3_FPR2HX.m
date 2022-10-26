%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G3P3 System Simulation
%
% Combines all of the component models to simulate G3P3 operation. This is
% a precursor to the eventual Simulink model that will combine the same
% models in the Simulink framework.
%
% Kaden Plewe
% 09/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all

%% setup cycle timing

% setup cycles (time units in seconds)
chargeDurration = 6*3600;
holdDurration = 10*3600;
dischargeDurration = 8*3600;
numCycles = 7;

dt1 = 10;
t1 = 0:dt1:(chargeDurration + holdDurration + dischargeDurration)*numCycles;
% t1 = 0:dt1:3600*2;
t1 = t1';

dt2 = 1200;
t2 = 0:dt2:(chargeDurration + holdDurration + dischargeDurration)*numCycles;


%% Define all system components (note that default model parameters are set to match the Modelica G3P3 model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heliostat Field
% Inputs: DNI
% Outputs: Qsolar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HF_ = HF();

Qsolar = zeros(length(t1), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Falling Particle Reciever
% Inputs: Ts_in, Tinf, mdot_s_in, Qsolar, t
% Outputs: Ts_out, mdot_s_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FPR_ = FPR();
FPR_.Ts0 = 600;
FPR_.tauLag = 1500;
FPR_.tauLead = 1e-6;
FPR_.dtLin = 10;

Ts_in_FPR = FPR_.Ts0*ones(length(t1), 1);
mdot_s_in_FPR = zeros(length(t1), 1);
Ts_out_FPR = FPR_.Ts0*ones(length(t1), 1);
Ts_out_lin_FPR = FPR_.Ts0*ones(length(t1), 1);
mdot_s_out_FPR = zeros(length(t1), 1);
Ts_FPR = cell(length(t1), 1); Ts_FPR{1} = FPR_.Ts0*ones(size(t1));
x_FPR = cell(length(t1), 1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hot Bin
% Inputs: Tin, mdot, t
% Outputs: Ts_out, T_bulk, Estored, ztop, ms, mdot_s_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hotTES = TES();
hotTES.T0 = 775;
hotTES.a = 0.065;
hotTES.a0 = 0.065;

Ts_in_hotTES = zeros(length(t2), 1);
mdot_hotTES = zeros(length(t2), 1);
Ts_out_hotTES = zeros(length(t2), 1);
Ts_bulk_hotTES = zeros(length(t2), 1);
Estored_hotTES = zeros(length(t2), 1);
ztop_hotTES = zeros(length(t2), 1);
ms_hotTES = zeros(length(t2), 1);
mdot_s_out_hotTES = zeros(length(t2), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particle-to-sCO2 Heat Exchanger
% Inputs: Ts_in, Tco2_in, mdot_s_in, mdot_CO2_in, t
% Outputs: Ts_out, Tco2_out, mdot_s_out, mdot_CO2_out, Ts, 
%          Tco2, Tm, Q_CO2, x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HX_ = HX();
HX_.n = 30;
HX_.Ts0 = 25;
HX_.Tco20 = 25;
HX_.h_s = 400;                   % (W/m2K) solid-wall heat transfer coefficient 
HX_.h_CO2 = 3000;                % (W/m2K) CO2 heat transfer coefficient 

Ts_in_HX = zeros(length(t1), 1); Ts_in_HX(1) = HX_.Ts0;
Tco2_in_HX = zeros(length(t1), 1); Tco2_in_HX(:) = 550;
mdot_s_in_HX = zeros(length(t1), 1); mdot_s_in_HX(:) = 5;
mdot_CO2_in_HX = zeros(length(t1), 1); mdot_CO2_in_HX(:) = 5;
Ts_out_HX = zeros(size(t1)); Ts_out(1) = HX_.Ts0;
Tco2_out_HX = zeros(size(t1)); Tco2_out(1) = HX_.Tco20;
Ts_out_lin_HX = Ts_out_HX;
Tco2_out_lin_HX = Tco2_out_HX;
mdot_s_out_HX = mdot_s_in_HX;
mdot_CO2_out_HX = mdot_CO2_in_HX;
Ts_HX = cell(size(t1)); Ts{1} = HX_.Ts0*ones(HX_.n, 1);
Tco2_HX = cell(size(t1)); Tco2{1} = HX_.Tco20*ones(HX_.n, 1);
Tm_HX = cell(size(t1)); Tm{1} = HX_.Tm0*ones(HX_.n, 1);
Q_CO2 = zeros(size(t1));
Q_s = zeros(size(t1));
x_HX = cell(length(t1), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Free Fall Ducts
% Inputs: Ts_in, mdot_s_in, Tinf, t
% Outputs: Ts_out, mdot_s_out, Ts, Tm, x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RecieverDownComer = FFD();
RecieverDownComer.L = 10;

Ts_in_RecieverDownComer = RecieverDownComer.Ts0*ones(length(t1), 1);
mdot_s_in_RecieverDownComer = zeros(length(t1), 1);
Ts_out_RecieverDownComer = RecieverDownComer.Ts0*ones(length(t1), 1);
mdot_s_out_RecieverDownComer = zeros(length(t1), 1);
Ts_RecieverDownComer = cell(length(t1), 1);
Tm_RecieverDownComer = cell(length(t1), 1);
x_RecieverDownComer = cell(length(t1), 1);

HotBinInlet = FFD();
HotBinInlet.L = 10;

Ts_in_HotBinInlet = HotBinInlet.Ts0*ones(length(t1), 1);
mdot_s_in_HotBinInlet = zeros(length(t1), 1);
Ts_out_HotBinInlet = HotBinInlet.Ts0*ones(length(t1), 1);
mdot_s_out_HotBinInlet = zeros(length(t1), 1);
Ts_HotBinInlet = cell(length(t1), 1);
Tm_HotBinInlet = cell(length(t1), 1);
x_HotBinInlet = cell(length(t1), 1);

HotBinDischarge = FFD();
HotBinDischarge.L = 10;

Ts_in_HotBinDischarge = HotBinDischarge.Ts0*ones(length(t1), 1);
mdot_s_in_HotBinDischarge = zeros(length(t1), 1);
Ts_out_HotBinDischarge = HotBinDischarge.Ts0*ones(length(t1), 1);
mdot_s_out_HotBinDischarge = zeros(length(t1), 1);
Ts_HotBinDischarge = cell(length(t1), 1);
Tm_HotBinDischarge = cell(length(t1), 1);
x_HotBinDischarge = cell(length(t1), 1);

HeatExchangerDownComer = FFD();
HeatExchangerDownComer.L = 10;

Ts_in_HeatExchangerDownComer = HeatExchangerDownComer.Ts0*ones(length(t1), 1);
mdot_s_in_HeatExchangerDownComer = zeros(length(t1), 1);
Ts_out_HeatExchangerDownComer = HeatExchangerDownComer.Ts0*ones(length(t1), 1);
mdot_s_out_HeatExchangerDownComer = zeros(length(t1), 1);
Ts_HeatExchangerDownComer = cell(length(t1), 1);
Tm_HeatExchangerDownComer = cell(length(t1), 1);
x_HeatExchangerDownComer = cell(length(t1), 1);

%% System data table

% High resolution component data
sysData1 = table(Qsolar, Ts_in_FPR, mdot_s_in_FPR, Ts_out_FPR, mdot_s_out_FPR, ...
    Ts_FPR, x_FPR, Ts_in_RecieverDownComer, mdot_s_in_RecieverDownComer, ...
    Ts_out_RecieverDownComer, mdot_s_out_RecieverDownComer, ...
    Ts_RecieverDownComer, Tm_RecieverDownComer, x_RecieverDownComer, ...
    Ts_in_HotBinInlet, mdot_s_in_HotBinInlet, Ts_out_HotBinInlet, ...
    mdot_s_out_HotBinInlet, Ts_HotBinInlet, Tm_HotBinInlet, ...
    x_HotBinInlet, Ts_in_HeatExchangerDownComer, mdot_s_in_HeatExchangerDownComer, ...
    Ts_out_HeatExchangerDownComer, mdot_s_out_HeatExchangerDownComer, ...
    Ts_HeatExchangerDownComer, Tm_HeatExchangerDownComer, ...
    x_HeatExchangerDownComer);

sysData3 = table(Ts_in_HX, Tco2_in_HX, mdot_s_in_HX, mdot_CO2_in_HX, Ts_out_HX, Tco2_out_HX, ...
    mdot_s_out_HX, mdot_CO2_out_HX, Ts_HX, Tco2_HX, Tm_HX, Q_CO2, ...
    Q_s, x_HX, Ts_out_lin_HX, Tco2_out_lin_HX, Ts_in_HotBinDischarge, mdot_s_in_HotBinDischarge, ...
    Ts_out_HotBinDischarge, mdot_s_out_HotBinDischarge, ...
    Ts_HotBinDischarge, Tm_HotBinDischarge, x_HotBinDischarge);

% Low resolution component data
sysData2 = table(Ts_in_hotTES, mdot_hotTES, Ts_out_hotTES, Ts_bulk_hotTES, Estored_hotTES, ...
    ztop_hotTES, ms_hotTES, mdot_s_out_hotTES);


%% load TMY3 data and set simulation dates
weather = readtable('ABQ_Weather_Lookup.xlsx');

% Winter Week (Dec 21 - Dec 28)
t_h_start = 8497 + 8; [~, ti_start] = min(abs(weather.time - t_h_start));
t_h_end = 8665 + 8; [~, ti_end] = min(abs(weather.time - t_h_end));


% Spring Week (March 21 - March 28)
% t_h_start = 2641 + 8; [~, ti_start] = min(abs(weather.time - t_h_start));
% t_h_end = 2809 + 8; [~, ti_end] = min(abs(weather.time - t_h_end));


% Summer Week (June 21 - June 28)
% t_h_start = 4105 + 8; [~, ti_start] = min(abs(weather.time - t_h_start));
% t_h_end = 4273 + 8; [~, ti_end] = min(abs(weather.time - t_h_end));


% Fall Week (September 23 - September 30)
% t_h_start = 6361 + 8; [~, ti_start] = min(abs(weather.time - t_h_start));
% t_h_end = 6529 + 8; [~, ti_end] = min(abs(weather.time - t_h_end));


%% Simulate over prescribed time frame

% additional simulation data

hour_day1 = zeros(size(t1));
hour_day2 = zeros(size(t2));
DNI = zeros(size(t1));
Tinf1 = zeros(size(t1));
Tinf2 = zeros(size(t2));

% time for system jump start
iStart1 = 1;
iStart2 = 2;
% [~, iStart2] = min(abs(9*3600 - t2));

% initialize system parameters


% profile on;

%% simulate reciever with high resolution

% for i = iStart1:length(t1)
%     
%     % load current weather conditions
%     hour_day1(i) = mod(t_h_start + t1(i)/3600, 24);
%     DNI(i) = interp1(weather.time, weather.DNI, t_h_start + t1(i)/3600, 'makima');
%     Tinf1(i) = interp1(weather.time, weather.Tinf, t_h_start + t1(i)/3600, 'makima');
%     
%     % set particle-to-sCO2 inlet temperature and flow rate
%     sysData1.Tco2_in(i) = 565;
%     sysData1.mdot_CO2_in_HX(i) = 5;
%     sysData1.mdot_CO2_out_HX(i) = sysData1.mdot_CO2_in_HX(i);
%     
%     % reciever collection limit
% %     sysData1.Qsolar(i) = step(HF_, DNI(i));
%     
%     % intermediate storage diverter valve control
%     sysData1.y1_intermediateStorageDiverter(i) = 1; % no flow to cold bin
%     
%     % simulate operation based on hour of the day
%     sysData1.Ts_in_FPR(i) = 600;
%     if hour_day1(i) >= 1 && hour_day1(i) < 9
%         sysData1.Qsolar(i) = 0;
%         sysData1.mdot_s_in_FPR(i) = 9;               
%         sysData1.mdot_hotTES(i) = -5;
%         sysData1.mdot_s_in_hotBinDischarge(i) = 5;        
%     elseif hour_day1(i) >= 9 && hour_day1(i) < 15
%         sysData1.Qsolar(i) = step(HF_, DNI(i));
%         sysData1.mdot_s_in_FPR(i) = 9;
%         sysData1.mdot_hotTES(i) = 9;
%         sysData1.mdot_s_in_hotBinDischarge(i) = 0;        
%     else
%         sysData1.Qsolar(i) = 0;
%         sysData1.mdot_s_in_FPR(i) = 9;
%         sysData1.mdot_hotTES(i) = 0;
%         sysData1.mdot_s_in_hotBinDischarge(i) = 0;
%     end
%     
%     % simulate system starting at the falling particle reciever
% %     [sysData1.Ts_out_FPR(i), sysData1.mdot_s_out_FPR(i), ...
% %         sysData1.Ts_FPR{i}, sysData1.x_FPR{1}] = ...
% %         step(FPR_, sysData1.Ts_in_FPR(i), Tinf1(i), sysData1.mdot_s_in_FPR(i), ...
% %          sysData1.Qsolar(i), t1(i) - t1(iStart1));
%     [sysData1.Ts_out_FPR(i), Ts_out_lin_FPR(i), sysData1.mdot_s_out_FPR(i)] = ...
%         step(FPR_, sysData1.Ts_in_FPR(i), Tinf1(i), ...
%         sysData1.Qsolar(i), t1(i) - t1(iStart1), 775);
%      
%     sysData1.Ts_in_RecieverDownComer(i) = sysData1.Ts_out_FPR(i);
%     sysData1.mdot_s_in_RecieverDownComer(i) = sysData1.mdot_s_out_FPR(i);
%     
%     [sysData1.Ts_out_RecieverDownComer(i), sysData1.mdot_s_out_RecieverDownComer(i), ...
%         sysData1.Ts_RecieverDownComer{i}, sysData1.Tm_RecieverDownComer{i}, ...
%         sysData1.x_RecieverDownComer{i}] = step(RecieverDownComer, ...
%         sysData1.Ts_in_RecieverDownComer(i), ...
%         sysData1.mdot_s_in_RecieverDownComer(i), Tinf1(i), t1(i) - t1(iStart1));
%     fprintf('Ts_out_FPR = %1.2f °C\n', sysData1.Ts_out_FPR(i));
% end
% 
% % save data table
% save('sysData1.mat', 'sysData1');


%% simulate storage with low resolution

% load data from receiver operation
% load('sysData1_Winter.mat', 'sysData1')
% load('sysData1_Summer.mat', 'sysData1')
% load('sysData1_Spring.mat', 'sysData1')
% load('sysData1_Fall.mat', 'sysData1')



% for i = iStart2:length(t2)
%     hour_day2(i) = mod(t_h_start + t2(i)/3600, 24);
%     Tinf2(i) = interp1(weather.time, weather.Tinf, t_h_start + t2(i)/3600, 'makima');
%     
%     if hour_day2(i) >= 1 && hour_day2(i) < 9             
%         sysData2.mdot_hotTES(i) = -5;        
%     elseif hour_day2(i) >= 9 && hour_day2(i) < 15
%         sysData2.mdot_hotTES(i) = 9;       
%     else
%         sysData2.mdot_s_in_hotBinDischarge(i) = 0;
%     end
%     
%     sysData2.Ts_in_hotTES(i) = interp1(t_h_start + t1/3600, sysData1.Ts_out_RecieverDownComer, ...
%         t_h_start + t2(i)/3600, 'makima');
% %     sysData2.Ts_in_hotTES(i) = 775;
%     
%     [sysData2.Ts_out_hotTES(i), sysData2.Ts_bulk_hotTES(i), sysData2.Estored_hotTES(i), ...
%         sysData2.ztop_hotTES(i), sysData2.ms_hotTES(i), sysData2.mdot_s_out_hotTES(i)] = ...
%         step(hotTES, sysData2.Ts_in_hotTES(i), Tinf2(i), sysData2.mdot_hotTES(i), t2(i));    
% end
% 
% 
% % save data table
% save('sysData2.mat', 'sysData2');

%% simulate heat exchanger with high resolution

% load data from storage operation
load('sysData2_Winter.mat', 'sysData2')
% load('sysData2_Summer.mat', 'sysData2')
% load('sysData2_Spring.mat', 'sysData2')
% load('sysData2_Fall.mat', 'sysData2')

for i = 1:length(t1)
    hour_day1(i) = mod(t_h_start + t1(i)/3600, 24);
    Tinf1(i) = interp1(weather.time, weather.Tinf, t_h_start + t1(i)/3600, 'makima');
    
    if hour_day1(i) >= 1 && hour_day1(i) < 9             
        sysData3.mdot_s_in_HotBinDischarge(i) = 5;
        sysData3.Ts_in_HotBinDischarge(i) = ...
                interp1(t2, sysData2.Ts_out_hotTES, t1(i), 'makima');               
    elseif hour_day1(i) >= 9 && hour_day1(i) < 15
        sysData3.mdot_s_in_HotBinDischarge(i) = 5;
        sysData3.Ts_in_HotBinDischarge(i) = 600;      
    else
        sysData3.mdot_s_in_HotBinDischarge(i) = 5;
        sysData3.Ts_in_HotBinDischarge(i) = 600; 
    end
    
    [sysData3.Ts_out_HotBinDischarge(i), sysData3.mdot_s_out_HotBinDischarge(i), ...
        sysData3.Ts_HotBinDischarge{i}, sysData3.Tm_HotBinDischarge{i}, ...
        sysData3.x_HotBinDischarge{i}] = step(HotBinDischarge, ...
        sysData3.Ts_in_HotBinDischarge(i), ...
        sysData3.mdot_s_in_HotBinDischarge(i), Tinf1(i), t1(i));
    
    sysData3.Ts_in_HX(i) = sysData3.Ts_out_HotBinDischarge(i);
    
    [sysData3.Ts_out_HX(i), sysData3.Tco2_out_HX(i), ...
        sysData3.mdot_s_out_HX(i), sysData3.mdot_CO2_out_HX(i), ...
        sysData3.Ts_HX{i}, sysData3.Tco2_HX{i}, sysData3.Tm_HX{i}, ...
        sysData3.Q_CO2(i), sysData3.Q_s(i), sysData3.x_HX{i}, ...
        sysData3.Ts_out_lin_HX(i), sysData3.Tco2_out_lin_HX(i)] ...
        = step(HX_, sysData3.Ts_in_HX(i), sysData3.Tco2_in_HX(i), sysData3.mdot_s_in_HX(i), ...
                sysData3.mdot_CO2_in_HX(i), t1(i));
    fprintf('Tco2_out_HX = %1.2f °C\n', sysData3.Tco2_out_HX(i));
    
end

% save data table
save('sysData3.mat', 'sysData3');

% profsave;



