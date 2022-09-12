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

%% Define all system components (note that default model parameters are set to match the Modelica G3P3 model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heliostat Field
% Inputs: DNI
% Outputs: Qsolar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HF_ = HF();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass Flow Hopper
% Inputs: Ts_in, Tinf, mdot_s_in, Qsolar, t
% Outputs: Ts_out, mdot_s_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HX_MFH = MFH();
Reciever_MFH = MFH();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Falling Particle Reciever
% Inputs: Ts_in, Tinf, mdot_s_in, Qsolar, t
% Outputs: Ts_out, mdot_s_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FPR_ = FPR();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hot and Cold Storage Bins
% Inputs: Tin, mdot, t
% Outputs: Ts_out, T_bulk, Estored, ztop, mdot_s_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hotTES = TES();
hotTES.T0 = 800;
coldTES = TES();
coldTES.T0 = 800;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intermediate Lumped Storage Bin
% Inputs: Tin, Tinf, mdot_s_in, mdot_s_out, t
% Outputs: Ts_out, qloss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intTES = LSB();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particle-to-sCO2 Heat Exchanger
% Inputs: Ts_in, Tco2_in, mdot_s_in, mdot_CO2_in, t
% Outputs: Ts_out, Tco2_out, mdot_s_out, mdot_CO2_out, Ts, 
%          Tco2, Tm, Q_CO2, x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HX_ = HX();
HX_.n = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bucket Elevator
% Inputs: Ts_in, mdot_s_in, Tinf, t
% Outputs: Ts_out, mdot_s_out, Ts, Tm, x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BE_ = BE();
BE_.n = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solid Splitters
% Inputs: Ts_in, mdot_in
% Outputs: Ts_out1, Ts_out2, mdot_out1, mdot_out2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HotBinDiverter = SS();
RecieverDiverter = SS();
IntermediateStorageDiverter = SS();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solid Junctions
% Inputs: Ts_in1, Ts_in2, mdot_in1, mdot_in2
% Outputs: Ts_out, mdot_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HeatExchangerJunction = SJ();
ColdBinJunction = SJ();
BucketElevatorJunction = SJ();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Free Fall Ducts
% Inputs: Ts_in, mdot_s_in, Tinf, t
% Outputs: Ts_out, mdot_s_out, Ts, Tm, x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RecieverDownComer = FFD();
RecieverDownComer.L = 10;
HotBinBypass = FFD();
HotBinBypass.L = 10;
RecieverBypass = FFD();
RecieverBypass.L = 10;
HeaterDischarge = FFD();
HeaterDischarge.L = 10;
HotBinInlet = FFD();
HotBinInlet.L = 10;
HotBinDischarge = FFD();
HotBinDischarge.L = 10;
IntStorageDownComer = FFD();
IntStorageDownComer.L = 10;
HeatExchangerDownComer = FFD();
HeatExchangerDownComer.L = 10;
ColdBinBypass = FFD();
ColdBinBypass.L = 10;
ColdBinDischarge = FFD();
ColdBinDischarge.L = 10;
BucketElevatorDownComer = FFD();
BucketElevatorDownComer.L = 10;

%% load TMY3 data
weather = readtable('ABQ_Weather_Lookup.xlsx');

% Winter Week (Dec 21 - Dec 28)
% t_h_start = 8497;
% t_h_end = 8665;


% Spring Week (March 21 - March 28)
% t_h_start = 2641;
% t_h_end = 2809;


% Summer Week (June 21 - June 28)
% t_h_start = 4105;
% t_h_end = 4273;


% Fall Week (September 23 - September 30)
% t_h_start = 6361;
% t_h_end = 6529;



%% Set simulation parameters

% setup cycles (time units in seconds)
chargeDurration = 6*3600;
holdDurration = 10*3600;
dischargeDurration = 8*3600;
numCycles = 7;

dt = 600;
t = 0:dt:(chargeDurration + holdDurration + dischargeDurration)*numCycles;
nc = ceil(chargeDurration/dt);
nh = ceil(holdDurration/dt);
nd = ceil(dischargeDurration/dt);
nCycle = nc + nh + nd;

mdot = zeros(1, length(t));
for i = 1:numCycles
    for j = 1:nCycle
        if j <= nc
            mdot(j + nCycle*(i-1)) = 9;
        elseif j <= nh + nc && j > nc
            mdot(j + nCycle*(i-1)) = 0;                   
        else
            mdot(j + nCycle*(i-1)) = -5;
        end
    end    
end

Tin = 800; % + 20*sin(pi*t/3600);

y1 = zeros(size(t)); y1(1) = TES_.T0;
y2 = zeros(size(t)); y2(1) = TES_.T0;
y3 = zeros(size(t));
y4 = zeros(size(t)); y4(1) = 0.1;
profile on;
for i = 2:length(t)
    [y1(i), y2(i), y3(i), y4(i)] = step(TES_, Tin, mdot(i), t(i));
end
profsave;



