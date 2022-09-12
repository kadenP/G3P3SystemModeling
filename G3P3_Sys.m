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

dt = 600;
t = 0:dt:(chargeDurration + holdDurration + dischargeDurration)*numCycles;
nc = ceil(chargeDurration/dt);
nh = ceil(holdDurration/dt);
nd = ceil(dischargeDurration/dt);
nCycle = nc + nh + nd;

%% Define all system components (note that default model parameters are set to match the Modelica G3P3 model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heliostat Field
% Inputs: DNI
% Outputs: Qsolar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HF_ = HF();

Qsolar = zeros(length(t), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass Flow Hopper
% Inputs: Ts_in, Tinf, mdot_s_in, Qsolar, t
% Outputs: Ts_out, mdot_s_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HX_MFH = MFH();

Ts_in_HX_MFH = zeros(length(t), 1);
Ts_out_HX_MFH = zeros(length(t), 1);
mdot_s_in_HX_MFH = zeros(length(t), 1);
mdot_s_out_HX_MFH = zeros(length(t), 1);

Reciever_MFH = MFH();

Ts_in_Reciever_MFH = zeros(length(t), 1);
Ts_out_Reciever_MFH = zeros(length(t), 1);
mdot_s_in_Reciever_MFH = zeros(length(t), 1);
mdot_s_out_Reciever_MFH = zeros(length(t), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Falling Particle Reciever
% Inputs: Ts_in, Tinf, mdot_s_in, Qsolar, t
% Outputs: Ts_out, mdot_s_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FPR_ = FPR();

Ts_in_FPR = zeros(length(t), 1);
mdot_s_in_FPR = zeros(length(t), 1);
Ts_out_FPR = zeros(length(t), 1);
mdot_s_out_FPR = zeros(length(t), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hot and Cold Storage Bins
% Inputs: Tin, mdot, t
% Outputs: Ts_out, T_bulk, Estored, ztop, mdot_s_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hotTES = TES();
hotTES.T0 = 800;

Ts_in_hotTES = zeros(length(t), 1);
mdot_hotTES = zeros(length(t), 1);
Ts_out_hotTES = zeros(length(t), 1);
T_bulk_hotTES = zeros(length(t), 1);
Estored_hotTES = zeros(length(t), 1);
ztop_hotTES = zeros(length(t), 1);

coldTES = TES();
coldTES.T0 = 800;

Ts_in_coldTES = zeros(length(t), 1);
mdot_coldTES = zeros(length(t), 1);
Ts_out_coldTES = zeros(length(t), 1);
T_bulk_coldTES = zeros(length(t), 1);
Estored_coldTES = zeros(length(t), 1);
ztop_coldTES = zeros(length(t), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intermediate Lumped Storage Bin
% Inputs: Tin, Tinf, mdot_s_in, mdot_s_out, t
% Outputs: Ts_out, qloss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intTES = LSB();

Ts_in_intTES = zeros(length(t), 1);
mdot_s_in = zeros(length(t), 1);
mdot_s_out = zeros(length(t), 1);
Ts_out_intTES = zeros(length(t), 1);
qloss_intTES = zeros(length(t), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particle-to-sCO2 Heat Exchanger
% Inputs: Ts_in, Tco2_in, mdot_s_in, mdot_CO2_in, t
% Outputs: Ts_out, Tco2_out, mdot_s_out, mdot_CO2_out, Ts, 
%          Tco2, Tm, Q_CO2, x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HX_ = HX();
HX_.n = 30;

Ts_in_HX = zeros(length(t), 1);
Tco2_in_HX = zeros(length(t), 1);
mdot_s_in_HX = zeros(length(t), 1);
mdot_CO2_HX = zeros(length(t), 1);
Ts_out_HX = zeros(length(t), 1);
Tco2_out_HX = zeros(length(t), 1);
mdot_s_out_HX = zeros(length(t), 1);
mdot_CO2_out_HX = zeros(length(t), 1);
Ts_HX = cell(length(t), 1);
Tco2_HX = cell(length(t), 1);
Tm_HX = cell(length(t), 1);
Q_CO2 = zeros(length(t), 1);
x_HX = zeros(length(t), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bucket Elevator
% Inputs: Ts_in, mdot_s_in, Tinf, t
% Outputs: Ts_out, mdot_s_out, Ts, Tm, x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BE_ = BE();
BE_.n = 30;

Ts_in_BE = zeros(length(t), 1);
mdot_s_in_BE = zeros(length(t), 1);
Ts_out_BE = zeros(length(t), 1);
mdot_s_out_BE = zeros(length(t), 1);
Ts_BE = cell(length(t), 1);
Tm_BE = cell(length(t), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solid Splitters
% Inputs: Ts_in, mdot_in
% Outputs: Ts_out1, Ts_out2, mdot_out1, mdot_out2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HotBinDiverter = SS();

Ts_in_HotBinDiverter = zeros(length(t), 1);
mdot_in_HotBinDiverter = zeros(length(t), 1);
Ts_out1_HotBinDiverter = zeros(length(t), 1);
Ts_out2_HotBinDiverter = zeros(length(t), 1);
mdot_out1_HotBinDiverter = zeros(length(t), 1);
mdot_out2_HotBinDiverter = zeros(length(t), 1);

RecieverDiverter = SS();

Ts_in_RecieverDiverter = zeros(length(t), 1);
mdot_in_RecieverDiverter = zeros(length(t), 1);
Ts_out1_RecieverDiverter = zeros(length(t), 1);
Ts_out2_RecieverDiverter = zeros(length(t), 1);
mdot_out1_RecieverDiverter = zeros(length(t), 1);
mdot_out2_RecieverDiverter = zeros(length(t), 1);

IntermediateStorageDiverter = SS();

Ts_in_IntermediateStorageDiverter = zeros(length(t), 1);
mdot_in_IntermediateStorageDiverter = zeros(length(t), 1);
Ts_out1_IntermediateStorageDiverter = zeros(length(t), 1);
Ts_out2_IntermediateStorageDiverter = zeros(length(t), 1);
mdot_out1_IntermediateStorageDiverter = zeros(length(t), 1);
mdot_out2_IntermediateStorageDiverter = zeros(length(t), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solid Junctions
% Inputs: Ts_in1, Ts_in2, mdot_in1, mdot_in2
% Outputs: Ts_out, mdot_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HeatExchangerJunction = SJ();

Ts_in1_HeatExchangerJunction = zeros(length(t), 1);
Ts_in2_HeatExchangerJunction = zeros(length(t), 1);
mdot_in1_HeatExchangerJunction = zeros(length(t), 1);
mdot_in2_HeatExchangerJunction = zeros(length(t), 1);
Ts_out_HeatExchangerJunction = zeros(length(t), 1);
mdot_out_HeatExchangerJunction = zeros(length(t), 1);

ColdBinJunction = SJ();

Ts_in1_ColdBinJunction = zeros(length(t), 1);
Ts_in2_ColdBinJunction = zeros(length(t), 1);
mdot_in1_ColdBinJunction = zeros(length(t), 1);
mdot_in2_ColdBinJunction = zeros(length(t), 1);
Ts_out_ColdBinJunction = zeros(length(t), 1);
mdot_out_ColdBinJunction = zeros(length(t), 1);

BucketElevatorJunction = SJ();

Ts_in1_BucketElevatorJunction = zeros(length(t), 1);
Ts_in2_BucketElevatorJunction = zeros(length(t), 1);
mdot_in1_BucketElevatorJunction = zeros(length(t), 1);
mdot_in2_BucketElevatorJunction = zeros(length(t), 1);
Ts_out_BucketElevatorJunction = zeros(length(t), 1);
mdot_out_BucketElevatorJunction = zeros(length(t), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Free Fall Ducts
% Inputs: Ts_in, mdot_s_in, Tinf, t
% Outputs: Ts_out, mdot_s_out, Ts, Tm, x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RecieverDownComer = FFD();
RecieverDownComer.L = 10;

Ts_in_RecieverDownComer = zeros(length(t), 1);
mdot_s_in_RecieverDownComer = zeros(length(t), 1);
Ts_out_RecieverDownComer = zeros(length(t), 1);
mdot_s_out_RecieverDownComer = zeros(length(t), 1);
Ts_RecieverDownComer = cell(length(t), 1);
Tm_RecieverDownComer = cell(length(t), 1);
x_RecieverDownComer = cell(length(t), 1);

HotBinBypass = FFD();
HotBinBypass.L = 10;

Ts_in_HotBinBypass = zeros(length(t), 1);
mdot_s_in_HotBinBypass = zeros(length(t), 1);
Ts_out_HotBinBypass = zeros(length(t), 1);
mdot_s_out_HotBinBypass = zeros(length(t), 1);
Ts_HotBinBypass = cell(length(t), 1);
Tm_HotBinBypass = cell(length(t), 1);
x_HotBinBypass = cell(length(t), 1);

RecieverBypass = FFD();
RecieverBypass.L = 10;

Ts_in_RecieverBypass = zeros(length(t), 1);
mdot_s_in_RecieverBypass = zeros(length(t), 1);
Ts_out_RecieverBypass = zeros(length(t), 1);
mdot_s_out_RecieverBypass = zeros(length(t), 1);
Ts_RecieverBypass = cell(length(t), 1);
Tm_RecieverBypass = cell(length(t), 1);
x_RecieverBypass = cell(length(t), 1);

HeaterDischarge = FFD();
HeaterDischarge.L = 10;

Ts_in_HeaterDischarge = zeros(length(t), 1);
mdot_s_in_HeaterDischarge = zeros(length(t), 1);
Ts_out_HeaterDischarge = zeros(length(t), 1);
mdot_s_out_HeaterDischarge = zeros(length(t), 1);
Ts_HeaterDischarge = cell(length(t), 1);
Tm_HeaterDischarge = cell(length(t), 1);
x_HeaterDischarge = cell(length(t), 1);

HotBinInlet = FFD();
HotBinInlet.L = 10;

Ts_in_HotBinInlet = zeros(length(t), 1);
mdot_s_in_HotBinInlet = zeros(length(t), 1);
Ts_out_HotBinInlet = zeros(length(t), 1);
mdot_s_out_HotBinInlet = zeros(length(t), 1);
Ts_HotBinInlet = cell(length(t), 1);
Tm_HotBinInlet = cell(length(t), 1);
x_HotBinInlet = cell(length(t), 1);

HotBinDischarge = FFD();
HotBinDischarge.L = 10;

Ts_in_HotBinDischarge = zeros(length(t), 1);
mdot_s_in_HotBinDischarge = zeros(length(t), 1);
Ts_out_HotBinDischarge = zeros(length(t), 1);
mdot_s_out_HotBinDischarge = zeros(length(t), 1);
Ts_HotBinDischarge = cell(length(t), 1);
Tm_HotBinDischarge = cell(length(t), 1);
x_HotBinDischarge = cell(length(t), 1);

IntStorageDownComer = FFD();
IntStorageDownComer.L = 10;

Ts_in_IntStorageDownComer = zeros(length(t), 1);
mdot_s_in_IntStorageDownComer = zeros(length(t), 1);
Ts_out_IntStorageDownComer = zeros(length(t), 1);
mdot_s_out_IntStorageDownComer = zeros(length(t), 1);
Ts_IntStorageDownComer = cell(length(t), 1);
Tm_IntStorageDownComer = cell(length(t), 1);
x_IntStorageDownComer = cell(length(t), 1);

HeatExchangerDownComer = FFD();
HeatExchangerDownComer.L = 10;

Ts_in_HeatExchangerDownComer = zeros(length(t), 1);
mdot_s_in_HeatExchangerDownComer = zeros(length(t), 1);
Ts_out_HeatExchangerDownComer = zeros(length(t), 1);
mdot_s_out_HeatExchangerDownComer = zeros(length(t), 1);
Ts_HeatExchangerDownComer = cell(length(t), 1);
Tm_HeatExchangerDownComer = cell(length(t), 1);
x_HeatExchangerDownComer = cell(length(t), 1);

ColdBinBypass = FFD();
ColdBinBypass.L = 10;

Ts_in_ColdBinBypass = zeros(length(t), 1);
mdot_s_in_ColdBinBypass = zeros(length(t), 1);
Ts_out_ColdBinBypass = zeros(length(t), 1);
mdot_s_out_ColdBinBypass = zeros(length(t), 1);
Ts_ColdBinBypass = cell(length(t), 1);
Tm_ColdBinBypass = cell(length(t), 1);
x_ColdBinBypass = cell(length(t), 1);

ColdBinDischarge = FFD();
ColdBinDischarge.L = 10;

Ts_in_ColdBinDischarge = zeros(length(t), 1);
mdot_s_in_ColdBinDischarge = zeros(length(t), 1);
Ts_out_ColdBinDischarge = zeros(length(t), 1);
mdot_s_out_ColdBinDischarge = zeros(length(t), 1);
Ts_ColdBinDischarge = cell(length(t), 1);
Tm_ColdBinDischarge = cell(length(t), 1);
x_ColdBinDischarge = cell(length(t), 1);

BucketElevatorDownComer = FFD();
BucketElevatorDownComer.L = 10;

Ts_in_BucketElevatorDownComer = zeros(length(t), 1);
mdot_s_in_BucketElevatorDownComer = zeros(length(t), 1);
Ts_out_BucketElevatorDownComer = zeros(length(t), 1);
mdot_s_out_BucketElevatorDownComer = zeros(length(t), 1);
Ts_BucketElevatorDownComer = cell(length(t), 1);
Tm_BucketElevatorDownComer = cell(length(t), 1);
x_BucketElevatorDownComer = cell(length(t), 1);

%% System data table

sysData = table(Qsolar, Ts_in_HX_MFH, Ts_out_HX_MFH, mdot_s_in_HX_MFH, ...
    mdot_s_out_HX_MFH, Ts_in_Reciever_MFH, Ts_out_Reciever_MFH, ...
    mdot_s_in_Reciever_MFH, mdot_s_out_Reciever_MFH, Ts_in_FPR, ...
    mdot_s_in_FPR, Ts_out_FPR, mdot_s_out_FPR, Ts_in_hotTES, ...
    mdot_hotTES, Ts_out_hotTES, T_bulk_hotTES, Estored_hotTES, ...
    ztop_hotTES, Ts_in_coldTES, mdot_coldTES, Ts_out_coldTES, ...
    T_bulk_coldTES, Estored_coldTES, ztop_coldTES, Ts_in_intTES, ...
    mdot_s_in, mdot_s_out, Ts_out_intTES, qloss_intTES, Ts_in_HX, ...
    Tco2_in_HX, mdot_s_in_HX, mdot_CO2_HX, Ts_out_HX, Tco2_out_HX, ...
    mdot_s_out_HX, mdot_CO2_out_HX, Ts_HX, Tco2_HX, Tm_HX, Q_CO2, ...
    x_HX, Ts_in_BE, mdot_s_in_BE, Ts_out_BE, mdot_s_out_BE, Ts_BE, ...
    Tm_BE, Ts_in_HotBinDiverter, mdot_in_HotBinDiverter, ...
    Ts_out1_HotBinDiverter, Ts_out2_HotBinDiverter, ...
    mdot_out1_HotBinDiverter, mdot_out2_HotBinDiverter, ...
    Ts_in_RecieverDiverter, mdot_in_RecieverDiverter, ...
    Ts_out1_RecieverDiverter, Ts_out2_RecieverDiverter, ...
    mdot_out1_RecieverDiverter, mdot_out2_RecieverDiverter, ...
    Ts_in_IntermediateStorageDiverter, mdot_in_IntermediateStorageDiverter, ...
    Ts_out1_IntermediateStorageDiverter, Ts_out2_IntermediateStorageDiverter, ...
    mdot_out1_IntermediateStorageDiverter, mdot_out2_IntermediateStorageDiverter, ...
    Ts_in1_HeatExchangerJunction, Ts_in2_HeatExchangerJunction, ...
    mdot_in1_HeatExchangerJunction, mdot_in2_HeatExchangerJunction, ...
    Ts_out_HeatExchangerJunction, mdot_out_HeatExchangerJunction, ...
    Ts_in1_ColdBinJunction, Ts_in2_ColdBinJunction, ...
    mdot_in1_ColdBinJunction, mdot_in2_ColdBinJunction, ...
    Ts_out_ColdBinJunction, mdot_out_ColdBinJunction, ...
    Ts_in1_BucketElevatorJunction, Ts_in2_BucketElevatorJunction, ...
    mdot_in1_BucketElevatorJunction, mdot_in2_BucketElevatorJunction, ...
    Ts_out_BucketElevatorJunction, mdot_out_BucketElevatorJunction, ...
    Ts_in_RecieverDownComer, mdot_s_in_RecieverDownComer, ...
    Ts_out_RecieverDownComer, mdot_s_out_RecieverDownComer, ...
    Ts_RecieverDownComer, Tm_RecieverDownComer, x_RecieverDownComer, ...
    Ts_in_HotBinBypass, mdot_s_in_HotBinBypass, Ts_out_HotBinBypass, ...
    mdot_s_out_HotBinBypass, Ts_HotBinBypass, Tm_HotBinBypass, ...
    x_HotBinBypass, Ts_in_RecieverBypass, mdot_s_in_RecieverBypass, ...
    Ts_out_RecieverBypass, mdot_s_out_RecieverBypass, ...
    Ts_RecieverBypass, Tm_RecieverBypass, x_RecieverBypass, ...
    Ts_in_HeaterDischarge, mdot_s_in_HeaterDischarge, ...
    Ts_out_HeaterDischarge, mdot_s_out_HeaterDischarge, ...
    Ts_HeaterDischarge, Tm_HeaterDischarge, x_HeaterDischarge, ...
    Ts_in_HotBinInlet, mdot_s_in_HotBinInlet, Ts_out_HotBinInlet, ...
    mdot_s_out_HotBinInlet, Ts_HotBinInlet, Tm_HotBinInlet, ...
    x_HotBinInlet, Ts_in_HotBinDischarge, mdot_s_in_HotBinDischarge, ...
    Ts_out_HotBinDischarge, mdot_s_out_HotBinDischarge, ...
    Ts_HotBinDischarge, Tm_HotBinDischarge, x_HotBinDischarge, ...
    Ts_in_IntStorageDownComer, mdot_s_in_IntStorageDownComer, ...
    Ts_out_IntStorageDownComer, mdot_s_out_IntStorageDownComer, ...
    Ts_IntStorageDownComer, Tm_IntStorageDownComer, x_IntStorageDownComer, ...
    Ts_in_HeatExchangerDownComer, mdot_s_in_HeatExchangerDownComer, ...
    Ts_out_HeatExchangerDownComer, mdot_s_out_HeatExchangerDownComer, ...
    Ts_HeatExchangerDownComer, Tm_HeatExchangerDownComer, ...
    x_HeatExchangerDownComer, Ts_in_ColdBinBypass, mdot_s_in_ColdBinBypass, ...
    Ts_out_ColdBinBypass, mdot_s_out_ColdBinBypass, Ts_ColdBinBypass, ...
    Tm_ColdBinBypass, x_ColdBinBypass, Ts_in_ColdBinDischarge, ...
    mdot_s_in_ColdBinDischarge, Ts_out_ColdBinDischarge, ...
    mdot_s_out_ColdBinDischarge, Ts_ColdBinDischarge, ...
    Tm_ColdBinDischarge, x_ColdBinDischarge, Ts_in_BucketElevatorDownComer, ...
    mdot_s_in_BucketElevatorDownComer, Ts_out_BucketElevatorDownComer, ...
    mdot_s_out_BucketElevatorDownComer, Ts_BucketElevatorDownComer, ...
    Ts_BucketElevatorDownComer, Tm_BucketElevatorDownComer, ...
    x_BucketElevatorDownComer);


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



