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
<<<<<<< HEAD
numCycles = 2;

dt = 60;
=======
numCycles = 7;

dt = 600;
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
t = 0:dt:(chargeDurration + holdDurration + dischargeDurration)*numCycles;

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
% Inputs: Ts_in, mdot_s_in, t
% Outputs: Ts_out, mdot_s_out, Ts, x_
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HX_MFH = MFH();

<<<<<<< HEAD
Ts_in_HX_MFH = HX_MFH.Ts0*ones(length(t), 1);
Ts_out_HX_MFH = HX_MFH.Ts0*ones(length(t), 1);
mdot_s_in_HX_MFH = zeros(length(t), 1);
mdot_s_out_HX_MFH = zeros(length(t), 1);
Ts_HX_MFH = cell(length(t), 1);
x_HX_MFH = cell(length(t), 1);

Reciever_MFH = MFH();

Ts_in_Reciever_MFH = Reciever_MFH.Ts0*ones(length(t), 1);
Ts_out_Reciever_MFH = Reciever_MFH.Ts0*ones(length(t), 1);
mdot_s_in_Reciever_MFH = zeros(length(t), 1);
mdot_s_out_Reciever_MFH = zeros(length(t), 1);
Ts_Reciever_MFH = cell(length(t), 1);
x_Reciever_MFH = cell(length(t), 1);
=======
Ts_in_HX_MFH = zeros(length(t), 1);
Ts_out_HX_MFH = zeros(length(t), 1);
mdot_s_in_HX_MFH = zeros(length(t), 1);
mdot_s_out_HX_MFH = zeros(length(t), 1);

Reciever_MFH = MFH();

Ts_in_Reciever_MFH = zeros(length(t), 1);
Ts_out_Reciever_MFH = zeros(length(t), 1);
mdot_s_in_Reciever_MFH = zeros(length(t), 1);
mdot_s_out_Reciever_MFH = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Falling Particle Reciever
% Inputs: Ts_in, Tinf, mdot_s_in, Qsolar, t
% Outputs: Ts_out, mdot_s_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FPR_ = FPR();

<<<<<<< HEAD
Ts_in_FPR = FPR_.Ts0*ones(length(t), 1);
mdot_s_in_FPR = zeros(length(t), 1);
Ts_out_FPR = FPR_.Ts0*ones(length(t), 1);
=======
Ts_in_FPR = zeros(length(t), 1);
mdot_s_in_FPR = zeros(length(t), 1);
Ts_out_FPR = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
mdot_s_out_FPR = zeros(length(t), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electric Heater
% Inputs: Ts_in, mdot_s_in, Qsolar, t
% Outputs: Ts_out, mdot_s_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Heater = EH();

Ts_in_Heater = zeros(length(t), 1);
mdot_in_Heater = zeros(length(t), 1);
Tset_Heater = zeros(length(t), 1);
Ts_out_Heater = zeros(length(t), 1);
mdot_out_Heater = zeros(length(t), 1);
Qin_Heater = zeros(length(t), 1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hot and Cold Storage Bins
% Inputs: Tin, mdot, t
% Outputs: Ts_out, T_bulk, Estored, ztop, ms, mdot_s_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hotTES = TES();
hotTES.T0 = 800;

Ts_in_hotTES = zeros(length(t), 1);
mdot_hotTES = zeros(length(t), 1);
Ts_out_hotTES = zeros(length(t), 1);
T_bulk_hotTES = zeros(length(t), 1);
Estored_hotTES = zeros(length(t), 1);
ztop_hotTES = zeros(length(t), 1);
ms_hotTES = zeros(length(t), 1);
mdot_s_out_hotTES = zeros(length(t), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intermediate Lumped Storage Bin
% Inputs: Tin, Tinf, mdot_s_in, mdot_s_out, t
% Outputs: Ts_out, ms, qloss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intTES = LSB();

<<<<<<< HEAD
Ts_in_intTES = intTES.Ts0*ones(length(t), 1);
mdot_s_in_intTES = zeros(length(t), 1);
mdot_s_out_intTES = zeros(length(t), 1);
Ts_out_intTES = intTES.Ts0*ones(length(t), 1);
=======
Ts_in_intTES = zeros(length(t), 1);
mdot_s_in_intTES = zeros(length(t), 1);
mdot_s_out_intTES = zeros(length(t), 1);
Ts_out_intTES = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
ms_intTES = zeros(length(t), 1);
qloss_intTES = zeros(length(t), 1);

coldTES = LSB();
coldTES.H = 7;
coldTES.D = 4.5;

<<<<<<< HEAD
Ts_in_coldTES = coldTES.Ts0*ones(length(t), 1);
mdot_s_in_coldTES = zeros(length(t), 1);
mdot_s_out_coldTES = zeros(length(t), 1);
Ts_out_coldTES = coldTES.Ts0*ones(length(t), 1);
=======
Ts_in_coldTES = zeros(length(t), 1);
mdot_s_in_coldTES = zeros(length(t), 1);
mdot_s_out_coldTES = zeros(length(t), 1);
Ts_out_coldTES = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
ms_coldTES = zeros(length(t), 1);
qloss_coldTES = zeros(length(t), 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particle-to-sCO2 Heat Exchanger
% Inputs: Ts_in, Tco2_in, mdot_s_in, mdot_CO2_in, t
% Outputs: Ts_out, Tco2_out, mdot_s_out, mdot_CO2_out, Ts, 
%          Tco2, Tm, Q_CO2, x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HX_ = HX();
HX_.n = 30;

<<<<<<< HEAD
Ts_in_HX = HX_.Ts0*ones(length(t), 1);
Tco2_in_HX = HX_.Tco20*ones(length(t), 1);
mdot_s_in_HX = zeros(length(t), 1);
mdot_CO2_in_HX = zeros(length(t), 1);
Ts_out_HX = HX_.Ts0*ones(length(t), 1);
Tco2_out_HX = HX_.Tco20*ones(length(t), 1);
=======
Ts_in_HX = zeros(length(t), 1);
Tco2_in_HX = zeros(length(t), 1);
mdot_s_in_HX = zeros(length(t), 1);
mdot_CO2_in_HX = zeros(length(t), 1);
Ts_out_HX = zeros(length(t), 1);
Tco2_out_HX = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
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

<<<<<<< HEAD
Ts_in_BE = BE_.Ts0*ones(length(t), 1);
mdot_s_in_BE = zeros(length(t), 1);
Ts_out_BE = BE_.Ts0*ones(length(t), 1);
=======
Ts_in_BE = zeros(length(t), 1);
mdot_s_in_BE = zeros(length(t), 1);
Ts_out_BE = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
mdot_s_out_BE = zeros(length(t), 1);
Ts_BE = cell(length(t), 1);
Tm_BE = cell(length(t), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solid Splitters
% Inputs: Ts_in, mdot_in, y1
% Outputs: Ts_out1, Ts_out2, mdot_out1, mdot_out2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HotBinDiverter = SS();

Ts_in_HotBinDiverter = zeros(length(t), 1);
mdot_in_HotBinDiverter = zeros(length(t), 1);
y1_HotBinDiverter = zeros(length(t), 1);
Ts_out1_HotBinDiverter = zeros(length(t), 1);
Ts_out2_HotBinDiverter = zeros(length(t), 1);
mdot_out1_HotBinDiverter = zeros(length(t), 1);
mdot_out2_HotBinDiverter = zeros(length(t), 1);

RecieverDiverter = SS();

Ts_in_RecieverDiverter = zeros(length(t), 1);
mdot_in_RecieverDiverter = zeros(length(t), 1);
y1_RecieverDiverter = zeros(length(t), 1);
Ts_out1_RecieverDiverter = zeros(length(t), 1);
Ts_out2_RecieverDiverter = zeros(length(t), 1);
mdot_out1_RecieverDiverter = zeros(length(t), 1);
mdot_out2_RecieverDiverter = zeros(length(t), 1);

IntermediateStorageDiverter = SS();

Ts_in_IntermediateStorageDiverter = zeros(length(t), 1);
mdot_in_IntermediateStorageDiverter = zeros(length(t), 1);
y1_IntermediateStorageDiverter = zeros(length(t), 1);
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

<<<<<<< HEAD
Ts_in_RecieverDownComer = RecieverDownComer.Ts0*ones(length(t), 1);
mdot_s_in_RecieverDownComer = zeros(length(t), 1);
Ts_out_RecieverDownComer = RecieverDownComer.Ts0*ones(length(t), 1);
=======
Ts_in_RecieverDownComer = zeros(length(t), 1);
mdot_s_in_RecieverDownComer = zeros(length(t), 1);
Ts_out_RecieverDownComer = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
mdot_s_out_RecieverDownComer = zeros(length(t), 1);
Ts_RecieverDownComer = cell(length(t), 1);
Tm_RecieverDownComer = cell(length(t), 1);
x_RecieverDownComer = cell(length(t), 1);

HotBinBypass = FFD();
HotBinBypass.L = 10;

<<<<<<< HEAD
Ts_in_HotBinBypass = HotBinBypass.Ts0*ones(length(t), 1);
mdot_s_in_HotBinBypass = zeros(length(t), 1);
Ts_out_HotBinBypass = HotBinBypass.Ts0*ones(length(t), 1);
=======
Ts_in_HotBinBypass = zeros(length(t), 1);
mdot_s_in_HotBinBypass = zeros(length(t), 1);
Ts_out_HotBinBypass = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
mdot_s_out_HotBinBypass = zeros(length(t), 1);
Ts_HotBinBypass = cell(length(t), 1);
Tm_HotBinBypass = cell(length(t), 1);
x_HotBinBypass = cell(length(t), 1);

RecieverBypass = FFD();
RecieverBypass.L = 10;

<<<<<<< HEAD
Ts_in_RecieverBypass = RecieverBypass.Ts0*ones(length(t), 1);
mdot_s_in_RecieverBypass = zeros(length(t), 1);
Ts_out_RecieverBypass = RecieverBypass.Ts0*ones(length(t), 1);
=======
Ts_in_RecieverBypass = zeros(length(t), 1);
mdot_s_in_RecieverBypass = zeros(length(t), 1);
Ts_out_RecieverBypass = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
mdot_s_out_RecieverBypass = zeros(length(t), 1);
Ts_RecieverBypass = cell(length(t), 1);
Tm_RecieverBypass = cell(length(t), 1);
x_RecieverBypass = cell(length(t), 1);

HeaterDischarge = FFD();
HeaterDischarge.L = 10;

<<<<<<< HEAD
Ts_in_HeaterDischarge = HeaterDischarge.Ts0*ones(length(t), 1);
mdot_s_in_HeaterDischarge = zeros(length(t), 1);
Ts_out_HeaterDischarge = HeaterDischarge.Ts0*ones(length(t), 1);
=======
Ts_in_HeaterDischarge = zeros(length(t), 1);
mdot_s_in_HeaterDischarge = zeros(length(t), 1);
Ts_out_HeaterDischarge = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
mdot_s_out_HeaterDischarge = zeros(length(t), 1);
Ts_HeaterDischarge = cell(length(t), 1);
Tm_HeaterDischarge = cell(length(t), 1);
x_HeaterDischarge = cell(length(t), 1);

HotBinInlet = FFD();
HotBinInlet.L = 10;

<<<<<<< HEAD
Ts_in_HotBinInlet = HotBinInlet.Ts0*ones(length(t), 1);
mdot_s_in_HotBinInlet = zeros(length(t), 1);
Ts_out_HotBinInlet = HotBinInlet.Ts0*ones(length(t), 1);
=======
Ts_in_HotBinInlet = zeros(length(t), 1);
mdot_s_in_HotBinInlet = zeros(length(t), 1);
Ts_out_HotBinInlet = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
mdot_s_out_HotBinInlet = zeros(length(t), 1);
Ts_HotBinInlet = cell(length(t), 1);
Tm_HotBinInlet = cell(length(t), 1);
x_HotBinInlet = cell(length(t), 1);

HotBinDischarge = FFD();
HotBinDischarge.L = 10;

<<<<<<< HEAD
Ts_in_HotBinDischarge = HotBinDischarge.Ts0*ones(length(t), 1);
mdot_s_in_HotBinDischarge = zeros(length(t), 1);
Ts_out_HotBinDischarge = HotBinDischarge.Ts0*ones(length(t), 1);
=======
Ts_in_HotBinDischarge = zeros(length(t), 1);
mdot_s_in_HotBinDischarge = zeros(length(t), 1);
Ts_out_HotBinDischarge = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
mdot_s_out_HotBinDischarge = zeros(length(t), 1);
Ts_HotBinDischarge = cell(length(t), 1);
Tm_HotBinDischarge = cell(length(t), 1);
x_HotBinDischarge = cell(length(t), 1);

IntStorageDownComer = FFD();
IntStorageDownComer.L = 10;

<<<<<<< HEAD
Ts_in_IntStorageDownComer = IntStorageDownComer.Ts0*ones(length(t), 1);
mdot_s_in_IntStorageDownComer = zeros(length(t), 1);
Ts_out_IntStorageDownComer = IntStorageDownComer.Ts0*ones(length(t), 1);
=======
Ts_in_IntStorageDownComer = zeros(length(t), 1);
mdot_s_in_IntStorageDownComer = zeros(length(t), 1);
Ts_out_IntStorageDownComer = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
mdot_s_out_IntStorageDownComer = zeros(length(t), 1);
Ts_IntStorageDownComer = cell(length(t), 1);
Tm_IntStorageDownComer = cell(length(t), 1);
x_IntStorageDownComer = cell(length(t), 1);

HeatExchangerDownComer = FFD();
HeatExchangerDownComer.L = 10;

<<<<<<< HEAD
Ts_in_HeatExchangerDownComer = HeatExchangerDownComer.Ts0*ones(length(t), 1);
mdot_s_in_HeatExchangerDownComer = zeros(length(t), 1);
Ts_out_HeatExchangerDownComer = HeatExchangerDownComer.Ts0*ones(length(t), 1);
=======
Ts_in_HeatExchangerDownComer = zeros(length(t), 1);
mdot_s_in_HeatExchangerDownComer = zeros(length(t), 1);
Ts_out_HeatExchangerDownComer = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
mdot_s_out_HeatExchangerDownComer = zeros(length(t), 1);
Ts_HeatExchangerDownComer = cell(length(t), 1);
Tm_HeatExchangerDownComer = cell(length(t), 1);
x_HeatExchangerDownComer = cell(length(t), 1);

ColdBinBypass = FFD();
ColdBinBypass.L = 10;

<<<<<<< HEAD
Ts_in_ColdBinBypass = ColdBinBypass.Ts0*ones(length(t), 1);
mdot_s_in_ColdBinBypass = zeros(length(t), 1);
Ts_out_ColdBinBypass = ColdBinBypass.Ts0*ones(length(t), 1);
=======
Ts_in_ColdBinBypass = zeros(length(t), 1);
mdot_s_in_ColdBinBypass = zeros(length(t), 1);
Ts_out_ColdBinBypass = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
mdot_s_out_ColdBinBypass = zeros(length(t), 1);
Ts_ColdBinBypass = cell(length(t), 1);
Tm_ColdBinBypass = cell(length(t), 1);
x_ColdBinBypass = cell(length(t), 1);

ColdBinDischarge = FFD();
ColdBinDischarge.L = 10;

<<<<<<< HEAD
Ts_in_ColdBinDischarge = ColdBinDischarge.Ts0*ones(length(t), 1);
mdot_s_in_ColdBinDischarge = zeros(length(t), 1);
Ts_out_ColdBinDischarge = ColdBinDischarge.Ts0*ones(length(t), 1);
=======
Ts_in_ColdBinDischarge = zeros(length(t), 1);
mdot_s_in_ColdBinDischarge = zeros(length(t), 1);
Ts_out_ColdBinDischarge = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
mdot_s_out_ColdBinDischarge = zeros(length(t), 1);
Ts_ColdBinDischarge = cell(length(t), 1);
Tm_ColdBinDischarge = cell(length(t), 1);
x_ColdBinDischarge = cell(length(t), 1);

BucketElevatorDownComer = FFD();
BucketElevatorDownComer.L = 10;

<<<<<<< HEAD
Ts_in_BucketElevatorDownComer = BucketElevatorDownComer.Ts0*ones(length(t), 1);
mdot_s_in_BucketElevatorDownComer = zeros(length(t), 1);
Ts_out_BucketElevatorDownComer = BucketElevatorDownComer.Ts0*ones(length(t), 1);
=======
Ts_in_BucketElevatorDownComer = zeros(length(t), 1);
mdot_s_in_BucketElevatorDownComer = zeros(length(t), 1);
Ts_out_BucketElevatorDownComer = zeros(length(t), 1);
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
mdot_s_out_BucketElevatorDownComer = zeros(length(t), 1);
Ts_BucketElevatorDownComer = cell(length(t), 1);
Tm_BucketElevatorDownComer = cell(length(t), 1);
x_BucketElevatorDownComer = cell(length(t), 1);

%% System data table

<<<<<<< HEAD
sysData = table(Qsolar, Ts_in_HX_MFH, Ts_out_HX_MFH, mdot_s_in_HX_MFH, ...
    mdot_s_out_HX_MFH, Ts_HX_MFH, x_HX_MFH, Ts_in_Reciever_MFH, Ts_out_Reciever_MFH, ...
    mdot_s_in_Reciever_MFH, mdot_s_out_Reciever_MFH, Ts_Reciever_MFH, x_Reciever_MFH, Ts_in_FPR, ...
=======
Ts_in_Heater = zeros(length(t), 1);
mdot_in_Heater = zeros(length(t), 1);
Tset_Heater = zeros(length(t), 1);
Ts_out_Heater = zeros(length(t), 1);
mdot_out_Heater = zeros(length(t), 1);
Qin_Heater = zeros(length(t), 1);

sysData = table(Qsolar, Ts_in_HX_MFH, Ts_out_HX_MFH, mdot_s_in_HX_MFH, ...
    mdot_s_out_HX_MFH, Ts_in_Reciever_MFH, Ts_out_Reciever_MFH, ...
    mdot_s_in_Reciever_MFH, mdot_s_out_Reciever_MFH, Ts_in_FPR, ...
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
    mdot_s_in_FPR, Ts_out_FPR, mdot_s_out_FPR, Ts_in_Heater, ...
    mdot_in_Heater, Tset_Heater, Ts_out_Heater, mdot_out_Heater, ...
    Qin_Heater, Ts_in_hotTES, mdot_hotTES, Ts_out_hotTES, ...
    T_bulk_hotTES, Estored_hotTES, ztop_hotTES, ms_hotTES,  ...
    mdot_s_out_hotTES, Ts_in_coldTES, mdot_s_in_coldTES, mdot_s_out_coldTES, ...
    Ts_out_coldTES, ms_coldTES, qloss_coldTES, Ts_in_intTES, ...
    mdot_s_in_intTES, mdot_s_out_intTES, Ts_out_intTES, ms_intTES, qloss_intTES, Ts_in_HX, ...
    Tco2_in_HX, mdot_s_in_HX, mdot_CO2_in_HX, Ts_out_HX, Tco2_out_HX, ...
    mdot_s_out_HX, mdot_CO2_out_HX, Ts_HX, Tco2_HX, Tm_HX, Q_CO2, ...
    x_HX, Ts_in_BE, mdot_s_in_BE, Ts_out_BE, mdot_s_out_BE, Ts_BE, ...
    Tm_BE, Ts_in_HotBinDiverter, mdot_in_HotBinDiverter, y1_HotBinDiverter, ...
    Ts_out1_HotBinDiverter, Ts_out2_HotBinDiverter, ...
    mdot_out1_HotBinDiverter, mdot_out2_HotBinDiverter, ...
    Ts_in_RecieverDiverter, mdot_in_RecieverDiverter, y1_RecieverDiverter, ...
    Ts_out1_RecieverDiverter, Ts_out2_RecieverDiverter, ...
    mdot_out1_RecieverDiverter, mdot_out2_RecieverDiverter, ...
    Ts_in_IntermediateStorageDiverter, mdot_in_IntermediateStorageDiverter, y1_IntermediateStorageDiverter, ...
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


%% load TMY3 data and set simulation dates
weather = readtable('ABQ_Weather_Lookup.xlsx');

% Winter Week (Dec 21 - Dec 28)
t_h_start = 8497; [~, ti_start] = min(abs(weather.time - t_h_start));
t_h_end = 8665; [~, ti_end] = min(abs(weather.time - t_h_end));


% Spring Week (March 21 - March 28)
% t_h_start = 2641;
% t_h_end = 2809;


% Summer Week (June 21 - June 28)
% t_h_start = 4105;
% t_h_end = 4273;


% Fall Week (September 23 - September 30)
% t_h_start = 6361;
% t_h_end = 6529;


%% Simulate over prescribed time frame

% additional simulation data

hour_day = zeros(size(t));
DNI = zeros(size(t));
Tinf = zeros(size(t));
HX_op = zeros(size(t));

% start in charge mode
[~, iStart] = min(abs(9*3600 - t));

% initialize system parameters


% profile on;
for i = iStart:length(t)
    
    % load current weather conditions
    hour_day(i) = mod(t(i)/3600, 24);
    DNI(i) = interp1(weather.time, weather.DNI, t_h_start + t(i)/3600, 'makima');
    Tinf(i) = interp1(weather.time, weather.Tinf, t_h_start + t(i)/3600, 'makima');
    
    % set particle-to-sCO2 inlet temperature and flow rate
    sysData.Tco2_in(i) = 565;
    sysData.mdot_CO2_in_HX(i) = 5;
    sysData.mdot_CO2_out_HX(i) = sysData.mdot_CO2_in_HX(i);
    
    % reciever collection limit
<<<<<<< HEAD
    sysData.Qsolar(i) = step(HF_, DNI(i));
%     if sysData.ms_coldTES(i-1) > 1000
%         sysData.Qsolar(i) = step(HF_, DNI(i));
%     else
%         sysData.Qsolar(i) = 0;
%     end
    
    % hot bin diverter valve control
=======
    if sysData.ms_coldTES(i-1) > 1000
        sysData.Qsolar(i) = step(HF_, DNI(i));
    else
        sysData.Qsolar(i) = 0;
    end
    
    % reciever diverter valve control
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
%     if sysData.Ts_in_HotBinDiverter(i-1) > 765 && ...
%         sysData.Ts_in_HotBinDiverter(i-1) < 800
%         y1_HotBinDiverter(i) = 1;        
%     else
%         y1_HotBinDiverter(i) = 0;        
%     end
    
    % intermediate storage diverter valve control
    sysData.y1_intermediateStorageDiverter(i) = 1; % no flow to cold bin
    
    % hot storage bin and heat exchanger flow rate
    if hour_day(i) >= 1 && hour_day(i) < 9
        sysData.mdot_hotTES(i) = -5;
        sysData.y1_HotBinDiverter(i) = 0;
<<<<<<< HEAD
        sysData.mdot_s_out_Heater(i) = 0;        
    elseif hour_day(i) >= 9 && hour_day(i) < 15
        sysData.mdot_hotTES(i) = 9;
        sysData.mdot_s_in_FPR(i) = 9;
        
=======
        sysData.mdot_s_out_Heater(i) = 0;
    elseif hour_day(i) >= 9 && hour_day(i) < 15
        sysData.mdot_hotTES(i) = 9;
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
        sysData.y1_HotBinDiverter(i) = 1;
        sysData.mdot_s_out_Heater(i) = 0.1;
    else
        sysData.mdot_hotTES(i) = 0;
        sysData.y1_HotBinDiverter(i) = 0;
        sysData.mdot_s_out_Heater(i) = 0.1;
    end
    
<<<<<<< HEAD
    % simulate system starting at the falling particle reciever
    [sysData.Ts_out_FPR(i), sysData.mdot_s_out_FPR(i)] = ...
        step(FPR_, sysData.Ts_in_FPR(i), Tinf(i), sysData.mdot_s_in_FPR(i), ...
         sysData.Qsolar(i), t(i));
     
    sysData.Ts_in_RecieverDownComer(i) = sysData.Ts_out_FPR(i);
    sysData.mdot_s_in_RecieverDownComer(i) = sysData.mdot_s_out_FPR(i);
    
    [sysData.Ts_out_RecieverDownComer(i), sysData.mdot_s_out_RecieverDownComer(i), ...
        sysData.Ts_RecieverDownComer{i}, sysData.Tm_RecieverDownComer{i}, ...
        sysData.x_RecieverDownComer{i}] = step(RecieverDownComer, ...
        sysData.Ts_in_RecieverDownComer(i), ...
        sysData.mdot_s_in_RecieverDownComer(i), Tinf(i), t(i));
        

    
    
    
    
=======
    % add connections to fully define the system
>>>>>>> 4363ca3f991cf6451b390d57e33299f67f9ff04f
    sysData.mdot_s_in_HeaterDischarge(i) = sysData.mdot_s_out_Heater(i);
    
    
        
    
end
% profsave;



