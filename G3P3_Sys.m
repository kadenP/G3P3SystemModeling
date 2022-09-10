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

%% Define all system components

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heliostat Field
% Inputs: DNI
% Outputs: Qsolar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HF_ = HF();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Falling Particle Reciever
% Inputs: Ts_in, Tinf, mdot_s_in, Qsolar, t
% Outputs: Ts_out, mdot_s_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FPR_ = FPR();
FPR_.Ts0 = 25;

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
SS1 = SS();
SS2 = SS();
SS3 = SS();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solid Junctions
% Inputs: Ts_in1, Ts_in2, mdot_in1, mdot_in2
% Outputs: Ts_out, mdot_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SJ1 = SJ();
SJ2 = SJ();
SJ3 = SJ();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Free Fall Ducts
% Inputs: Ts_in, mdot_s_in, Tinf, t
% Outputs: Ts_out, mdot_s_out, Ts, Tm, x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FFD1 = FFD();
FFD2 = FFD();
FFD3 = FFD();
FFD4 = FFD();
FFD5 = FFD();
FFD6 = FFD();
FFD7 = FFD();
FFD8 = FFD();
FFD9 = FFD();
FFD10 = FFD();
FFD11 = FFD();

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



