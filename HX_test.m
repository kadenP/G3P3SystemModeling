

chargeDurration = 6*3600;
holdDurration = 10*3600;
dischargeDurration = 8*3600;
numCycles = 7;

dt = 10;
t = 0:dt:(chargeDurration + holdDurration + dischargeDurration)*numCycles;
t = t';
nc = ceil(chargeDurration/dt);
nh = ceil(holdDurration/dt);
nd = ceil(dischargeDurration/dt);
nCycle = nc + nh + nd;

HX_ = HX();
HX_.n = 30;
Tinf = 20;
HX_.Ts0 = 25;
HX_.Tco20 = 25;
HX_.h_s = 400;                   % (W/m2K) solid-wall heat transfer coefficient 
HX_.h_CO2 = 3000;                % (W/m2K) CO2 heat transfer coefficient 
mdot_s_in_HX = zeros(size(t)); mdot_s_in_HX(:) = 5;
mdot_CO2_in_HX = zeros(size(t)); mdot_CO2_in_HX(:) = 5;
Ts_in_HX = zeros(size(t)); Ts_in_HX(1) = HX_.Ts0;
Tco2_in_HX = zeros(size(t)); Tco2_in_HX(:) = 550;
Ts_out_HX = zeros(size(t)); Ts_out(1) = HX_.Ts0;
Tco2_out_HX = zeros(size(t)); Tco2_out(1) = HX_.Tco20;
mdot_s_out_HX = mdot_s_in_HX;
mdot_CO2_out_HX = mdot_CO2_in_HX;
Ts_HX = cell(size(t)); Ts{1} = HX_.Ts0*ones(HX_.n, 1);
Tco2_HX = cell(size(t)); Tco2{1} = HX_.Tco20*ones(HX_.n, 1);
Tm_HX = cell(size(t)); Tm{1} = HX_.Tm0*ones(HX_.n, 1);
Q_CO2 = zeros(size(t));
Q_s = zeros(size(t));
x_HX = cell(length(t), 1);

HotBinDischarge = FFD();
HotBinDischarge.L = 10;

Ts_in_HotBinDischarge = HotBinDischarge.Ts0*ones(length(t), 1);
mdot_s_in_HotBinDischarge = zeros(length(t), 1);
Ts_out_HotBinDischarge = HotBinDischarge.Ts0*ones(length(t), 1);
mdot_s_out_HotBinDischarge = zeros(length(t), 1);
Ts_HotBinDischarge = cell(length(t), 1);
Tm_HotBinDischarge = cell(length(t), 1);
x_HotBinDischarge = cell(length(t), 1);

load('sysData2_Summer.mat', 'TES_Summer');

sysData3 = table(Ts_in_HX, Tco2_in_HX, mdot_s_in_HX, mdot_CO2_in_HX, Ts_out_HX, Tco2_out_HX, ...
    mdot_s_out_HX, mdot_CO2_out_HX, Ts_HX, Tco2_HX, Tm_HX, Q_CO2, ...
    Q_s, x_HX, Ts_in_HotBinDischarge, mdot_s_in_HotBinDischarge, ...
    Ts_out_HotBinDischarge, mdot_s_out_HotBinDischarge, ...
    Ts_HotBinDischarge, Tm_HotBinDischarge, x_HotBinDischarge);

for i = 1:numCycles
    for j = 1:nCycle
        if j <= nc
            sysData3.mdot_s_in_HotBinDischarge(j + nCycle*(i-1)) = 5;
            sysData3.Ts_in_HotBinDischarge(j + nCycle*(i-1)) = 600;
        elseif j <= nh + nc && j > nc
            sysData3.mdot_s_in_HotBinDischarge(j + nCycle*(i-1)) = 5;
            sysData3.Ts_in_HotBinDischarge(j + nCycle*(i-1)) = 600;
        else
            sysData3.mdot_s_in_HotBinDischarge(j + nCycle*(i-1)) = 5;
            sysData3.Ts_in_HotBinDischarge(j + nCycle*(i-1)) = ...
                interp1(TES_Summer.t, TES_Summer.Ts_out, t(j + nCycle*(i-1)), 'makima');
        end
    end    
end

for i = 2:length(t)
    
    [sysData3.Ts_out_HotBinDischarge(i), sysData3.mdot_s_out_HotBinDischarge(i), ...
        sysData3.Ts_HotBinDischarge{i}, sysData3.Tm_HotBinDischarge{i}, ...
        sysData3.x_HotBinDischarge{i}] = step(HotBinDischarge, ...
        sysData3.Ts_in_HotBinDischarge(i), ...
        sysData3.mdot_s_in_HotBinDischarge(i), Tinf, t(i));
    
    sysData3.Ts_in_HX(i) = sysData3.Ts_out_HotBinDischarge(i);
    
    [sysData3.Ts_out_HX(i), sysData3.Tco2_out_HX(i), ...
        sysData3.mdot_s_out_HX(i), sysData3.mdot_CO2_out_HX(i), ...
        sysData3.Ts_HX{i}, sysData3.Tco2_HX{i}, sysData3.Tm_HX{i}, ...
        sysData3.Q_CO2(i), sysData3.Q_s(i), sysData3.x_HX{i}] ...
        = step(HX_, sysData3.Ts_in_HX(i), sysData3.Tco2_in_HX(i), sysData3.mdot_s_in_HX(i), ...
                sysData3.mdot_CO2_in_HX(i), t(i));
    fprintf('Tco2_out_HX = %1.2f �C\n', sysData3.Tco2_out_HX(i));
end