TES_ = TES();
TES_.T0 = 775;

% setup cycles (time units in seconds)
chargeDurration = 6*3600;
holdDurration = 10*3600;
dischargeDurration = 8*3600;
numCycles = 1;

dt = 600;
t = 0:dt:(chargeDurration + holdDurration + dischargeDurration)*numCycles;
nc = ceil(chargeDurration/dt);
nh = ceil(holdDurration/dt);
nd = ceil(dischargeDurration/dt);
nCycle = nc + nh + nd;

mdot = zeros(length(t), 1);
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

Tin = 775; % + 20*sin(pi*t/3600);
Tinf = 0;
Ts_out = zeros(size(t)); Ts_out(1) = TES_.T0;
Ts_bulk = zeros(size(t)); Ts_bulk(1) = TES_.T0;
Estored = zeros(size(t));
ztop_ = zeros(size(t)); ztop_(1) = 0.1;
ms = zeros(size(t));
mdot_s_out = zeros(size(t));
% profile on;
for i = 1:length(t)
    [Ts_out(i), Ts_bulk(i), Estored(i), ztop_(i), ms(i), mdot_s_out(i)] = ...
        step(TES_, Tin, Tinf, mdot(i), t(i));
end
% profsave;



