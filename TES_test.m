TES_ = TES();
TES_.T0 = 800;

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



