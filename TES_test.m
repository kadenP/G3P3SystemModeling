TES_ = TES();
TES_.T0 = 800;
t = 0:600:3600*6;
Tin = 800; % + 20*sin(pi*t/3600);
mdot = zeros(1, length(t));
mdot(1:end) = 9; %9 + 2*sin(pi*t/3600);
y1 = zeros(size(t)); y1(1) = TES_.T0;
y2 = zeros(size(t)); y2(1) = TES_.T0;
y3 = zeros(size(t));
y4 = zeros(size(t));
for i = 2:length(t)
    [y1(i), y2(i), y3(i), y4(i)] = step(TES_, Tin, mdot(i), t(i));
end