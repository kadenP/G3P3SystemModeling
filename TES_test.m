TES_ = TES();
TES_.T0 = 800;
t = 0:1200:3600*6;
Tin = 700; % + 20*sin(pi*t/3600);
mdot = zeros(1, length(t));
mdot(1:end) = -5; %9 + 2*sin(pi*t/3600);
y1 = zeros(size(t));
y2 = zeros(size(t));
y3 = zeros(size(t));
y4 = zeros(size(t));
for i = 2:length(t)
    [y1(i), y2(i), y3(i), y4(i)] = step(TES_, Tin, mdot(i), t(i));
end