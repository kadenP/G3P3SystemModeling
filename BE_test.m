BE_ = BE();
BE_.n = 30;
t = 0:10:3600;
Ts_in = 775;
Tinf = 20;
mdot_s_in = 9;
Ts_out = zeros(size(t)); Ts_out(1) = BE_.Ts0;
mdot_s_out = zeros(size(t)); mdot_s_out(1) = mdot_s_in;
Ts = cell(size(t)); Ts{1} = BE_.Ts0*ones(BE_.n, 1);
Tm = cell(size(t)); Tm{1} = BE_.Tm0*ones(BE_.n, 1);

for i = 2:length(t)
    [Ts_out(i), mdot_s_out(i), Ts{i}, Tm{i}, x] = step(BE_, Ts_in, mdot_s_in, ...
                Tinf, t(i));
end