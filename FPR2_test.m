FPR_ = FPR2();
FPR_.Ts0 = 600;
FPR_.n = 250;
t = 0:10:3600;
Ts_in = 600;
Tset = 775;
Tinf = 25;
Ts_out = zeros(size(t)); Ts_out(1) = FPR_.Ts0;
Ts = cell(size(t)); Ts{1} = FPR_.Ts0*ones(FPR_.n, 1);
mdot_s_in = 9*ones(size(t)); 
mdot_s_out = zeros(size(t));
Qsolar = 2.25e6;
for i = 2:length(t)
    [Ts_out(i), mdot_s_out(i), Ts{i}, x_] = step(FPR_, Ts_in, Tinf, mdot_s_in(i), Qsolar, t(i));
end