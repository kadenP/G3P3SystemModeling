FPR_ = FPR();
FPR_.Ts0 = 25;
t = 0:1:3600;
Ts_in = 400;
Tinf = 25;
mdot_s_in = 4;
Ts_out = zeros(size(t)); Ts_out(1) = FPR_.Ts0;
mdot_s_out = zeros(size(t)); mdot_s_out(1) = mdot_s_in;
Qsolar = 2.25e6;
for i = 2:length(t)
    [Ts_out(i), mdot_s_out(i)] = step(FPR_, Ts_in, Tinf, mdot_s_in, Qsolar, t(i));
end