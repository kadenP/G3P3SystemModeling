FPR_ = FPR();
FPR_.Ts0 = 25;
t = 0:1:3600;
Ts_in = 600;
Tset = 775;
Tinf = 25;
Ts_out = zeros(size(t)); Ts_out(1) = FPR_.Ts0;
mdot_s_in = zeros(size(t)); 
mdot_s_out = zeros(size(t));
Qsolar = 2.25e6;
for i = 2:length(t)
    [Ts_out(i), mdot_s_out(i)] = step(FPR_, Ts_in, Tinf, Tset, Qsolar, t(i));
end