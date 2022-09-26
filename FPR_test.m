FPR_ = FPR();
FPR_.Ts0 = 600;
t = 0:1:3600;
Ts_in = 600;
Tset = 775;
Tinf = 25;
Qsolar = linspace(1e6, 2.25e6, 10);
Qsolar = Qsolar(10)
Ts_out = zeros(length(t), 1); Ts_out(1, :) = FPR_.Ts0;
mdot_s_in = 9; 
mdot_s_out = zeros(size(t));
for i = 2:length(t)
    [Ts_out(i), mdot_s_out(i)] = step(FPR_, Ts_in, Tinf, mdot_s_in, Qsolar, t(i));
end
Ts_out10 = Ts_out;