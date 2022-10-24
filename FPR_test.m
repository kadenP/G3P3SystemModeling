FPR_ = FPR();
FPR_.Ts0 = 600;
FPR_.Tset = 775;
t = 0:1:3600;
Ts_in = 600;
Tinf = 25;
Qsolar = linspace(1e6, 2.25e6, 10);
Qsolar = Qsolar(10);
Ts_out = zeros(length(t), 1); Ts_out(1, :) = FPR_.Ts0;
Ts_out_lin = zeros(length(t), 1); Ts_out_lin(1, :) = FPR_.Ts0;
mdot_s_in = 9; 
mdot_s_out = zeros(size(t));
for i = 2:length(t)
    [Ts_out(i), Ts_out_lin(i), mdot_s_out(i)] = ...
        step(FPR_, Ts_in, Tinf, Qsolar, t(i));
end
Ts_out10 = Ts_out;