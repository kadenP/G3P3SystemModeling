LSB_ = LSB();
LSB_.ms0 = 30;
t = 0:10:3600;
Ts_in = 800;
Tinf = 20;
mdot_s_in = 9;
mdot_s_out = 8;
Ts_out = zeros(size(t)); Ts_out(1) = LSB_.Ts0;
qloss = zeros(size(t)); 

for i = 2:length(t)
    [Ts_out(i), qloss(i)] = step(LSB_, Ts_in, Tinf, mdot_s_in, mdot_s_out, t(i));
end