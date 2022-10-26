FPR_ = FPR();
FPR_.Ts0 = 600;
FPR_.tauLag = 600;
FPR_.tauLead = 1e-6;
FPR_.dtLin = 10;
t = 0:10:3600;
Tset_ = zeros(size(t));
% Tset_(1:60) = linspace(600, 775, 60);
Tset_(:) = 775;
Ts_in = 600;
Tinf = 25;
Qsolar = linspace(1e6, 2.25e6, 10);
Qsolar = Qsolar(10);
Ts_out = zeros(length(t), 1); Ts_out(1, :) = FPR_.Ts0;
Ts_out_lin = zeros(length(t), 1); Ts_out_lin(1, :) = FPR_.Ts0;
mdot_s_in = 9; 
mdot_s_out = zeros(size(t));
for i = 1:length(t)    
    [Ts_out(i), Ts_out_lin(i), mdot_s_out(i)] = ...
        step(FPR_, Ts_in, Tinf, Qsolar, t(i), Tset_(i));
end