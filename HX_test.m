HX_ = HX();
HX_.n = 30;
t = 0:10:3600;
HX_.Ts0 = 25;
HX_.Tco20 = 25;
HX_.h_s = 400;                   % (W/m2K) solid-wall heat transfer coefficient 
HX_.h_CO2 = 3000;                % (W/m2K) CO2 heat transfer coefficient 
Ts_in = 775;
Tco2_in = 550;
mdot_s_in = 5;
mdot_CO2_in = 5;
Ts_out = zeros(size(t)); Ts_out(1) = HX_.Ts0;
Tco2_out = zeros(size(t)); Tco2_out(1) = HX_.Tco20;
mdot_s_out = zeros(size(t)); mdot_s_out(1) = mdot_s_in;
mdot_CO2_out = zeros(size(t)); mdot_CO2_out(1) = mdot_CO2_in;
Ts = cell(size(t)); Ts{1} = HX_.Ts0*ones(HX_.n, 1);
Tco2 = cell(size(t)); Tco2{1} = HX_.Tco20*ones(HX_.n, 1);
Tm = cell(size(t)); Tm{1} = HX_.Tm0*ones(HX_.n, 1);
Q_CO2 = zeros(size(t));
Q_s = zeros(size(t));

for i = 2:length(t)
    [Ts_out(i), Tco2_out(i), mdot_s_out(i), mdot_CO2_out(i), Ts{i}, Tco2{i}, ...
                Tm{i}, Q_CO2(i), Q_s(i), x] = step(HX_, Ts_in, Tco2_in, mdot_s_in, ...
                mdot_CO2_in, t(i));
end