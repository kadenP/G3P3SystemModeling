BE_ = BE();
BE_.n = 30;
t = 0:60:3600*8;
BE_.H = 12.192;
BE_.W = 1.1049; 
BE_.L = 0.4064; 
BE_.Wb = 1.1049/37.75;
BE_.Lb = 0.4064/37.75;
BE_.hinf = 10;
BE_.hsw = 1e6;

% insulation
BE_.wallInsulation{1, 1} = 'mineral wool';
% BE_.wallInsulation{1, 2} = [0, 0]; % m (0")
% BE_.wallInsulation{1, 2} = [0, 0.5*0.1016]; % m (2")
BE_.wallInsulation{1, 2} = [0, 0.1016]; % m (4")
% BE_.wallInsulation{1, 2} = [0, 2*0.1016]; % m (8")
BE_.wallInsulation{1, 3} = 0.1;         % W/mK
BE_.wallInsulation{1, 4} = 100;         % kg/m3
BE_.wallInsulation{1, 5} = 1000;         % J/kgK

% ramp inlet temperature
tauR = 30*60;
Tss = 275;
Tinf = 20;
Tin_ = @(t) (Tss - Tinf)*(1 - exp(-t/tauR)) + Tinf;
Ts_in = Tin_(t);

% Ts_in = 600;

vs = 0.5;
mdot_s_in = vs*BE_.Wb*BE_.Lb*BE_.rho_s*BE_.phi_s;
Ts_out = zeros(size(t)); Ts_out(1) = BE_.Ts0;
qLoss = zeros(size(t));
mdot_s_out = zeros(size(t)); mdot_s_out(1) = mdot_s_in;
Ts = cell(size(t)); Ts{1} = BE_.Ts0*ones(BE_.n, 1);
Tm = cell(size(t)); Tm{1} = BE_.Tm0*ones(BE_.n, 1);

for i = 2:length(t)
    [Ts_out(i), mdot_s_out(i), Ts{i}, Tm{i}, x, qLoss(i)] = step(BE_, Ts_in(i), mdot_s_in, ...
                Tinf, t(i));
end


% residence time computation
Vdot = mdot_s_in/(BE_.rho_s*BE_.phi_s);
vs = Vdot/(BE_.Wb*BE_.Lb);
tres = BE_.H/vs;



