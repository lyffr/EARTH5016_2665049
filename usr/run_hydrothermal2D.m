% run_hydrothermal2D.m
% ***** PRACTICAL 30 RUN SCRIPT: KD0 sensitivity sweep ********************

clear; close all;

% --- Domain / grid ---
W  = 1e3;      % width  [m]
D  = 1e3;      % depth  [m]
Nz = 100;      % grid in z
h  = D/Nz;     % spacing (dx=dz=h)

% enforce periodic closure in x
Nx = round(W/h);
W  = Nx*h;

fprintf('Grid: Nx=%d, Nz=%d, h=%.4g m, W=%.4g m, D=%.4g m\n', Nx, Nz, h, W, D);

% --- Heat transport parameters ---
kT0  = 2;        % thermal conductivity [W/m/K]
rho0 = 2700;     % reference density [kg/m3]
cP0  = 1100;     % heat capacity [J/kg/K]
Qr0  = 0;        % volumetric heat source [W/m3]

% --- Buoyancy / EOS ---
g     = 9.81;
alpha = 3e-5;
Tref  = 0;

% --- Isothermal boundaries ---
Ttop = 0;
Tbot = 1000;

% --- Numerical controls ---
CFL_T   = 0.25;     % smaller is safer (avoid T overshoot)
CFL_P   = 0.2;
tol_div = 1e-8;
maxPIt  = 500;
maxStep = 2000;
plotEvery = 50;

yr = 365*24*3600;

% --- KD0 sensitivity sweep (as in practical slide) ---
KD0_list = logspace(-8, -6, 5);

results = struct([]);

for n = 1:numel(KD0_list)
    KD0 = KD0_list(n);

    fprintf('\n==== Running KD0 = %.3e ====\n', KD0);

    out = main_hydrothermal2D( ...
        W,D,Nx,Nz,h, ...
        kT0,rho0,cP0,Qr0, ...
        g,alpha,Tref,Ttop,Tbot, ...
        CFL_T,CFL_P,tol_div,maxPIt,maxStep,plotEvery,yr, ...
        KD0);

    results(n).KD0  = KD0;
    results(n).Nu   = out.Nu;
    results(n).vmax = out.vmax;
    results(n).t_yr = out.t/yr;

    fprintf('KD0=%.3e | Nu=%.4f | vmax=%.3e m/s | t=%.3g yr\n', ...
        KD0, out.Nu, out.vmax, out.t/yr);
end

figure('Name','Sensitivity Summary'); clf;
semilogx([results.KD0],[results.Nu],'-o'); grid on;
xlabel('KD0');
ylabel('Nusselt number');
title('Sensitivity to Darcy mobility parameter KD0');
