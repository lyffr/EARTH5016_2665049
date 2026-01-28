% run_advdiff2D.m
%***** RUN 2D ADVECTION DIFFUSION MODEL  *********************************
clear; close all;

%--- Model Parameters (Slide 48) ---
W     = 1e3;           % domain width x-direction [m]
D     = 1e3;           % domain depth z-direction [m]
Nz    = 100;           % grid size z-direction
Nx    = Nz * (W/D);    % grid size x-direction (ensure square cells)
h     = D/Nz;          % grid spacing (h = dx = dz) [m]

%--- Physics Parameters ---
kT0   = 2;             % thermal conductivity [W/m/K]
rho0  = 2700;          % density [kg/m3]
cP0   = 1100;          % heat capacity [J/kg/K]
Qr0   = 1e-6;          % heat productivity [W/m3]
u0    = 1e-6;          % advection x-speed [m/s]
w0    = 1e-6;          % advection z-speed [m/s]

%--- Initial Condition Parameters ---
T0    = 100;           % background temperature
dT    = 1000;          % peak temperature anomaly
sgm0  = 50;            % initial pulse width [m]

%--- Numerical Controls ---
BC    = 'periodic';    % boundary condition
ADVN  = 'UPW3';        % advection scheme
TINT  = 'RK2';         % time integration
CFL   = 0.5;           % stability limiter
tend  = W/u0/2;        % stop time
nop   = 50;            % plot every 'nop' steps
yr    = 365*24*3600;   % seconds in a year

%--- Run Solver ---
run('main_solution2D.m');
