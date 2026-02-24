%% ========================================================================
% EARTH5016 - Helmsdale Geothermal Feasibility (Individual Final Report)
% File: project1_final_report_v2.m
%
% Purpose:
%   2D variable-coefficient heat transport model (diffusion + advection + source)
%   with report-ready diagnostics, data-pack-aligned geometry (right-side fault),
%   and improved figure layout / geometry tuning.
%
% Notes:
%   - Geometry revised to better match the project data pack conceptual section:
%       * He1 granite (left)
%       * He2 granite (central)
%       * right-side sedimentary wedge (Ms/Cm)
%       * right-side dipping fault
%       * basement gneiss at depth
%   - Advection remains a scenario-based fault-focused flow field
%     (not a full groundwater flow solver).
%   - y is positive downward.
%
% Author: [Your Name / Student Number]
% ========================================================================

clear; clc; close all;

%% =========================
% 1) MODEL CONFIGURATION
% ==========================
cfg = struct();

% --- Domain (0-16 km to align with data pack x-axis) ---
cfg.Lx = 16000;     % [m]
cfg.Ly = 6000;      % [m] numerical depth (figure focuses on upper 4.5 km)
cfg.Nx = 161;       % 100 m spacing
cfg.Ny = 61;        % 100 m spacing

cfg.x = linspace(0, cfg.Lx, cfg.Nx);
cfg.y = linspace(0, cfg.Ly, cfg.Ny);   % positive downward
[cfg.X, cfg.Y] = meshgrid(cfg.x, cfg.y);
cfg.dx = cfg.x(2) - cfg.x(1);
cfg.dy = cfg.y(2) - cfg.y(1);

% --- Thermal BC / IC ---
cfg.T_surf      = 10;      % [°C] top fixed temperature
cfg.grad_bottom = 0.035;   % [°C/m] bottom geothermal gradient
cfg.grad_init   = 0.030;   % [°C/m] initial linear gradient

% --- Fluid properties (advective heat capacity factor) ---
cfg.rho_f = 1000;          % [kg/m^3]
cfg.cp_f  = 4200;          % [J/kg/K]

% --- Time integration ---
cfg.time_scheme    = 'euler';  % 'euler' or 'rk2'
cfg.dt_safety      = 0.20;
cfg.runtime_factor = 1.30;     % runtime = factor * diffusive equilibration timescale

% --- Flow scenario (important for report physics discussion) ---
% 'upflow'   : Vy < 0
% 'downflow' : Vy > 0
% 'none'     : no advection
cfg.flow_mode = 'upflow';

% --- Output controls ---
cfg.save_figures  = true;
cfg.save_tables   = true;
cfg.out_prefix    = 'helmsdale_finalreport_v2';
cfg.show_progress = true;

%% =========================
% 2) PARAMETERISATION (BASELINE)
% ==========================
% Units:
%   k   [W m^-1 K^-1]
%   rho [kg m^-3]
%   cp  [J kg^-1 K^-1]
%   Qr  [W m^-3]
%   KD  [m^2 Pa^-1 s^-1] (proxy parameter in scenario advection)

par = struct();

% He1 granite
par.He1.k   = 2.67;
par.He1.rho = 2630;
par.He1.cp  = 836;
par.He1.Qr  = 6.53e-6;
par.He1.KD  = 1.0e-10;

% He2 granite
par.He2.k   = 2.67;
par.He2.rho = 2630;
par.He2.cp  = 836;
par.He2.Qr  = 5.20e-6;
par.He2.KD  = 1.0e-10;

% Sediments (Ms, Cm)
par.Ms.k   = 1.78;  par.Ms.rho = 2073;  par.Ms.cp = 1361;  par.Ms.Qr = 0.5e-6;  par.Ms.KD = 1.0e-8;
par.Cm.k   = 1.78;  par.Cm.rho = 2073;  par.Cm.cp = 1361;  par.Cm.Qr = 0.5e-6;  par.Cm.KD = 1.0e-8;

% Basement gneiss
par.Bg.k   = 2.50;
par.Bg.rho = 2700;
par.Bg.cp  = 900;
par.Bg.Qr  = 1.0e-6;
par.Bg.KD  = 1.0e-11;

% Fault zone (overprint)
par.Fz.k   = 2.70;
par.Fz.rho = 2299;
par.Fz.cp  = 1031;
par.Fz.Qr  = 1.0e-6;
par.Fz.KD  = 2.5e-7;   % key sensitivity parameter later

%% =========================
% 3) GEOLOGICAL GEOMETRY (DATA-PACK-ALIGNED CONCEPT)
% ==========================
% Type codes:
%   1 = He1
%   2 = He2
%   3 = Ms
%   4 = Cm
%   5 = Bg
%   6 = Fault
Type = zeros(cfg.Ny, cfg.Nx);

% --- Basement top (curved) ---
y_bg_top_1d = 4000 - 800 * ((cfg.x - 8000) / 8000).^2;   % [m]
y_bg_top_1d = max(3000, min(4200, y_bg_top_1d));
Ybg = repmat(y_bg_top_1d, cfg.Ny, 1);

% --- He1/He2 internal boundary (near-vertical, mild curvature) ---
% tuned to look closer to data pack (not overly curved)
y_ref = min(cfg.y, 4500);                      % shape mainly in mapped interval
x_he1_he2_1d = 5600 - 300 * (y_ref / 4500).^1.4;
Xhe = repmat(x_he1_he2_1d(:), 1, cfg.Nx);

% --- Fault geometry (right side; tuned) ---
fault_surf_x    = 11850;   % [m]
fault_dip_deg   = 55;      % from horizontal (larger = steeper)
fault_half_width = 140;    % [m], thinner visual fault zone
Xfault_center   = fault_surf_x + cfg.Y ./ tand(fault_dip_deg);

% --- Right-side sediment wedge (Ms over Cm), only on hanging-wall side ---
xw = (cfg.X - fault_surf_x) / (cfg.Lx - fault_surf_x);
xw = max(0, min(1, xw));

% tuned to better resemble data pack wedge thickness
y_ms_bottom = 550 + 450 * xw;    % ~0.55 km near fault to ~1.0 km far right
y_cm_bottom = 1300 + 700 * xw;   % ~1.3 km near fault to ~2.0 km far right

% --- Build masks ---
mask_bedrock_above_bg = (cfg.Y < Ybg);
mask_basement         = (cfg.Y >= Ybg);

mask_he1 = mask_bedrock_above_bg & (cfg.X < Xhe);
mask_he2 = mask_bedrock_above_bg & (cfg.X >= Xhe);

mask_right_of_fault = (cfg.X > (Xfault_center + fault_half_width));
mask_ms = mask_right_of_fault & (cfg.Y < y_ms_bottom) & mask_bedrock_above_bg;
mask_cm = mask_right_of_fault & (cfg.Y >= y_ms_bottom) & (cfg.Y < y_cm_bottom) & mask_bedrock_above_bg;

% Sediments replace granite where present
mask_he1(mask_ms | mask_cm) = false;
mask_he2(mask_ms | mask_cm) = false;

% Fault overprint
mask_fault = abs(cfg.X - Xfault_center) <= fault_half_width;

% Assign geology types
Type(mask_basement) = 5;
Type(mask_he1)      = 1;
Type(mask_he2)      = 2;
Type(mask_ms)       = 3;
Type(mask_cm)       = 4;
Type(mask_fault)    = 6;

%% =========================
% 4) BUILD PROPERTY FIELDS
% ==========================
K   = nan(cfg.Ny, cfg.Nx);
Rho = nan(cfg.Ny, cfg.Nx);
Cp  = nan(cfg.Ny, cfg.Nx);
Qr  = nan(cfg.Ny, cfg.Nx);
KD  = nan(cfg.Ny, cfg.Nx);

idx = (Type == 5); [K(idx), Rho(idx), Cp(idx), Qr(idx), KD(idx)] = deal(par.Bg.k,  par.Bg.rho,  par.Bg.cp,  par.Bg.Qr,  par.Bg.KD);
idx = (Type == 1); [K(idx), Rho(idx), Cp(idx), Qr(idx), KD(idx)] = deal(par.He1.k, par.He1.rho, par.He1.cp, par.He1.Qr, par.He1.KD);
idx = (Type == 2); [K(idx), Rho(idx), Cp(idx), Qr(idx), KD(idx)] = deal(par.He2.k, par.He2.rho, par.He2.cp, par.He2.Qr, par.He2.KD);
idx = (Type == 3); [K(idx), Rho(idx), Cp(idx), Qr(idx), KD(idx)] = deal(par.Ms.k,  par.Ms.rho,  par.Ms.cp,  par.Ms.Qr,  par.Ms.KD);
idx = (Type == 4); [K(idx), Rho(idx), Cp(idx), Qr(idx), KD(idx)] = deal(par.Cm.k,  par.Cm.rho,  par.Cm.cp,  par.Cm.Qr,  par.Cm.KD);
idx = (Type == 6); [K(idx), Rho(idx), Cp(idx), Qr(idx), KD(idx)] = deal(par.Fz.k,  par.Fz.rho,  par.Fz.cp,  par.Fz.Qr,  par.Fz.KD);

if any(isnan(K(:))) || any(isnan(Rho(:))) || any(isnan(Cp(:))) || any(isnan(Qr(:))) || any(isnan(KD(:)))
    error('Property assignment incomplete: some cells remain NaN.');
end

%% =========================
% 5) FLOW FIELD (SCENARIO-BASED ADVECTION)
% ==========================
% y positive downward:
%   Vy < 0 => upward flow (upflow)
%   Vy > 0 => downward flow (downflow)

Vx = zeros(cfg.Ny, cfg.Nx);
Vy = zeros(cfg.Ny, cfg.Nx);

flow_scale = 2.0e-2;  % chosen to give fault Vy ~1e-9 to 1e-8 m/s for KD ~1e-7-1e-8

switch lower(cfg.flow_mode)
    case 'upflow'
        Vy(mask_fault) = -flow_scale * KD(mask_fault);
    case 'downflow'
        Vy(mask_fault) = +flow_scale * KD(mask_fault);
    case 'none'
        Vy(:) = 0;
    otherwise
        error('Unknown cfg.flow_mode. Use upflow/downflow/none.');
end

%% =========================
% 6) NUMERICAL SETUP (VERIFICATION-ORIENTED)
% ==========================
D = K ./ (Rho .* Cp);                     % thermal diffusivity [m^2/s]
AdvCoeff = (cfg.rho_f * cfg.cp_f) ./ (Rho .* Cp);

% Runtime from diffusive equilibration timescale
kappa_rep   = median(D(:));
t_diff_sec  = cfg.Ly^2 / (pi^2 * kappa_rep);
sec_per_yr  = 365.25 * 24 * 3600;
t_diff_yr   = t_diff_sec / sec_per_yr;
runtime_yr  = cfg.runtime_factor * t_diff_yr;

% Explicit stability limits
Dmax = max(D(:));
dt_diff = cfg.dt_safety / (2 * Dmax * (1/cfg.dx^2 + 1/cfg.dy^2));

vx_max = max(abs(Vx(:)));
vy_max = max(abs(Vy(:)));

dt_adv = inf;
if vx_max > 0
    dt_adv = min(dt_adv, cfg.dt_safety * cfg.dx / vx_max);
end
if vy_max > 0
    dt_adv = min(dt_adv, cfg.dt_safety * cfg.dy / vy_max);
end

dt = min(dt_diff, dt_adv);
nt = ceil(runtime_yr * sec_per_yr / dt);

% Initial condition
T = cfg.T_surf + cfg.grad_init * cfg.Y;
T = apply_temperature_BC(T, cfg);

%% =========================
% 7) SOLVER LOOP (EXPLICIT)
% ==========================
fprintf('\n=== Helmsdale 2D Heat Transport Model (Final Report v2) ===\n');
fprintf('Domain: %.1f km x %.1f km | Grid: %d x %d | dx=%.0f m, dy=%.0f m\n', ...
    cfg.Lx/1000, cfg.Ly/1000, cfg.Nx, cfg.Ny, cfg.dx, cfg.dy);
fprintf('Flow mode: %s | Time scheme: %s (explicit)\n', cfg.flow_mode, upper(cfg.time_scheme));
fprintf('kappa_rep = %.3e m^2/s\n', kappa_rep);
fprintf('Diffusive timescale = %.1f kyr | Runtime factor = %.2f | Runtime = %.1f kyr\n', ...
    t_diff_yr/1e3, cfg.runtime_factor, runtime_yr/1e3);
fprintf('dt_diff = %.3e s | dt_adv = %.3e s | dt = %.3e s | nt = %d\n', dt_diff, dt_adv, dt, nt);
fprintf('Max |Vy| in fault = %.3e m/s\n', max(abs(Vy(:))));

for n = 1:nt
    switch lower(cfg.time_scheme)
        case 'euler'
            RHS = compute_rhs_temperature(T, K, Rho, Cp, Qr, Vx, Vy, AdvCoeff, cfg);
            Tnew = T + dt * RHS;
            Tnew = apply_temperature_BC(Tnew, cfg);

        case 'rk2'
            k1 = compute_rhs_temperature(T, K, Rho, Cp, Qr, Vx, Vy, AdvCoeff, cfg);
            Tstar = apply_temperature_BC(T + dt * k1, cfg);

            k2 = compute_rhs_temperature(Tstar, K, Rho, Cp, Qr, Vx, Vy, AdvCoeff, cfg);
            Tnew = apply_temperature_BC(T + 0.5 * dt * (k1 + k2), cfg);

        otherwise
            error('Unsupported time scheme.');
    end

    T = Tnew;

    if cfg.show_progress
        stride = max(1, floor(nt/10));
        if mod(n, stride) == 0 || n == nt
            fprintf('Progress: %3.0f%%\n', 100*n/nt);
        end
    end
end

%% =========================
% 8) VERIFICATION & VALIDATION DIAGNOSTICS
% ==========================
verify = struct();
verify.scheme       = cfg.time_scheme;
verify.dt           = dt;
verify.dt_diff      = dt_diff;
verify.dt_adv       = dt_adv;
verify.kappa_rep    = kappa_rep;
verify.t_diff_kyr   = t_diff_yr/1e3;
verify.runtime_kyr  = runtime_yr/1e3;
verify.dx_m         = cfg.dx;
verify.dy_m         = cfg.dy;

% He1 observational tie-in (partial validation)
He1_depth_m = [51, 150, 249, 351, 450, 550, 649, 750, 850]';
He1_temp_C  = [15.1, 15.9, 21.5, 21.0, 25.0, 28.3, 31.6, 34.2, 38.5]';

% Representative profiles (aligned with revised geometry)
idx_he1   = round(3000  / cfg.dx) + 1;
idx_he2   = round(9000  / cfg.dx) + 1;
idx_fault = round(13000 / cfg.dx) + 1;
idx_right = round(15000 / cfg.dx) + 1;

T_model_he1 = interp1(cfg.y(:), T(:, idx_he1), He1_depth_m, 'linear', 'extrap');
rmse_he1 = sqrt(mean((T_model_he1 - He1_temp_C).^2));

validate = struct();
validate.He1_RMSE_C = rmse_he1;
validate.note = 'He1 is a shallow observational tie-in only; deep validation remains limited.';

fprintf('\n--- Verification diagnostics (numerical) ---\n');
fprintf('Scheme: %s | dt = %.3e s\n', upper(verify.scheme), verify.dt);
fprintf('dt_diff = %.3e s | dt_adv = %.3e s\n', verify.dt_diff, verify.dt_adv);
fprintf('Diffusive timescale = %.1f kyr | Runtime = %.1f kyr\n', verify.t_diff_kyr, verify.runtime_kyr);

fprintf('\n--- Validation tie-in (partial) ---\n');
fprintf('He1 shallow-profile RMSE = %.2f °C\n', validate.He1_RMSE_C);

%% =========================
% 9) REPORT METRICS: ISOTHERM DEPTHS
% ==========================
targets = [50 70 100 120];

d_he1   = estimate_isotherm_depths(T(:, idx_he1),   cfg.y, targets);
d_he2   = estimate_isotherm_depths(T(:, idx_he2),   cfg.y, targets);
d_fault = estimate_isotherm_depths(T(:, idx_fault), cfg.y, targets);
d_right = estimate_isotherm_depths(T(:, idx_right), cfg.y, targets);

tbl_isotherms = table(targets(:), d_he1(:)/1000, d_he2(:)/1000, d_fault(:)/1000, d_right(:)/1000, ...
    'VariableNames', {'Target_C','He1_km','He2_km','FaultCorridor_km','RightWedge_km'});

fprintf('\n--- Preliminary isotherm depths (km) ---\n');
disp(tbl_isotherms);

if cfg.save_tables
    writetable(tbl_isotherms, [cfg.out_prefix '_isotherm_depths.csv']);
end

%% =========================
% 10) VISUALISATION (IMPROVED LAYOUT)
% ==========================

% ---------- Figure 1: Thermal structure + profiles (tiledlayout) ----------
f1 = figure('Color','w', 'Position',[60 60 1400 700]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% --- Left panel: 2D thermal structure ---
ax1 = nexttile(tl,1); hold(ax1,'on');
imagesc(ax1, cfg.x/1000, cfg.y/1000, T);
set(ax1,'YDir','reverse');
xlim(ax1,[0 cfg.Lx/1000]);
ylim(ax1,[0 4.5]);                  % report-friendly view
pbaspect(ax1,[16 4.5 1]);           % fixes "squashed" plot issue
box(ax1,'on');
colormap(ax1, jet);
caxis(ax1,[10 160]);
cb = colorbar(ax1);
ylabel(cb, 'Temperature (°C)');

% Geology overlays
plot(ax1, cfg.x/1000, y_bg_top_1d/1000, 'k-', 'LineWidth',1.2);  % basement top
plot(ax1, x_he1_he2_1d/1000, cfg.y/1000, 'k--', 'LineWidth',1.0); % He1/He2 boundary

xfault = (fault_surf_x + cfg.y./tand(fault_dip_deg))/1000;
plot(ax1, xfault, cfg.y/1000, 'r-', 'LineWidth',2.2);                          % fault center
plot(ax1, xfault - fault_half_width/1000, cfg.y/1000, 'r--', 'LineWidth',1.0); % fault edges
plot(ax1, xfault + fault_half_width/1000, cfg.y/1000, 'r--', 'LineWidth',1.0);

% Sedimentary layer boundaries (defined across full x for convenience)
plot(ax1, cfg.x/1000, y_ms_bottom(1,:)/1000, 'c-', 'LineWidth',1.0);
plot(ax1, cfg.x/1000, y_cm_bottom(1,:)/1000, 'm-', 'LineWidth',1.0);

% Isotherms
[Ciso, hIso] = contour(ax1, cfg.x/1000, cfg.y/1000, T, targets, 'w-', 'LineWidth',1.3);
clabel(Ciso, hIso, 'Color','w', 'FontSize',9, 'FontWeight','bold');

% Labels
text(ax1, 1.5, 1.2, 'He1 granite', 'FontWeight','bold', 'Color','k');
text(ax1, 8.0, 1.2, 'He2 granite', 'FontWeight','bold', 'Color','k');
text(ax1, 14.2, 0.9, 'Ms', 'FontWeight','bold', 'Color','k');
text(ax1, 14.2, 1.8, 'Cm', 'FontWeight','bold', 'Color','k');
text(ax1, 13.8, 3.9, 'Basement', 'FontWeight','bold', 'Color','k');
text(ax1, 12.1, 3.2, sprintf('Fault (%s)', cfg.flow_mode), ...
    'Color','r', 'FontWeight','bold', 'Rotation', -55);

xlabel(ax1,'Distance (km)');
ylabel(ax1,'Depth (km)');
title(ax1,'A. Preliminary 2D Thermal Structure (Report Geometry)');

% --- Right panel: profiles + He1 tie-in ---
ax2 = nexttile(tl,2); hold(ax2,'on');
plot(ax2, T(:, idx_he1),   cfg.y/1000, 'r-',  'LineWidth',2, 'DisplayName','Model: He1 profile');
plot(ax2, T(:, idx_he2),   cfg.y/1000, 'b--', 'LineWidth',2, 'DisplayName','Model: He2 profile');
plot(ax2, T(:, idx_fault), cfg.y/1000, 'g-.', 'LineWidth',2, 'DisplayName','Model: Fault corridor');
plot(ax2, T(:, idx_right), cfg.y/1000, 'k:',  'LineWidth',2, 'DisplayName','Model: Right wedge');
plot(ax2, He1_temp_C, He1_depth_m/1000, 'ko', 'MarkerFaceColor','y', 'MarkerSize',7, ...
    'DisplayName','Data: Drill Hole He1');

set(ax2,'YDir','reverse');
ylim(ax2,[0 4.5]);
grid(ax2,'on');
xlabel(ax2,'Temperature (°C)');
ylabel(ax2,'Depth (km)');
title(ax2, sprintf('B. Profile Comparison and He1 Tie-in (RMSE = %.2f °C)', rmse_he1));
legend(ax2,'Location','southwest');

title(tl, 'Helmsdale Geothermal Feasibility Study - Final Report Model (WIP)', ...
    'FontWeight','bold','FontSize',14);

if cfg.save_figures
    exportgraphics(f1, [cfg.out_prefix '_thermal_and_profiles.png'], 'Resolution', 300);
end

% ---------- Figure 2: Geology zoning check ----------
f2 = figure('Color','w', 'Position',[120 120 1100 460]); hold on;
imagesc(cfg.x/1000, cfg.y/1000, Type);
set(gca,'YDir','reverse');
axis tight; box on;
xlim([0 cfg.Lx/1000]);
ylim([0 4.5]);
pbaspect([16 4.5 1]);
xlabel('Distance (km)');
ylabel('Depth (km)');
title('Geological Zoning Used in the Numerical Model (Geometry Check)');

% Geology colormap
cmap = [
    0.95 0.55 0.50;  % He1
    0.86 0.68 0.80;  % He2
    0.70 0.88 0.95;  % Ms
    0.55 0.82 0.70;  % Cm
    0.80 0.73 0.65;  % Bg
    0.90 0.10 0.10   % Fault
];
colormap(gca, cmap);
caxis([1 6]);
cb2 = colorbar('Ticks',1:6, 'TickLabels',{'He1','He2','Ms','Cm','Bg','Fault'});
ylabel(cb2,'Geology Type');

% overlays
plot(cfg.x/1000, y_bg_top_1d/1000, 'k-', 'LineWidth',1.2);
plot(x_he1_he2_1d/1000, cfg.y/1000, 'k--', 'LineWidth',1.0);
plot(xfault, cfg.y/1000, 'w-', 'LineWidth',2.0);

if cfg.save_figures
    exportgraphics(f2, [cfg.out_prefix '_geology_zoning_check.png'], 'Resolution', 300);
end

%% =========================
% 11) REPORT NOTES (CONSOLE)
% ==========================
fprintf('\nReport-writing reminders:\n');
fprintf('1) Verification = numerical correctness (dt limits, scheme, BC implementation).\n');
fprintf('2) Validation  = real-system appropriateness (geometry consistency + He1 tie-in).\n');
fprintf('3) Advection can enhance OR suppress shallow temperatures depending on flow direction.\n');
fprintf('4) Current advection is a scenario field, not a full groundwater flow solution.\n');
fprintf('5) For the main report parameter study, vary granite Qr and fault KD.\n');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function Tbc = apply_temperature_BC(Tin, cfg)
% Enforce boundary conditions on cell-centre values:
%   Top    : fixed T = T_surf
%   Bottom : fixed gradient dT/dy = grad_bottom
%   Left   : insulated (dT/dx = 0)
%   Right  : insulated (dT/dx = 0)

    Tbc = Tin;

    % Top Dirichlet
    Tbc(1,:) = cfg.T_surf;

    % Bottom Neumann
    Tbc(end,:) = Tbc(end-1,:) + cfg.grad_bottom * cfg.dy;

    % Left/right insulated
    Tbc(:,1)   = Tbc(:,2);
    Tbc(:,end) = Tbc(:,end-1);
end

function RHS = compute_rhs_temperature(T, K, Rho, Cp, Qr, Vx, Vy, AdvCoeff, cfg)
% RHS of:
% dT/dt = (1/rho/cp) * div(k grad T) - AdvCoeff*(v·gradT) + Qr/(rho cp)
%
% Boundary-consistent finite differences via ghost cells.

    % Enforce BCs before stencil
    T = apply_temperature_BC(T, cfg);

    % Ghost-padded fields
    [Tpad, Kpad] = build_ghost_fields(T, K, cfg);

    % Interior aliases
    Tc = Tpad(2:end-1, 2:end-1);
    Te = Tpad(2:end-1, 3:end);
    Tw = Tpad(2:end-1, 1:end-2);
    Tn = Tpad(1:end-2, 2:end-1);   % shallower
    Ts = Tpad(3:end,   2:end-1);   % deeper

    Kc = Kpad(2:end-1, 2:end-1);
    Ke = Kpad(2:end-1, 3:end);
    Kw = Kpad(2:end-1, 1:end-2);
    Kn = Kpad(1:end-2, 2:end-1);
    Ks = Kpad(3:end,   2:end-1);

    % ---- Diffusion (conservative variable-coefficient form) ----
    k_e = 0.5 * (Kc + Ke);
    k_w = 0.5 * (Kc + Kw);
    k_n = 0.5 * (Kc + Kn);
    k_s = 0.5 * (Kc + Ks);

    div_k_gradT = ...
        (k_e .* (Te - Tc) - k_w .* (Tc - Tw)) / cfg.dx^2 + ...
        (k_s .* (Ts - Tc) - k_n .* (Tc - Tn)) / cfg.dy^2;

    DiffTerm = div_k_gradT ./ (Rho .* Cp);

    % ---- Advection (upwind) ----
    dTdx = zeros(size(T));
    dTdy = zeros(size(T));

    % x-direction upwind
    mask_vx_pos = (Vx > 0);
    mask_vx_neg = (Vx < 0);
    dTdx(mask_vx_pos) = (Tc(mask_vx_pos) - Tw(mask_vx_pos)) / cfg.dx;
    dTdx(mask_vx_neg) = (Te(mask_vx_neg) - Tc(mask_vx_neg)) / cfg.dx;

    % y-direction upwind (y positive downward)
    mask_vy_pos = (Vy > 0);   % downward flow: upwind from Tn
    mask_vy_neg = (Vy < 0);   % upward flow: upwind from Ts
    dTdy(mask_vy_pos) = (Tc(mask_vy_pos) - Tn(mask_vy_pos)) / cfg.dy;
    dTdy(mask_vy_neg) = (Ts(mask_vy_neg) - Tc(mask_vy_neg)) / cfg.dy;

    AdvTerm = AdvCoeff .* (Vx .* dTdx + Vy .* dTdy);

    % ---- Source ----
    SourceTerm = Qr ./ (Rho .* Cp);

    RHS = DiffTerm - AdvTerm + SourceTerm;
end

function [Tpad, Kpad] = build_ghost_fields(T, K, cfg)
% Ghost cells for:
%   top Dirichlet, bottom Neumann, left/right Neumann

    [Ny, Nx] = size(T);

    Tpad = zeros(Ny+2, Nx+2);
    Kpad = zeros(Ny+2, Nx+2);

    % Interior
    Tpad(2:end-1, 2:end-1) = T;
    Kpad(2:end-1, 2:end-1) = K;

    % Left/right ghosts (insulated)
    Tpad(2:end-1, 1)   = T(:,1);
    Tpad(2:end-1, end) = T(:,end);
    Kpad(2:end-1, 1)   = K(:,1);
    Kpad(2:end-1, end) = K(:,end);

    % Top ghost (Dirichlet)
    Tpad(1, 2:end-1) = 2*cfg.T_surf - T(1,:);
    Kpad(1, 2:end-1) = K(1,:);

    % Bottom ghost (Neumann)
    Tpad(end, 2:end-1) = T(end,:) + cfg.grad_bottom * cfg.dy;
    Kpad(end, 2:end-1) = K(end,:);

    % Corners
    Tpad(1,1)     = Tpad(1,2);
    Tpad(1,end)   = Tpad(1,end-1);
    Tpad(end,1)   = Tpad(end,2);
    Tpad(end,end) = Tpad(end,end-1);

    Kpad(1,1)     = Kpad(2,2);
    Kpad(1,end)   = Kpad(2,end-1);
    Kpad(end,1)   = Kpad(end-1,2);
    Kpad(end,end) = Kpad(end-1,end-1);
end

function depths_m = estimate_isotherm_depths(Tcol, y, targets)
% Linear interpolation of isotherm depths.
% Returns NaN if target not reached.

    depths_m = nan(size(targets));

    for i = 1:numel(targets)
        Ti = targets(i);
        idx = find(Tcol >= Ti, 1, 'first');

        if isempty(idx)
            depths_m(i) = NaN;
            continue;
        end

        if idx == 1
            depths_m(i) = y(1);
            continue;
        end

        y1 = y(idx-1); y2 = y(idx);
        T1 = Tcol(idx-1); T2 = Tcol(idx);

        if abs(T2 - T1) < eps
            depths_m(i) = y2;
        else
            depths_m(i) = y1 + (Ti - T1) * (y2 - y1) / (T2 - T1);
        end
    end
end