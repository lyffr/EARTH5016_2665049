%% ========================================================================
% EARTH5016 - Helmsdale Geothermal Feasibility
% File: run_parameter_study.m
%
% Purpose:
%   3x3 parameter study for the individual report:
%   - Parameter A: He2 granite radiogenic heat production (Qr)
%   - Parameter B: Fault-zone KD (advective transport proxy)
%
% Outputs:
%   - CSV table with scenario metrics
%   - Sensitivity plots for 100C and 120C isotherm depths (fault corridor)
%
% Notes:
%   - Standalone script (does NOT require project1_final_report_v2.m)
%   - Reuses the same numerical model structure and geometry logic
% ========================================================================

clear; clc; close all;

%% =========================
% 0) STUDY SETTINGS
% ==========================
study = struct();

% --- Which granite Qr to vary ---
% Recommended: vary He2 only (He1 is tied to shallow observations)
study.vary_which_Qr = 'He2';   % 'He2' or 'BothHe1He2'

% --- 3x3 parameter grid ---
study.He2_Qr_factors = [0.8, 1.0, 1.2];     % multiplier on baseline He2.Qr
study.Fz_KD_values   = [1.0e-7, 2.5e-7, 5.0e-7];  % [m^2 Pa^-1 s^-1]

% --- Flow scenario for the study ---
study.flow_mode = 'upflow';   % 'upflow', 'downflow', 'none'

% --- Save outputs ---
study.out_prefix = 'helmsdale_paramstudy';
study.save_csv   = true;
study.save_figs  = true;

%% =========================
% 1) BASE MODEL CONFIGURATION (same style as final_report_v2)
% ==========================
cfg = struct();

cfg.Lx = 16000;     % [m]
cfg.Ly = 6000;      % [m]
cfg.Nx = 161;       % 100 m spacing
cfg.Ny = 61;        % 100 m spacing

cfg.x = linspace(0, cfg.Lx, cfg.Nx);
cfg.y = linspace(0, cfg.Ly, cfg.Ny);
[cfg.X, cfg.Y] = meshgrid(cfg.x, cfg.y);
cfg.dx = cfg.x(2) - cfg.x(1);
cfg.dy = cfg.y(2) - cfg.y(1);

cfg.T_surf      = 10;      % [°C]
cfg.grad_bottom = 0.035;   % [°C/m]
cfg.grad_init   = 0.030;   % [°C/m]

cfg.rho_f = 1000;          % [kg/m^3]
cfg.cp_f  = 4200;          % [J/kg/K]

cfg.time_scheme    = 'euler';   % 'euler' or 'rk2'
cfg.dt_safety      = 0.20;
cfg.runtime_factor = 1.30;

cfg.flow_mode = study.flow_mode;
cfg.show_progress = false;      % keep off for batch runs

%% =========================
% 2) BASELINE PARAMETERS
% ==========================
par0 = struct();

% He1 granite
par0.He1.k   = 2.67;
par0.He1.rho = 2630;
par0.He1.cp  = 836;
par0.He1.Qr  = 6.53e-6;
par0.He1.KD  = 1.0e-10;

% He2 granite
par0.He2.k   = 2.67;
par0.He2.rho = 2630;
par0.He2.cp  = 836;
par0.He2.Qr  = 5.20e-6;
par0.He2.KD  = 1.0e-10;

% Sediments
par0.Ms.k   = 1.78;  par0.Ms.rho = 2073;  par0.Ms.cp = 1361;  par0.Ms.Qr = 0.5e-6;  par0.Ms.KD = 1.0e-8;
par0.Cm.k   = 1.78;  par0.Cm.rho = 2073;  par0.Cm.cp = 1361;  par0.Cm.Qr = 0.5e-6;  par0.Cm.KD = 1.0e-8;

% Basement
par0.Bg.k   = 2.50;
par0.Bg.rho = 2700;
par0.Bg.cp  = 900;
par0.Bg.Qr  = 1.0e-6;
par0.Bg.KD  = 1.0e-11;

% Fault zone
par0.Fz.k   = 2.70;
par0.Fz.rho = 2299;
par0.Fz.cp  = 1031;
par0.Fz.Qr  = 1.0e-6;
par0.Fz.KD  = 2.5e-7;

%% =========================
% 3) RUN 3x3 PARAMETER STUDY
% ==========================
nQ = numel(study.He2_Qr_factors);
nK = numel(study.Fz_KD_values);
nCases = nQ * nK;

results(nCases,1) = struct();
case_id = 0;

fprintf('\n=== Running 3x3 parameter study (%s, flow=%s) ===\n', study.vary_which_Qr, cfg.flow_mode);

for iQ = 1:nQ
    for iK = 1:nK
        case_id = case_id + 1;

        par = par0;

        % Apply Qr variation
        qf = study.He2_Qr_factors(iQ);
        switch lower(study.vary_which_Qr)
            case 'he2'
                par.He2.Qr = par0.He2.Qr * qf;
            case 'bothhe1he2'
                par.He1.Qr = par0.He1.Qr * qf;
                par.He2.Qr = par0.He2.Qr * qf;
            otherwise
                error('Unknown study.vary_which_Qr option.');
        end

        % Apply fault KD variation
        par.Fz.KD = study.Fz_KD_values(iK);

        fprintf('Case %d/%d: Qr factor = %.2f, Fz.KD = %.2e\n', case_id, nCases, qf, par.Fz.KD);

        % Run model
        out = run_helmsdale_case(cfg, par);

        % Store flat results
        results(case_id).CaseID = case_id;
        results(case_id).FlowMode = cfg.flow_mode;
        results(case_id).QrFactor = qf;
        results(case_id).He1_Qr = par.He1.Qr;
        results(case_id).He2_Qr = par.He2.Qr;
        results(case_id).Fz_KD = par.Fz.KD;

        results(case_id).He1_RMSE_C = out.rmse_he1;
        results(case_id).dt_s = out.verify.dt;
        results(case_id).dt_diff_s = out.verify.dt_diff;
        results(case_id).dt_adv_s = out.verify.dt_adv;
        results(case_id).runtime_kyr = out.verify.runtime_kyr;

        % Isotherm depths [km]
        results(case_id).He1_50_km   = out.iso.He1_km(1);
        results(case_id).He1_70_km   = out.iso.He1_km(2);
        results(case_id).He1_100_km  = out.iso.He1_km(3);
        results(case_id).He1_120_km  = out.iso.He1_km(4);

        results(case_id).He2_50_km   = out.iso.He2_km(1);
        results(case_id).He2_70_km   = out.iso.He2_km(2);
        results(case_id).He2_100_km  = out.iso.He2_km(3);
        results(case_id).He2_120_km  = out.iso.He2_km(4);

        results(case_id).Fault_50_km  = out.iso.Fault_km(1);
        results(case_id).Fault_70_km  = out.iso.Fault_km(2);
        results(case_id).Fault_100_km = out.iso.Fault_km(3);
        results(case_id).Fault_120_km = out.iso.Fault_km(4);

        results(case_id).Right_50_km  = out.iso.Right_km(1);
        results(case_id).Right_70_km  = out.iso.Right_km(2);
        results(case_id).Right_100_km = out.iso.Right_km(3);
        results(case_id).Right_120_km = out.iso.Right_km(4);

        results(case_id).MaxT_C = out.maxT;
    end
end

%% =========================
% 4) RESULTS TABLE + CSV
% ==========================
TBL = struct2table(results);

% Sort for neat output
TBL = sortrows(TBL, {'QrFactor','Fz_KD'});

disp(' ');
disp('=== Parameter study summary table ===');
disp(TBL(:, {'CaseID','QrFactor','Fz_KD','He1_RMSE_C','Fault_100_km','Fault_120_km','He2_100_km','Right_100_km'}));

if study.save_csv
    writetable(TBL, [study.out_prefix '_summary.csv']);
end

%% =========================
% 5) SENSITIVITY PLOTS (for report)
% ==========================
% Build matrices for plotting (rows = Qr factor, cols = Fz KD)
M_Fault100 = nan(nQ, nK);
M_Fault120 = nan(nQ, nK);
M_He1RMSE  = nan(nQ, nK);

for iQ = 1:nQ
    for iK = 1:nK
        row = TBL(TBL.QrFactor == study.He2_Qr_factors(iQ) & TBL.Fz_KD == study.Fz_KD_values(iK), :);
        M_Fault100(iQ,iK) = row.Fault_100_km;
        M_Fault120(iQ,iK) = row.Fault_120_km;
        M_He1RMSE(iQ,iK)  = row.He1_RMSE_C;
    end
end

% --- Plot 1: Fault 100C isotherm depth vs Fz KD ---
f1 = figure('Color','w','Position',[100 100 820 560]); hold on;
for iQ = 1:nQ
    semilogx(study.Fz_KD_values, M_Fault100(iQ,:), '-o', 'LineWidth', 2, ...
        'DisplayName', sprintf('He2 Q_r factor = %.1f', study.He2_Qr_factors(iQ)));
end
grid on;
xlabel('Fault-zone KD (m^2 Pa^{-1} s^{-1})');
ylabel('100°C isotherm depth in fault corridor (km)');
title('Sensitivity of 100°C Isotherm Depth to Fault KD and He2 Q_r');
legend('Location','best');

if study.save_figs
    exportgraphics(f1, [study.out_prefix '_fault100_depth_sensitivity.png'], 'Resolution', 300);
end

% --- Plot 2: Fault 120C isotherm depth vs Fz KD ---
f2 = figure('Color','w','Position',[120 120 820 560]); hold on;
for iQ = 1:nQ
    semilogx(study.Fz_KD_values, M_Fault120(iQ,:), '-s', 'LineWidth', 2, ...
        'DisplayName', sprintf('He2 Q_r factor = %.1f', study.He2_Qr_factors(iQ)));
end
grid on;
xlabel('Fault-zone KD (m^2 Pa^{-1} s^{-1})');
ylabel('120°C isotherm depth in fault corridor (km)');
title('Sensitivity of 120°C Isotherm Depth to Fault KD and He2 Q_r');
legend('Location','best');

if study.save_figs
    exportgraphics(f2, [study.out_prefix '_fault120_depth_sensitivity.png'], 'Resolution', 300);
end

% --- Plot 3: He1 RMSE check (to track calibration consistency) ---
f3 = figure('Color','w','Position',[140 140 820 560]); hold on;
for iQ = 1:nQ
    semilogx(study.Fz_KD_values, M_He1RMSE(iQ,:), '-d', 'LineWidth', 2, ...
        'DisplayName', sprintf('He2 Q_r factor = %.1f', study.He2_Qr_factors(iQ)));
end
grid on;
xlabel('Fault-zone KD (m^2 Pa^{-1} s^{-1})');
ylabel('He1 RMSE (°C)');
title('He1 RMSE across Parameter Study (Validation Tie-in Check)');
legend('Location','best');

if study.save_figs
    exportgraphics(f3, [study.out_prefix '_He1_RMSE_check.png'], 'Resolution', 300);
end

%% =========================
% 6) CONSOLE SUMMARY (best / practical cases)
% ==========================
% Example ranking: prioritize shallow 120C in fault corridor, but avoid poor He1 RMSE
validRows = TBL.He1_RMSE_C <= min(TBL.He1_RMSE_C) + 2.0;   % loose filter
TBL_valid = TBL(validRows,:);

if ~isempty(TBL_valid)
    TBL_rank = sortrows(TBL_valid, {'Fault_120_km','Fault_100_km','He1_RMSE_C'}, {'ascend','ascend','ascend'});
    disp(' ');
    disp('=== Top candidate scenarios (shallow hot fault corridor, with acceptable He1 RMSE) ===');
    disp(TBL_rank(1:min(5,height(TBL_rank)), {'CaseID','QrFactor','Fz_KD','He1_RMSE_C','Fault_100_km','Fault_120_km'}));
end

fprintf('\nDone. Outputs saved with prefix: %s\n', study.out_prefix);

%% ========================================================================
% LOCAL FUNCTION: RUN ONE CASE
% ========================================================================
function out = run_helmsdale_case(cfg, par)
% Runs one model realization and returns key metrics only (no heavy plotting)

    % -------------------------
    % 1) GEOLOGICAL GEOMETRY (same logic as final_report_v2)
    % -------------------------
    Type = zeros(cfg.Ny, cfg.Nx);

    % Basement top
    y_bg_top_1d = 4000 - 800 * ((cfg.x - 8000) / 8000).^2;
    y_bg_top_1d = max(3000, min(4200, y_bg_top_1d));
    Ybg = repmat(y_bg_top_1d, cfg.Ny, 1);

    % He1/He2 boundary (mild curvature)
    y_ref = min(cfg.y, 4500);
    x_he1_he2_1d = 5600 - 300 * (y_ref / 4500).^1.4;
    Xhe = repmat(x_he1_he2_1d(:), 1, cfg.Nx);

    % Fault geometry (right side)
    fault_surf_x = 11850;
    fault_dip_deg = 55;
    fault_half_width = 140;
    Xfault_center = fault_surf_x + cfg.Y ./ tand(fault_dip_deg);

    % Sediment wedge
    xw = (cfg.X - fault_surf_x) / (cfg.Lx - fault_surf_x);
    xw = max(0, min(1, xw));
    y_ms_bottom = 550 + 450 * xw;
    y_cm_bottom = 1300 + 700 * xw;

    % Masks
    mask_bedrock_above_bg = (cfg.Y < Ybg);
    mask_basement         = (cfg.Y >= Ybg);

    mask_he1 = mask_bedrock_above_bg & (cfg.X < Xhe);
    mask_he2 = mask_bedrock_above_bg & (cfg.X >= Xhe);

    mask_right_of_fault = (cfg.X > (Xfault_center + fault_half_width));
    mask_ms = mask_right_of_fault & (cfg.Y < y_ms_bottom) & mask_bedrock_above_bg;
    mask_cm = mask_right_of_fault & (cfg.Y >= y_ms_bottom) & (cfg.Y < y_cm_bottom) & mask_bedrock_above_bg;

    mask_he1(mask_ms | mask_cm) = false;
    mask_he2(mask_ms | mask_cm) = false;

    mask_fault = abs(cfg.X - Xfault_center) <= fault_half_width;

    % Assign types
    Type(mask_basement) = 5;
    Type(mask_he1)      = 1;
    Type(mask_he2)      = 2;
    Type(mask_ms)       = 3;
    Type(mask_cm)       = 4;
    Type(mask_fault)    = 6;

    % -------------------------
    % 2) PROPERTY FIELDS
    % -------------------------
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

    if any(isnan(K(:)))
        error('NaN in property assignment.');
    end

    % -------------------------
    % 3) FLOW FIELD (scenario-based)
    % -------------------------
    Vx = zeros(cfg.Ny, cfg.Nx);
    Vy = zeros(cfg.Ny, cfg.Nx);

    flow_scale = 2.0e-2;
    switch lower(cfg.flow_mode)
        case 'upflow'
            Vy(mask_fault) = -flow_scale * KD(mask_fault);
        case 'downflow'
            Vy(mask_fault) = +flow_scale * KD(mask_fault);
        case 'none'
            Vy(:) = 0;
        otherwise
            error('Unknown flow mode.');
    end

    % -------------------------
    % 4) NUMERICAL SETUP
    % -------------------------
    D = K ./ (Rho .* Cp);
    AdvCoeff = (cfg.rho_f * cfg.cp_f) ./ (Rho .* Cp);

    kappa_rep = median(D(:));
    t_diff_sec = cfg.Ly^2 / (pi^2 * kappa_rep);
    sec_per_yr = 365.25 * 24 * 3600;
    runtime_yr = cfg.runtime_factor * t_diff_sec / sec_per_yr;

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

    T = cfg.T_surf + cfg.grad_init * cfg.Y;
    T = apply_temperature_BC(T, cfg);

    % -------------------------
    % 5) SOLVER LOOP
    % -------------------------
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
                Tnew = apply_temperature_BC(T + 0.5*dt*(k1 + k2), cfg);

            otherwise
                error('Unsupported time scheme.');
        end
        T = Tnew;
    end

    % -------------------------
    % 6) METRICS
    % -------------------------
    % Profile indices
    idx_he1   = round(3000  / cfg.dx) + 1;
    idx_he2   = round(9000  / cfg.dx) + 1;
    idx_fault = round(13000 / cfg.dx) + 1;
    idx_right = round(15000 / cfg.dx) + 1;

    % He1 observational tie-in
    He1_depth_m = [51, 150, 249, 351, 450, 550, 649, 750, 850]';
    He1_temp_C  = [15.1, 15.9, 21.5, 21.0, 25.0, 28.3, 31.6, 34.2, 38.5]';
    T_model_he1 = interp1(cfg.y(:), T(:, idx_he1), He1_depth_m, 'linear', 'extrap');
    rmse_he1 = sqrt(mean((T_model_he1 - He1_temp_C).^2));

    % Isotherm depths
    targets = [50 70 100 120];
    d_he1   = estimate_isotherm_depths(T(:, idx_he1),   cfg.y, targets) / 1000;
    d_he2   = estimate_isotherm_depths(T(:, idx_he2),   cfg.y, targets) / 1000;
    d_fault = estimate_isotherm_depths(T(:, idx_fault), cfg.y, targets) / 1000;
    d_right = estimate_isotherm_depths(T(:, idx_right), cfg.y, targets) / 1000;

    % Outputs
    out = struct();
    out.T = T;
    out.rmse_he1 = rmse_he1;
    out.maxT = max(T(:));

    out.verify = struct();
    out.verify.dt = dt;
    out.verify.dt_diff = dt_diff;
    out.verify.dt_adv = dt_adv;
    out.verify.runtime_kyr = runtime_yr / 1e3;

    out.iso = struct();
    out.iso.targets = targets;
    out.iso.He1_km = d_he1;
    out.iso.He2_km = d_he2;
    out.iso.Fault_km = d_fault;
    out.iso.Right_km = d_right;
end

%% ========================================================================
% LOCAL NUMERICAL HELPER FUNCTIONS
% ========================================================================

function Tbc = apply_temperature_BC(Tin, cfg)
% Boundary conditions on cell-centre temperatures:
%   Top: fixed T
%   Bottom: fixed geothermal gradient
%   Left/Right: insulated

    Tbc = Tin;
    Tbc(1,:)   = cfg.T_surf;
    Tbc(end,:) = Tbc(end-1,:) + cfg.grad_bottom * cfg.dy;
    Tbc(:,1)   = Tbc(:,2);
    Tbc(:,end) = Tbc(:,end-1);
end

function RHS = compute_rhs_temperature(T, K, Rho, Cp, Qr, Vx, Vy, AdvCoeff, cfg)
% dT/dt = (1/rho/cp) div(k grad T) - AdvCoeff*(v.gradT) + Qr/(rho cp)
% Ghost-cell implementation for BC-consistent stencils

    T = apply_temperature_BC(T, cfg);
    [Tpad, Kpad] = build_ghost_fields(T, K, cfg);

    Tc = Tpad(2:end-1, 2:end-1);
    Te = Tpad(2:end-1, 3:end);
    Tw = Tpad(2:end-1, 1:end-2);
    Tn = Tpad(1:end-2, 2:end-1);  % shallower
    Ts = Tpad(3:end,   2:end-1);  % deeper

    Kc = Kpad(2:end-1, 2:end-1);
    Ke = Kpad(2:end-1, 3:end);
    Kw = Kpad(2:end-1, 1:end-2);
    Kn = Kpad(1:end-2, 2:end-1);
    Ks = Kpad(3:end,   2:end-1);

    % Diffusion (conservative variable-coefficient form)
    k_e = 0.5 * (Kc + Ke);
    k_w = 0.5 * (Kc + Kw);
    k_n = 0.5 * (Kc + Kn);
    k_s = 0.5 * (Kc + Ks);

    div_k_gradT = ...
        (k_e .* (Te - Tc) - k_w .* (Tc - Tw)) / cfg.dx^2 + ...
        (k_s .* (Ts - Tc) - k_n .* (Tc - Tn)) / cfg.dy^2;

    DiffTerm = div_k_gradT ./ (Rho .* Cp);

    % Advection (upwind)
    dTdx = zeros(size(T));
    dTdy = zeros(size(T));

    mask_vx_pos = (Vx > 0);
    mask_vx_neg = (Vx < 0);
    dTdx(mask_vx_pos) = (Tc(mask_vx_pos) - Tw(mask_vx_pos)) / cfg.dx;
    dTdx(mask_vx_neg) = (Te(mask_vx_neg) - Tc(mask_vx_neg)) / cfg.dx;

    % y positive downward
    mask_vy_pos = (Vy > 0);  % downward: upwind from shallower (Tn)
    mask_vy_neg = (Vy < 0);  % upward:   upwind from deeper   (Ts)
    dTdy(mask_vy_pos) = (Tc(mask_vy_pos) - Tn(mask_vy_pos)) / cfg.dy;
    dTdy(mask_vy_neg) = (Ts(mask_vy_neg) - Tc(mask_vy_neg)) / cfg.dy;

    AdvTerm = AdvCoeff .* (Vx .* dTdx + Vy .* dTdy);

    % Source
    SourceTerm = Qr ./ (Rho .* Cp);

    RHS = DiffTerm - AdvTerm + SourceTerm;
end

function [Tpad, Kpad] = build_ghost_fields(T, K, cfg)
% Ghost cells for:
%   top Dirichlet, bottom Neumann, left/right Neumann

    [Ny, Nx] = size(T);
    Tpad = zeros(Ny+2, Nx+2);
    Kpad = zeros(Ny+2, Nx+2);

    Tpad(2:end-1, 2:end-1) = T;
    Kpad(2:end-1, 2:end-1) = K;

    % Left/right insulated
    Tpad(2:end-1,1)   = T(:,1);
    Tpad(2:end-1,end) = T(:,end);
    Kpad(2:end-1,1)   = K(:,1);
    Kpad(2:end-1,end) = K(:,end);

    % Top Dirichlet
    Tpad(1,2:end-1) = 2*cfg.T_surf - T(1,:);
    Kpad(1,2:end-1) = K(1,:);

    % Bottom Neumann
    Tpad(end,2:end-1) = T(end,:) + cfg.grad_bottom * cfg.dy;
    Kpad(end,2:end-1) = K(end,:);

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
% Linear interpolation of depth at target temperature
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