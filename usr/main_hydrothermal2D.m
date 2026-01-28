% main_hydrothermal2D.m
% ***** PRACTICAL 30: 2D DARCY FLOW + HEAT TRANSPORT (PSEUDO-TRANSIENT) *****
% BCs:
%   x-sides: periodic
%   top/bottom: closed (w = 0) for Darcy flow
%   top/bottom: isothermal for T  -> T(z=0)=Ttop, T(z=D)=Tbot (imposed at faces via ghost)

function out = main_hydrothermal2D( ...
    W,D,Nx,Nz,h, ...
    kT0,rho0,cP0,Qr0, ...
    g,alpha,Tref,Ttop,Tbot, ...
    CFL_T,CFL_P,tol_div,maxPIt,maxStep,plotEvery,yr, ...
    KD0)

% -------------------- sanity --------------------
if Nx < 5
    error('Nx must be >= 5 (needed for periodic stencils).');
end
if Nz < 4
    error('Nz must be >= 4 (needed for z ghost layers).');
end

% -------------------- grid ----------------------
xc = h/2:h:W-h/2;
zc = h/2:h:D-h/2;
[Xc,Zc] = meshgrid(xc,zc);

% periodic x index maps
ix3 = [Nx, 1:Nx, 1];              % Nx+2
ix5 = [Nx-1, Nx, 1:Nx, 1, 2];     % Nx+4 (for UPW3 option)

% -------------------- parameters ----------------
kT  = kT0  * ones(Nz,Nx);   % W/m/K
rho = rho0 * ones(Nz,Nx);   % kg/m3 (reference used for cp)
cP  = cP0  * ones(Nz,Nx);   % J/kg/K
Qr  = Qr0  * ones(Nz,Nx);   % W/m3

% -------------------- initial condition ----------
% linear gradient from Ttop to Tbot + small perturbation
T = Ttop + (Tbot - Ttop) * (Zc / D);
T = T + 1.0 * sin(2*pi*Xc/W) .* sin(pi*Zc/D);

% pressure initial (gauge)
P = zeros(Nz,Nx);

t = 0;

figure(1); clf;
set(gcf,'Name',sprintf('Hydrothermal Convection | KD0=%.2e',KD0));


% Choose advection scheme (stable default)
ADVN = 'UPW1';  % 'UPW1' (recommended) or 'UPW3' (may overshoot without limiter)

w=zeros(Nz+1,Nx); u=zeros(Nz,Nx+1); P=zeros(Nz,Nx);

for step = 1:maxStep

    % -------- 3) adaptive dt for temperature ----
    umax = max(abs(u(:)));
    wmax = max(abs(w(:)));

    dt_adv = (h/2) / (umax + wmax + eps);

    % diffusion stability using kappa = k/(rho*cp)
    kappa_eff = max( kT(:) ./ (rho(:).*cP(:)) );
    dt_dff = (h/2)^2 / (kappa_eff + eps);

    dt = CFL_T * min(dt_adv, dt_dff);

    % -------- 4) update temperature (RK2) -------
    k1 = T_rate(T, u, w, kT, rho, cP, Qr, h, ix3, ix5, Ttop, Tbot, ADVN);
    k2 = T_rate(T + 0.5*dt*k1, u, w, kT, rho, cP, Qr, h, ix3, ix5, Ttop, Tbot, ADVN);
    T  = T + dt*k2;

    % -------- 1) T-dependent density ------------
    rhoT = rho0 .* (1 - alpha .* (T - Tref));

    % -------- 2) pseudo-transient pressure solve: div(vD)->0
    [P, u, w, divR, pit] = pressure_pseudotransient(P, rhoT, KD0, g, h, CFL_P, tol_div, maxPIt);

    t = t + dt;

    % -------- diagnostics -----------------------
    % cell-centre velocities for plotting only
    ucc = 0.5*(u(:,1:end-1) + u(:,2:end));
    wcc = 0.5*(w(1:end-1,:) + w(2:end,:));

    vmag = sqrt(ucc.^2 + wcc.^2);
    vmax = max(vmag(:));

    Nu = nusselt_top(T, Ttop, Tbot, kT0, D, h);

    if mod(step, plotEvery)==0 || step==1
        makefig(xc, zc, T, P, ucc, wcc, divR, t/yr, KD0, Nu, vmax, pit);
        drawnow;
    end

    % simple early stop
    rmsDiv = sqrt(mean(divR(:).^2));
    if (rmsDiv < tol_div) && (max(abs(k2(:)))*dt < 1e-8) && (pit < maxPIt)
        fprintf('Converged at step %d, t=%.3g yr\n', step, t/yr);
        break;
    end
end

out.T = T;
out.P = P;
out.u = u;
out.w = w;
out.t = t;
out.Nu = Nu;
out.vmax = vmax;

end

% ===================== Temperature rate ==========================
function dTdt = T_rate(T, u, w, kT, rho, cP, Qr, h, ix3, ix5, Ttop, Tbot, ADVN)
adv = advection_T(T, u, w, h, ix5, Ttop, Tbot, ADVN);   % K/s
dif = diffusion_T(T, kT, h, ix3, Ttop, Tbot);           % W/m3
dTdt = adv + dif./(rho.*cP) + Qr./(rho.*cP);
end

% ===================== Diffusion (periodic x, Dirichlet z at faces) ==========================
function dfdt = diffusion_T(f, k, h, ix, Ttop, Tbot)

% z: Dirichlet via ghost (boundary at faces z=0 and z=D)
fext = [2*Ttop - f(1,:); f; 2*Tbot - f(end,:)];
kext = [k(1,:); k; k(end,:)];

kfz = 0.5*(kext(1:end-1,:) + kext(2:end,:));
qz  = -kfz .* diff(fext,1,1) / h;    % W/m2
dfz = -diff(qz,1,1) / h;             % W/m3

% x: periodic
kfx = 0.5*(k(:,ix(1:end-1)) + k(:,ix(2:end)));
qx  = -kfx .* diff(f(:,ix),1,2) / h;
dfx = -diff(qx,1,2) / h;

dfdt = dfz + dfx;
end

% ===================== Advection (conservative flux form) ==========================
function adv = advection_T(f, u, w, h, ix5, Ttop, Tbot, ADVN)

Nz = size(f,1);
Nx = size(f,2);

% ---- x direction (periodic) ----
% stencil in x (works for UPW1 and UPW3)
f_imm = f(:, ix5(1:end-4));
f_im  = f(:, ix5(2:end-3));
f_ic  = f(:, ix5(3:end-2));
f_ip  = f(:, ix5(4:end-1));
f_ipp = f(:, ix5(5:end));

switch ADVN
    case 'UPW1'
        % values at faces: right (i+1/2) and left (i-1/2)
        f_ip_pos = f_ic;  f_ip_neg = f_ip;
        f_im_pos = f_im;  f_im_neg = f_ic;
    case 'UPW3'
        f_ip_pos = (2*f_ip + 5*f_ic - f_im )./6;
        f_ip_neg = (2*f_ic + 5*f_ip - f_ipp)./6;
        f_im_pos = (2*f_ic + 5*f_im - f_imm)./6;
        f_im_neg = (2*f_im + 5*f_ic - f_ip )./6;
    otherwise
        error('Unknown ADVN scheme.');
end

% face velocities: left/right faces per cell
uL = u(:,1:end-1);      % Nz x Nx  (i-1/2)
uR = u(:,2:end);        % Nz x Nx  (i+1/2)
uL_pos = max(0,uL); uL_neg = min(0,uL);
uR_pos = max(0,uR); uR_neg = min(0,uR);

qR = uR_pos.*f_ip_pos + uR_neg.*f_ip_neg;   % (uT) at i+1/2
qL = uL_pos.*f_im_pos + uL_neg.*f_im_neg;   % (uT) at i-1/2
adv_x = -(qR - qL) / h;

% ---- z direction (Dirichlet at faces, ghost padding) ----
% two ghost layers so UPW3 is possible; UPW1 also works
fpad = [ ...
    2*Ttop - f(2,:); ...
    2*Ttop - f(1,:); ...
    f; ...
    2*Tbot - f(end,:); ...
    2*Tbot - f(end-1,:) ...
];

f_jmm = fpad(1:Nz,     :);
f_jm  = fpad(2:Nz+1,   :);
f_jc  = fpad(3:Nz+2,   :);
f_jp  = fpad(4:Nz+3,   :);
f_jpp = fpad(5:Nz+4,   :);

switch ADVN
    case 'UPW1'
        f_jp_pos = f_jc;  f_jp_neg = f_jp;
        f_jm_pos = f_jm;  f_jm_neg = f_jc;
    case 'UPW3'
        f_jp_pos = (2*f_jp + 5*f_jc - f_jm )./6;
        f_jp_neg = (2*f_jc + 5*f_jp - f_jpp)./6;
        f_jm_pos = (2*f_jc + 5*f_jm - f_jmm)./6;
        f_jm_neg = (2*f_jm + 5*f_jc - f_jp )./6;
end

% face velocities: top/bottom faces per cell
wT = w(1:end-1,:);      % Nz x Nx (j-1/2)
wB = w(2:end  ,:);      % Nz x Nx (j+1/2)
wT_pos = max(0,wT); wT_neg = min(0,wT);
wB_pos = max(0,wB); wB_neg = min(0,wB);

qB = wB_pos.*f_jp_pos + wB_neg.*f_jp_neg;   % (wT) at j+1/2
qT = wT_pos.*f_jm_pos + wT_neg.*f_jm_neg;   % (wT) at j-1/2
adv_z = -(qB - qT) / h;

adv = adv_x + adv_z;
end

% ===================== Pressure pseudo-transient solve ==========================
function [P, u, w, divR, it] = pressure_pseudotransient(P, rhoT, KD0, g, h, CFL_P, tol_div, maxIt)

beta = 0.9;
rmsDiv = 1;
it = 1;
dP = 0;
while rmsDiv > tol_div
    [u, w] = darcy_flux(P, rhoT, KD0, g, h);
    divR   = divergence(u, w, h);

    rmsDiv = sqrt(mean(divR(:).^2));

    % NOTE: KD0 treated as Darcy mobility (as per practical wording).
    % Stability: dtP ~ h^2 / KD0
    dtP = CFL_P * h^2 / (KD0 + eps);
    dP  = - dtP * divR + beta*dP;
    P = P + dP;

    if ~ mod(it,100)
    fprintf(' - it = %d;  res = %1.4e \n',it,rmsDiv)
    end
    it =it+1;
end

[u, w] = darcy_flux(P, rhoT, KD0, g, h);
divR   = divergence(u, w, h);

end

% ===================== Darcy flux: vD = -KD0*(grad P - rho g ez) ==========================
function [u, w] = darcy_flux(P, rhoT, KD0, g, h)

[Nz,Nx] = size(P);
mob = KD0;

% u at x-faces (periodic): Nz x (Nx+1)
Pe   = P(:, [Nx, 1:Nx, 1]);       % Nz x (Nx+2)
dPdx = diff(Pe, 1, 2) / h;        % Nz x (Nx+1)
u    = -mob * dPdx;

% w at z-faces: (Nz+1) x Nx, enforce w=0 at top/bottom by ghost P
Pext = [P(1,:) - h*rhoT(1,:)*g; P; P(end,:) + h*rhoT(end,:)*g]; % (Nz+2) x Nx
rhoext = [rhoT(1,:); rhoT; rhoT(end,:)];

dPdz = diff(Pext, 1, 1) / h;      % (Nz+1) x Nx
rho_face = 0.5*(rhoext(1:end-1,:) + rhoext(2:end,:));

w = -mob * (dPdz - rho_face*g);
end

% ===================== divergence at cell centres ==========================
function div = divergence(u, w, h)
div = (u(:,2:end) - u(:,1:end-1))/h + (w(2:end,:) - w(1:end-1,:))/h;
end

% ===================== Nusselt number at top ==========================
function Nu = nusselt_top(T, Ttop, Tbot, kT0, D, h)
q_top = -kT0 * (T(1,:) - Ttop) / (h/2);
q_ref = -kT0 * (Tbot - Ttop) / D;
Nu = mean(q_top) / (q_ref + eps);
end

% ===================== Plotting ==========================
function makefig(xc, zc, T, P, ucc, wcc, divR, t_yr, KD0, Nu, vmax, pit)

subplot(2,2,1);
imagesc(xc, zc, T); axis equal tight; colorbar;
title(sprintf('T [C], t=%.3g yr', t_yr));
xlabel('x [m]'); ylabel('z [m]');

subplot(2,2,2);
imagesc(xc, zc, P); axis equal tight; colorbar;
title(sprintf('P [Pa] (gauge) | KD0=%g', KD0));
xlabel('x [m]'); ylabel('z [m]');

subplot(2,2,3);
imagesc(xc, zc, sqrt(ucc.^2 + wcc.^2)); axis equal tight; colorbar;
title(sprintf('|v_D| [m/s], vmax=%.2e', vmax));
xlabel('x [m]'); ylabel('z [m]');

subplot(2,2,4);
imagesc(xc, zc, divR); axis equal tight; colorbar;
rmsDiv = sqrt(mean(divR(:).^2));
title(sprintf('div(v_D) rms=%.2e | P-it=%d | Nu=%.3f', rmsDiv, pit, Nu));
xlabel('x [m]'); ylabel('z [m]');

end
