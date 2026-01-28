% main_solution2D.m
%***** 2D ADVECTION DIFFUSION MODEL  *************************************

%--- 1. Initialise Coordinates (Slide 49) ---
xc = h/2:h:W-h/2;     % x-cell centres
zc = h/2:h:D-h/2;     % z-cell centres
xf = 0:h:W;           % x-cell faces
zf = 0:h:D;           % z-cell faces
[Xc, Zc] = meshgrid(xc, zc); % Create 2D coordinate arrays

%--- 2. Set Time Step ---
kappa  = kT0./rho0./cP0;
dt_adv = (h/2) / max(abs(u0), abs(w0));
dt_dff = (h/2)^2 / kappa;
dt     = CFL * min(dt_adv, dt_dff);

%--- 3. Initialise Ghost Index Maps (Slide 50) ---
switch BC
    case 'periodic'
        % x-direction (Columns)
        ix3 = [Nx, 1:Nx, 1];
        ix5 = [Nx-1, Nx, 1:Nx, 1, 2];
        % z-direction (Rows)
        iz3 = [Nz, 1:Nz, 1];
        iz5 = [Nz-1, Nz, 1:Nz, 1, 2];
    case 'insulating'
        % Zero-gradient (repeat boundaries)
        ix3 = [1, 1:Nx, Nx];
        ix5 = [1, 1, 1:Nx, Nx, Nx];
        iz3 = [1, 1:Nz, Nz];
        iz5 = [1, 1, 1:Nz, Nz, Nz];
end

%--- 4. Initialise Fields (Slide 51) ---
% Coefficients as 2D arrays
kT  = kT0  .* ones(Nz, Nx);
cP  = cP0  .* ones(Nz, Nx);
rho = rho0 .* ones(Nz, Nx);
Qr  = Qr0  .* ones(Nz, Nx);

% Velocity Fields (Slide 51)
w = w0 .* ones(Nz+1, Nx);
u = u0 .* ones(Nz, Nx+1);

% Initial Temperature (use analytical solution in 2D)
T  = analytical(T0, dT, sgm0, kappa, u0, w0, Xc, Zc, D, W, 0);

Tin = T; % Store initial for plotting
Ta  = T; % Analytical solution
t   = 0;
k   = 0;

% Initial Plot (Slide 54)
figure(1); clf;
makefig(xc, zc, T, T-Ta, t/yr);

%--- 5. Main Loop ---
while t <= tend

    t = t + dt;
    k = k + 1;

    % Select Time Integration (Explicit)
    switch TINT
        case 'FE1'
            dTdt = get_rate(T, u, w, kT, rho, cP, Qr, h, ix3, ix5, iz3, iz5, ADVN);

        case 'RK2'
            k1 = get_rate(T,            u, w, kT, rho, cP, Qr, h, ix3, ix5, iz3, iz5, ADVN);
            k2 = get_rate(T + k1*dt/2,  u, w, kT, rho, cP, Qr, h, ix3, ix5, iz3, iz5, ADVN);
            dTdt = k2;
    end

    % Update Temperature
    T = T + dTdt * dt;

    % Update Analytical Solution (2D)
    Ta = analytical(T0, dT, sgm0, kappa, u0, w0, Xc, Zc, D, W, t);

    % Plotting (Slide 54)
    if ~mod(k, nop)
        makefig(xc, zc, T, T-Ta, t/yr);
    end
end

%--- 6. Helper Functions ---

% Get rate of change (Slide 53)
function dTdt = get_rate(T, u, w, kT, rho, cP, Qr, h, ix3, ix5, iz3, iz5, ADVN)
    dTdt = advection(T, u, w, h, ix5, iz5, ADVN) ...
         + diffusion(T, kT, h, ix3, iz3)./(rho.*cP) ...
         + Qr./(rho.*cP);
end

% 2D Diffusion Function (Slide 52)
function dfdt = diffusion(f, k, h, ix, iz)

    % calculate diffusive flux coefficient at cell faces
    kfz = (k(iz(1:end-1), :) + k(iz(2:end), :)) / 2;
    kfx = (k(:, ix(1:end-1)) + k(:, ix(2:end))) / 2;

    % calculate diffusive flux of scalar field f
    qz = -kfz .* diff(f(iz, :), 1, 1) / h;
    qx = -kfx .* diff(f(:, ix), 1, 2) / h;

    % calculate diffusion flux balance for rate of change
    dfdt = - diff(qz, 1, 1) / h ...
           - diff(qx, 1, 2) / h;
end

% 2D Advection Function (UPW1 / UPW3) in both x and z directions
% u is (Nz, Nx+1), w is (Nz+1, Nx) per Slide 51
function adv = advection(f, u, w, h, ix, iz, ADVN)

    % convert staggered velocities to cell-centred for this scheme
    ucc = 0.5 * (u(:, 1:end-1) + u(:, 2:end));      % Nz x Nx
    wcc = 0.5 * (w(1:end-1, :) + w(2:end, :));      % Nz x Nx

    % --- X-Direction Advection (cols) ---
    u_pos = max(0, ucc);
    u_neg = min(0, ucc);

    f_imm = f(:, ix(1:end-4));
    f_im  = f(:, ix(2:end-3));
    f_ic  = f(:, ix(3:end-2));
    f_ip  = f(:, ix(4:end-1));
    f_ipp = f(:, ix(5:end));

    switch ADVN
        case 'UPW1'
            f_ip_pos = f_ic;  f_im_pos = f_im;
            f_ip_neg = f_ip;  f_im_neg = f_ic;
        case 'UPW3'
            f_ip_pos = (2*f_ip + 5*f_ic - f_im )./6;
            f_im_pos = (2*f_ic + 5*f_im - f_imm)./6;
            f_ip_neg = (2*f_ic + 5*f_ip - f_ipp)./6;
            f_im_neg = (2*f_im + 5*f_ic - f_ip )./6;
    end

    qx_pos = u_pos .* (f_ip_pos - f_im_pos) / h;
    qx_neg = u_neg .* (f_ip_neg - f_im_neg) / h;
    adv_x  = -(qx_pos + qx_neg);

    % --- Z-Direction Advection (rows) ---
    w_pos = max(0, wcc);
    w_neg = min(0, wcc);

    f_jmm = f(iz(1:end-4), :);
    f_jm  = f(iz(2:end-3), :);
    f_jc  = f(iz(3:end-2), :);
    f_jp  = f(iz(4:end-1), :);
    f_jpp = f(iz(5:end), :);

    switch ADVN
        case 'UPW1'
            f_jp_pos = f_jc;  f_jm_pos = f_jm;
            f_jp_neg = f_jp;  f_jm_neg = f_jc;
        case 'UPW3'
            f_jp_pos = (2*f_jp + 5*f_jc - f_jm )./6;
            f_jm_pos = (2*f_jc + 5*f_jm - f_jmm)./6;
            f_jp_neg = (2*f_jc + 5*f_jp - f_jpp)./6;
            f_jm_neg = (2*f_jm + 5*f_jc - f_jp )./6;
    end

    qz_pos = w_pos .* (f_jp_pos - f_jm_pos) / h;
    qz_neg = w_neg .* (f_jp_neg - f_jm_neg) / h;
    adv_z  = -(qz_pos + qz_neg);

    adv = adv_x + adv_z;
end

% Analytical solution in 2D (moving + diffusing Gaussian, periodic wrap)
function Ta = analytical(T0, dT, sgm0, kappa, u0, w0, Xc, Zc, D, W, t)

    x_shift = Xc - (W/2 + u0*t);
    x_shift = x_shift - W * round(x_shift/W);

    z_shift = Zc - (D/2 + w0*t);
    z_shift = z_shift - D * round(z_shift/D);

    sgmt = sqrt(sgm0^2 + 2*kappa*t);
    Ta   = T0 + dT .* (sgm0./sgmt).^2 .* exp( -(x_shift.^2 + z_shift.^2) ./ (2*sgmt.^2) );
end

% Plotting Function (Slide 54)
function makefig(x, z, T, Err, t)
    subplot(2,1,1);
    imagesc(x, z, T); axis equal tight; colorbar
    ylabel('z [m]','FontSize',15)
    title(['Temperature [C]; time = ',num2str(t),' [yr]'],'FontSize',17)

    subplot(2,1,2)
    imagesc(x, z, Err); axis equal tight; colorbar
    xlabel('x [m]','FontSize',15)
    ylabel('z [m]','FontSize',15)
    title('Num. Error [C]','FontSize',17)

    drawnow;
end
