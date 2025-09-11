function Sept25_In_adaptiveRK4
% Adaptive RK4 for the hyperbolic-parabolic chemotaxis model

tic; clear; close all; clc
set(0, 'defaultaxesfontsize', 20, 'defaultaxesfontWeight', 'bold', 'defaultaxesLineWidth', 1)
format long g;

% Gridspace
xlim = 10; x0 = -xlim; xf = xlim; dx = 0.1; x = x0:dx:xf; N = numel(x);
t0 = 0; tf = 0.1; t = t0; dt = 0.01; dt_min = 1e-5; dt_max = 0.1; tol = 1e-2; % tolerance for adaptive step

% Parameters
alpha0 = 0.25; alpha1 = 0.05; alpha2 = 0.02; D2 = 1e2; tau = 1e3; m0 = 10;
beta0 = 3.5; beta = 100; gamma = 100; chi0 = 80; u = 0.05; r = 0.1; m = 10;

% Initial conditions
n0 = beta0; c0 = beta0*n0/(1 + gamma*n0); eps = 0.01;
k = 10*pi/(xf - x0); gain = 0.1; vg = 0.1; omegai = 0.1;
n = n0 + eps*sech(x - vg*t0).*cos(k*x - omegai*t0).*exp(gain*t0);
c = c0 + eps*sech(x - vg*t0).*cos(k*x - omegai*t0).*exp(gain*t0);
J = zeros(1, N);

% Storage
T = t0; Nf = n; Cf = c; Jf = J;
while t < tf
    if t + dt > tf
        dt = tf - t;
    end
    % RK4 full step
    [n1, c1, J1] = rk4_step(n, c, J, dt, x, dx, alpha0, alpha1, alpha2, D2, tau, beta0, beta, gamma, chi0, u, r, m, m0);

    % RK4 two half-steps
    [nh, ch, Jh] = rk4_step(n, c, J, dt/2, x, dx, alpha0, alpha1, alpha2, D2, tau, beta0, beta, gamma, chi0, u, r, m, m0);
    [n2, c2, J2] = rk4_step(nh, ch, Jh, dt/2, x, dx, alpha0, alpha1, alpha2, D2, tau, beta0, beta, gamma, chi0, u, r, m, m0);

    % Error estimate
    err = max([max(abs(n1 - n2)), max(abs(c1 - c2)), max(abs(J1 - J2))]);

    % Adaptive step size control
    if err < tol
        % Accept step
        t = t + dt;
        n = n2; c = c2; J = J2;
        T(end+1) = t;
        Nf(end+1, :) = n; Cf(end+1, :) = c; Jf(end+1, :) = J;
        % Increase dt for next step
        dt = min(dt * min(2, 0.9*(tol/err)^0.2), dt_max);
    else
        % Reject step, decrease dt
        dt = max(dt * max(0.1, 0.9*(tol/err)^0.25), dt_min);
    end
end

% Plotting
figure; mesh(x, T, Nf); grid on
xlabel('x'); ylabel('t'); zlabel('n(x,t)');
title('Cell density n(x,t)');

figure; mesh(x, T, Cf); grid on
xlabel('x'); ylabel('t'); zlabel('c(x,t)');
title('Chemoattractant c(x,t)');

toc
end

function [n_new, c_new, J_new] = rk4_step(n, c, J, dt, x, dx, alpha0, alpha1, alpha2, D2, tau, beta0, beta, gamma, chi0, u, r, m, m0)
    % One RK4 step for the system
    [dn1, dc1, dJ1] = chemotaxis_rhs(n, c, J, x, dx, alpha0, alpha1, alpha2, D2, tau, beta0, beta, gamma, chi0, u, r, m, m0);
    [dn2, dc2, dJ2] = chemotaxis_rhs(n + 0.5*dt*dn1, c + 0.5*dt*dc1, J + 0.5*dt*dJ1, x, dx, alpha0, alpha1, alpha2, D2, tau, beta0, beta, gamma, chi0, u, r, m, m0);
    [dn3, dc3, dJ3] = chemotaxis_rhs(n + 0.5*dt*dn2, c + 0.5*dt*dc2, J + 0.5*dt*dJ2, x, dx, alpha0, alpha1, alpha2, D2, tau, beta0, beta, gamma, chi0, u, r, m, m0);
    [dn4, dc4, dJ4] = chemotaxis_rhs(n + dt*dn3, c + dt*dc3, J + dt*dJ3, x, dx, alpha0, alpha1, alpha2, D2, tau, beta0, beta, gamma, chi0, u, r, m, m0);

    n_new = n + 1 * (dn1 + 2*dn2 + 2*dn3 + dn4) / 6;
    c_new = c + 1 * (dc1 + 2*dc2 + 2*dc3 + dc4) / 6;
    J_new = J + 1 * (dJ1 + 2*dJ2 + 2*dJ3 + dJ4) / 6;
end

function [dn, dc, dJ] = chemotaxis_rhs(n, c, J, x, dx, alpha0, alpha1, alpha2, D2, tau, beta0, beta, gamma, chi0, u, r, m, m0)
    N = numel(x);
    dn = zeros(1, N); dc = zeros(1, N); dJ = zeros(1, N);

    % Neumann BCs (zero-flux)
    n_ext = [n(1), n(1), n, n(end), n(end)];
    c_ext = [c(1), c(1), c, c(end), c(end)];
    for j = 1:N
        jp2 = j + 2; jp1 = j + 1; jm1 = j - 1 + 2; jm2 = j - 2 + 2;
        dn(j) = J(j);
        % Example: simple diffusion and reaction (customize as needed)
        dJ(j) = -D2 * (n_ext(jm2) - 4*n_ext(jm1) + 6*n_ext(jp1) - 4*n_ext(jp2) + n_ext(jp2+1)) / (tau * dx^4) ...
                + alpha0 / (tau * dx^2) * (n_ext(jm1) - 2*n_ext(jp1) + n_ext(jp2)) ...
                - alpha1 * ((n_ext(jp2))^(1+m) - 2*(n_ext(jp1))^(1+m) + (n_ext(jp1))^(1+m)) / (tau * dx^2 * (1+m)) ...
                - chi0 * ((n_ext(jp1) - n_ext(jm1)) * (c_ext(jp1) - c_ext(jm1))) / (tau * dx^2 * (1 + beta * c(j))^2) ...
                + beta0 * n(j) / (1 + gamma * n(j)) ...
                + r * n(j) / tau;
        dc(j) = c_ext(jp2) / dx^2 - 2 * c_ext(jp1) / dx^2 + c_ext(j) / dx^2 ...
                - alpha2 * ((c_ext(jp2))^(1+m0) - 2*(c_ext(jp1))^(1+m0) + (c_ext(j))^(1+m0)) / (dx^2 * (1+m0)) ...
                + beta0 * n(j) / (1 + gamma * n(j));
    end
end

function y = sech(x)
    y = 1 ./ cosh(x);
end
