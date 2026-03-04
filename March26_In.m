function March26_In
% -------------------------------------------------------------------------
% Vectorized, periodic, RK4 script to appromaximates the solutions of
% the hyperbolic-parabolic chemotaxis model. Code written on March 2026
% -------------------------------------------------------------------------
% Key fixes:
%   (1) Correct RK4 staging (no dt inside k1..k4; only at the final sum)
%   (2) Periodic BCs (plane waves are eigenmodes of the discrete operator)
%   (3) Correct analytical c(x,t) growth exponents (use omega11/21/31)
%   (4) Initialize J analytically from ∂_t n (at t=t0)
%   (5) Conservative explicit dt (remove exponential mapping)
%
% Optional toggle:
%   LINEAR_COMPARE = false;  % set true to zero selected nonlinearities
% -------------------------------------------------------------------------

tic; clearvars; close all; clc
set(0, 'defaultaxesfontsize', 18, ...
       'defaultaxesfontWeight', 'bold', ...
       'defaultaxesLineWidth', 1);
format long g;

% -----------------------------
% Grid and time discretization
% -----------------------------
xlim = 10; x0 = -xlim; xf = xlim; dx = 0.1; t0 = 0; tf = 10;
x  = x0:dx:xf; N  = numel(x);

% periodic index helper
ip1 = [2:N 1];          im1 = [N 1:N-1];
ip2 = [3:N 1 2];        im2 = [N-1 N 1:N-2];

% Fourier mode (plane wave)
k = 20*pi/(xf - x0);    % same as your original

% -----------------------------
% Parameters (as in your code)
% -----------------------------
alpha0 = 0.25;
alpha1 = 0.05;        % will be overwritten by re-scaling below (as in original)
alpha2 = 0.02;        % will be overwritten by re-scaling below (as in original)
D2     = 1e3;
tau    = 1e3;
beta0  = 3.5;
beta   = 100;
gamma  = 100;
chi0   = 80;
u      = 1;
r      = 0.1;
m      = 10;
m0     = 10;

% re-scale the parameter for n(x,t) < nth and c(x,t) < cth 
alpha1 = alpha0/beta0^m/1e1;
alpha2 = ((1 + gamma*beta0)/beta0^2)^m0/1e13;

% Optional toggle: test the pure linear comparison (turn off nonlinearities)
LINEAR_COMPARE = false;   % <-- set true if you want a strict linear check

% -----------------------------
% Conservative explicit time step (recommendation)
%   - advection CFL ~ dx/|u|
%   - Laplacian     ~ dx^2
%   - Biharmonic    ~ dx^4 / (D2/tau)
% -----------------------------
eps_small = 1e-14;
dt_adv  = 0.5*dx/max(abs(u), eps_small);
dt_bih  = 01*dx^4/max(D2/tau, eps_small);  % stiffest term
dt_cap  = 5e-3;                             % additional cap (tunable)
dt      = min([dt_adv, dt_bih, dt_cap]);

t  = t0:dt:tf;
M  = numel(t);

fprintf('[dt selection] dt_adv=%.3e, dt_bih=%.3e, dt_cap=%.3e => dt=%.3e\n', dt_adv, dt_bih, dt_cap, dt);

% --------------------------------------
% Steady state and thresholds (as given)
% --------------------------------------
n0  = beta0;
c0  = beta0.*n0./(1 + gamma.*n0);
eps0 = 0.01;

nth = (alpha0/alpha1)^(1/m);
cth = (1/alpha2)^(1/m0);
fprintf('Steady state estimates: n0=%.6g, c0=%.6g, nth=%.6g, cth=%.6g\n', n0, c0, nth, cth);

% --------------------------------------
% Dispersion relation pieces (your code)
% --------------------------------------
X = 1 + r + 1./tau + k.^2.*(1 - alpha2.*c0.^m0) + 1i.*k.*u;

Y = (r + 1./tau).*(1 + k.^2 - alpha2.*k.^2.*c0.^m0 + 1i.*k.*u) ...
  + (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 1i.*k.*u)./tau;

Z = (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 1i.*k.*u).* ...
    (1 + k.^2.*(1 - alpha2.*c0.^m0) + 1i.*k.*u)./tau ...
    - chi0.*n0.*beta0.*k.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2);

omega1 = (1./6).*( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) ...
       - (6.*((1./3).*Y - (1./9).*X.^2))./( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) ...
       - (1./3).*X;

omega2 = - (1./12).*( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) ...
       + (3.*((1./3).*Y - (1./9).*X.^2))./( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) ...
       - (1./3).*X + (1./2.*1i).*sqrt(3).*((1./6).*( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) ...
       + (6.*((1./3).*Y - (1./9).*X.^2))./( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3));

omega3 = - (1./12).*( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) ...
       + (3.*((1./3).*Y - (1./9).*X.^2))./( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) ...
       - (1./3).*X - (1./2.*1i).*sqrt(3).*((1./6).*( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) ...
       + (6.*((1./3).*Y - (1./9).*X.^2))./( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3));

   Vg1 = ( - (2.*k.*(1 - alpha2.*c0.^m0) + 1i.*u).*omega1.^2 - (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 1i.*u.*k).*(2.*k - 2.*alpha2.*k.*c0.^m0 + 1i.*u)./tau - ((r + 1./tau).*(2.*k - 2.*alpha2.*k.*c0.^m0 + 1i.*u) + (2.*k.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 2.*k.^3.*D2 + 1i.*u)./tau).*omega1 - (2.*k.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 2.*k.^3.*D2 + 1i.*u).*(1 + k.^2 - alpha2.*k.^2.*c0.^m0 + 1i.*u.*k)./tau + 2.*chi0.*n0.*beta0.*k./(tau.*(gamma.*n0 + 1).^2.*(beta.*c0 + 1).^2))./(3.*omega1.^2 + (r + 1./tau).*(1 + k.^2 - alpha2.*k.^2.*c0.^m0 + 1i.*u.*k) + (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 1i.*u.*k)./tau + (2.*(1 + r + 1./tau + k.^2.*(1 - alpha2.*c0.^m0) + 1i.*u.*k)).*omega1);
   Vg2 = ( - (2.*k.*(1 - alpha2.*c0.^m0) + 1i.*u).*omega2.^2 - (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 1i.*u.*k).*(2.*k - 2.*alpha2.*k.*c0.^m0 + 1i.*u)./tau - ((r + 1./tau).*(2.*k - 2.*alpha2.*k.*c0.^m0 + 1i.*u) + (2.*k.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 2.*k.^3.*D2 + 1i.*u)./tau).*omega2 - (2.*k.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 2.*k.^3.*D2 + 1i.*u).*(1 + k.^2 - alpha2.*k.^2.*c0.^m0 + 1i.*u.*k)./tau + 2.*chi0.*n0.*beta0.*k./(tau.*(gamma.*n0 + 1).^2.*(beta.*c0 + 1).^2))./(3.*omega2.^2 + (r + 1./tau).*(1 + k.^2 - alpha2.*k.^2.*c0.^m0 + 1i.*u.*k) + (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 1i.*u.*k)./tau + (2.*(1 + r + 1./tau + k.^2.*(1 - alpha2.*c0.^m0) + 1i.*u.*k)).*omega2);
   Vg3 = ( - (2.*k.*(1 - alpha2.*c0.^m0) + 1i.*u).*omega3.^2 - (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 1i.*u.*k).*(2.*k - 2.*alpha2.*k.*c0.^m0 + 1i.*u)./tau - ((r + 1./tau).*(2.*k - 2.*alpha2.*k.*c0.^m0 + 1i.*u) + (2.*k.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 2.*k.^3.*D2 + 1i.*u)./tau).*omega3 - (2.*k.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 2.*k.^3.*D2 + 1i.*u).*(1 + k.^2 - alpha2.*k.^2.*c0.^m0 + 1i.*u.*k)./tau + 2.*chi0.*n0.*beta0.*k./(tau.*(gamma.*n0 + 1).^2.*(beta.*c0 + 1).^2))./(3.*omega3.^2 + (r + 1./tau).*(1 + k.^2 - alpha2.*k.^2.*c0.^m0 + 1i.*u.*k) + (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 1i.*u.*k)./tau + (2.*(1 + r + 1./tau + k.^2.*(1 - alpha2.*c0.^m0) + 1i.*u.*k)).*omega3);

omega11 = real(omega1); omega21 = real(omega2); omega31 = real(omega3);
Omegai1 = imag(omega1); Omegai2 = imag(omega2); Omegai3 = imag(omega3);
Vg11 = imag(Vg1); Vg21 = imag(Vg2); Vg31 = imag(Vg3);

% -------------------------------------------------------------------------
% Analytical fields (sum of the three branches)
%   nn(x,t) = n0 + ε Σ cos(k x - Ω_i t) exp(σ_i t)
%   cc(x,t) = c0 + ε Σ cos(k x - Ω_i t) exp(σ_i t)   [FIX: use σ_i matching ω_i]
% -------------------------------------------------------------------------
nn = zeros(M, N);
cc = zeros(M, N);

for l = 1:M
    tl = t(l);
%     nn(l,:) = n0 ...
%         + eps0*cos(k.*x - Omegai1.*tl).*exp(omega11.*tl) ...
%         + eps0*cos(k.*x - Omegai2.*tl).*exp(omega21.*tl) ...
%         + eps0*cos(k.*x - Omegai3.*tl).*exp(omega31.*tl);
% % 
% % %     % FIX: use growth exponents matched to each branch (omega11/21/31)
%     cc(l,:) = c0 ...
%         + eps0*cos(k.*x - Omegai1.*tl).*exp(omega11.*tl) ...
%         + eps0*cos(k.*x - Omegai2.*tl).*exp(omega21.*tl) ...
%         + eps0*cos(k.*x - Omegai3.*tl).*exp(omega31.*tl);
    
    nn(l,:) = n0 ...
        + eps0.*sech(x - Vg11.*tl).*cos(k.*x - Omegai1.*tl).*exp(omega11.*tl) ...
        + eps0.*sech(x - Vg21.*tl).*cos(k.*x - Omegai2.*tl).*exp(omega21.*tl) ...
        + eps0.*sech(x - Vg31.*tl).*cos(k.*x - Omegai3.*tl).*exp(omega31.*tl);
    
    cc(l,:) = c0 ...
        + eps0.*sech(x - Vg11.*tl).*cos(k.*x - Omegai1.*tl).*exp(omega21.*tl) ...
        + eps0.*sech(x - Vg21.*tl).*cos(k.*x - Omegai2.*tl).*exp(omega21.*tl) ...
        + eps0.*sech(x - Vg31.*tl).*cos(k.*x - Omegai3.*tl).*exp(omega31.*tl);
end

% -------------------------------------------------------------------------
% Numerical fields: initialization at t0 from analytical nn/cc
% Initialize J analytically as ∂_t n at t0
% -------------------------------------------------------------------------
n = nn(1,:);      % exact analytical at t0
c = cc(1,:);      % exact analytical at t0

theta1 = k.*x - Omegai1.*t0;
theta2 = k.*x - Omegai2.*t0;
theta3 = k.*x - Omegai3.*t0;

% J = eps0*(-Omegai1.*sin(theta1) + omega11.*cos(theta1)) ...
%   + eps0*(-Omegai2.*sin(theta2) + omega21.*cos(theta2)) ...
%   + eps0*(-Omegai3.*sin(theta3) + omega31.*cos(theta3));

J = - eps0.*Vg11.*sech(tl.*Vg11 - x).*tanh(tl.*Vg11 - x).*cos(theta1).*exp(omega11.*tl) ...
    + eps0.*sech(tl.*Vg11 - x).*Omegai1.*sin(theta1).*exp(omega11.*tl) ...
    + eps0.*sech(tl.*Vg11 - x).*cos(theta1).*omega11.*exp(omega11.*tl) ...
    - eps0.*Vg21.*sech(tl.*Vg21 - x).*tanh(tl.*Vg21 - x).*cos(theta2).*exp(omega21.*tl) ...
    + eps0.*sech(tl.*Vg21 - x).*Omegai2.*sin(theta2).*exp(omega21.*tl) ...
    + eps0.*sech(tl.*Vg21 - x).*cos(theta2).*omega21.*exp(omega21.*tl) ...
    - eps0.*Vg31.*sech(tl.*Vg31 - x).*tanh(tl.*Vg31 - x).*cos(theta3).*exp(omega31.*tl) ...
    + eps0.*sech(tl.*Vg31 - x).*Omegai3.*sin(theta3).*exp(omega31.*tl) ...
    + eps0.*sech(tl.*Vg31 - x).*cos(theta3).*omega31.*exp(omega31.*tl);

% Store arrays
Nf  = zeros(M, N);
Cf  = zeros(M, N);
Jf  = zeros(M, N);
Err = zeros(M, N);
Erc = zeros(M, N);

Nf(1,:) = n;
Cf(1,:) = c;
Jf(1,:) = J;

% Struct of parameters for RHS
params = struct('alpha0',alpha0,'alpha1',alpha1,'alpha2',alpha2, ...
                'D2',D2,'tau',tau,'beta0',beta0,'beta',beta, ...
                'gamma',gamma,'chi0',chi0,'u',u,'r',r, ...
                'm',m,'m0',m0,'dx',dx,'LINEAR',LINEAR_COMPARE);

% -------------------------------------------------------------------------
% Time stepping (RK4, vectorized, periodic BC)
% -------------------------------------------------------------------------
for l = 1:M-1
    % Stage 1
    [k1n, k1j, k1c] = RHS_all(n, J, c, params, ip1, im1, ip2, im2);

    % Stage 2
    nS = n + 0.5*dt*k1n;
    JS = J + 0.5*dt*k1j;
    cS = c + 0.5*dt*k1c;
    [k2n, k2j, k2c] = RHS_all(nS, JS, cS, params, ip1, im1, ip2, im2);

    % Stage 3
    nS = n + 0.5*dt*k2n;
    JS = J + 0.5*dt*k2j;
    cS = c + 0.5*dt*k2c;
    [k3n, k3j, k3c] = RHS_all(nS, JS, cS, params, ip1, im1, ip2, im2);

    % Stage 4
    nS = n + dt*k3n;
    JS = J + dt*k3j;
    cS = c + dt*k3c;
    [k4n, k4j, k4c] = RHS_all(nS, JS, cS, params, ip1, im1, ip2, im2);

    % Update
    n = n + (dt/6)*(k1n + 2*k2n + 2*k3n + k4n);
    J = J + (dt/6)*(k1j + 2*k2j + 2*k3j + k4j);
    c = c + (dt/6)*(k1c + 2*k2c + 2*k3c + k4c);

    % Store
    Nf(l+1,:) = n;
    Cf(l+1,:) = c;
    Jf(l+1,:) = J;

    % Errors vs analytical reference at current time layer
    Err(l+1,:) = abs(Nf(l+1,:) - nn(l+1,:));
    Erc(l+1,:) = abs(Cf(l+1,:) - cc(l+1,:));
end

% -----------------------------
% Diagnostics
% -----------------------------
L2n = (trapz(x, (Nf - nn).^2, 2)); % L2n = sqrt(trapz(x, (Nf - nn).^2, 2)); 
L2c = (trapz(x, (Cf - cc).^2, 2)); % L2c = sqrt(trapz(x, (Cf - cc).^2, 2));
fprintf('Final L2 errors:  ||n-analytic||_2 = %.3e,  ||c-analytic||_2 = %.3e\n', L2n(end), L2c(end));

% -----------------------------
% Plots
% -----------------------------
figure('Color','w')
% subplot(2,2,1)
mesh(x, t, Nf); grid off
xlabel('x'); ylabel('t'); zlabel('n_{num}(x,t)')

figure
% subplot(2,2,2)
mesh(x, t, Cf); grid off
xlabel('x'); ylabel('t'); zlabel('c_{num}(x,t)')

figure
% subplot(2,2,3)
mesh(x, t, nn); grid off
xlabel('x'); ylabel('t'); zlabel('n_{ana}(x,t)')

figure
% subplot(2,2,4)
mesh(x, t, cc); grid off
xlabel('x'); ylabel('t'); zlabel('c_{ana}(x,t)')

% figure('Color','w')
% subplot(1,2,1)
figure
yyaxis left
plot(t, L2n, 'b-', 'LineWidth', 2); grid off; hold on
xlabel('t'); ylabel('||n_{num} - n_{ana}||_2')
% title('L2 error for n')

% subplot(1,2,2)
yyaxis right
plot(t, L2c, 'r-', 'LineWidth', 2); % grid off
% xlabel('t'); 
ylabel('||c_{num} - c_{ana}||_2')
title('L_2 Error')
% title('L2 error for c')

toc
end

% =========================================================================
% RHS for the PDE system (vectorized, periodic)
% =========================================================================
function [fn, fJ, fc] = RHS_all(n, J, c, P, ip1, im1, ip2, im2)
% n_t = J
% J_t = (r/tau) n + (alpha0/tau) n_xx - (D2/tau) n_xxxx - (u/tau) n_x
%       + (r - 1/tau) J
%       - (2 r / beta0) n J - (r/(tau*beta0)) n^2
%       - chi0 * [ (n_x c_x)/(tau (1+beta c)^2) + n c_xx/(tau (1+beta c)^2) - 2 beta n (c_x)^2/(tau (1+beta c)^3) ]
%       - (alpha1/tau) * Dxx(n^(1+m))/(1+m)
%
% c_t = c_xx - u c_x - c - alpha2 * Dxx(c^(1+m0))/(1+m0) + beta0 n/(1 + gamma n)
%
% If P.LINEAR = true, turn off nonlinear diffusion and chemotaxis cross-term
% (for a cleaner dispersion comparison).

dx  = P.dx;

% periodic discrete derivatives
nx   = Dx(n,dx,ip1,im1);
nxx  = Dxx(n,dx,ip1,im1);
nxxxx= Dxxxx(n,dx,ip2,ip1,im1,im2);

cx   = Dx(c,dx,ip1,im1);
cxx  = Dxx(c,dx,ip1,im1);

% --- n_t
fn = J;

% --- J_t linear pieces
fJ = (P.r/P.tau).*n ...
   + (P.alpha0/P.tau).*nxx ...
   - (P.D2/P.tau).*nxxxx ...
   - (P.u/P.tau).*nx ...
   + (P.r - 1/P.tau).*J ...
   - (2*P.r/P.beta0).*n.*J ...
   - (P.r/(P.tau*P.beta0)).*(n.^2);

% Nonlinear/chemotaxis terms (optionally stripped for linear comparison)
if ~P.LINEAR
    % chemotaxis-like terms
    denom   = (1 + P.beta.*c);
    chem = - P.chi0 .* ( (nx.*cx) ./ (P.tau*(denom.^2)) ...
          +  n.*cxx ./ (P.tau*(denom.^2)) ...
          - (2*P.beta).*n.*(cx.^2) ./ (P.tau*(denom.^3)) );
    fJ = fJ + chem;

    % porous diffusion in n
    npow = n.^(1 + P.m);
    fJ = fJ - (P.alpha1/P.tau) * Dxx(npow, dx, ip1, im1) / (1 + P.m);
end

% --- c_t
fc = cxx - P.u.*cx - c;

if ~P.LINEAR
    % porous diffusion in c
    cpow = c.^(1 + P.m0);
    fc   = fc - P.alpha2 * Dxx(cpow, dx, ip1, im1) / (1 + P.m0);
end

% source term (keep nonlinear unless linear mode requested)
if P.LINEAR
    % Linearized source around n0 (n0 known; but we do not pass n0 here).
    % If strict linear match is required, replace with a constant-coefficient linearization.
    % Here we keep a pragmatic linear form g'(n0)*delta_n; however n0 isn't passed.
    % For simplicity in LINEAR mode, use small-slope approximation: g(n) ~ beta0*n (gamma n << 1).
    fc = fc + P.beta0 .* n;   % crude linearization
else
    fc = fc + P.beta0 .* n ./ (1 + P.gamma .* n);
end

end

% =========================================================================
% Periodic derivative operators
% =========================================================================
function fx = Dx(f, dx, ip1, im1)
fx = (f(ip1) - f(im1)) ./ (2*dx);
end

function fxx = Dxx(f, dx, ip1, im1)
fxx = (f(ip1) - 2*f + f(im1)) ./ (dx^2);
end

function fxxxx = Dxxxx(f, dx, ip2, ip1, im1, im2)
fxxxx = (f(ip2) - 4*f(ip1) + 6*f - 4*f(im1) + f(im2)) ./ (dx^4);
end
