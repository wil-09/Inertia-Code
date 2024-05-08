function Comb_taucrit
% This code plots the critical parameters determined by investigating Hopf
% Bifurcations in the model with Inertia in Porous media

clear all; close all; clc
set(0, 'defaultaxesfontsize', 20, 'defaultaxesfontWeight', 'bold', 'defaultaxesLineWidth', 1)
format long g; 

Klim = 10; K0 = -Klim; Kf = Klim; dK = 0.001; K = 0:dK:Kf;
% Parameters of the model
alpha0 = 0.25; beta = 100; beta0 = 3.5; gamma = 100; 
chi0 = 80; D2 = 1; r = 0.1;

% Parameters characterizing the mediu0m porosity
alpha1 = 0.05; alpha2 = 0.02; m = 12; m0 = 10;
% alpha1 = alpha1.^m; alpha2 = alpha2.^m0;

% The steady state of the system reads
% for D2 = 10^(0)
% for r = 10^(-3)
for alpha1 = 0.05
n0 = beta0; c0 = beta0.*n0./(1 + gamma.*n0);
taucrit = (-((1 + K.^2 + r - alpha2.*c0.^m0.*K.^2) + ((r.*(1 + K.^2 - alpha2.*K.^2.*c0.^m0)) - ((r + alpha0.*K.^2 + D2.*K.^4 - alpha1.*n0.^m.*K.^2).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0) - (chi0.*n0.*beta0.*K.^2)./((1 + gamma.*n0.^2).^2.*(1 + beta.*c0).^2)))./(1 + r + K.^2.*(1 + alpha0 + D2.*K.^2 - alpha1.*n0.^m - alpha2.*c0.^m0))) + sqrt(((1 + K.^2 + r - alpha2.*c0.^m0.*K.^2) + ((r.*(1 + K.^2 - alpha2.*K.^2.*c0.^m0)) - ((r + alpha0.*K.^2 + D2.*K.^4 - alpha1.*n0.^m.*K.^2).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0) - (chi0.*n0.*beta0.*K.^2)./((1 + gamma.*n0.^2).^2.*(1 + beta.*c0).^2)))./(1 + r + K.^2.*(1 + alpha0 + D2.*K.^2 - alpha1.*n0.^m - alpha2.*c0.^m0))).^2 - 4.*(r.*(1 + K.^2 - alpha2.*K.^2.*c0.^m0)).*(1 + K.^2 + r - alpha2.*c0.^m0.*K.^2)./(1 + r + K.^2.*(1 + alpha0 + D2.*K.^2 - alpha1.*n0.^m - alpha2.*c0.^m0))))./2;
omegai = (1 + r + r./tau + K.^2.*((1 - alpha2.*c0.^m0).*(r + 1./tau) + (alpha0 - alpha1.*n0^m + D2.*K.^2)./tau) )
plot(K, 1./taucrit, 'linewidth', 3); hold on
end

% for D2 = 10^3
% for r = 1
for alpha1 = 0.1
n0 = beta0; c0 = beta0.*n0./(1 + gamma.*n0);
taucrit = (-((1 + K.^2 + r - alpha2.*c0.^m0.*K.^2) + ((r.*(1 + K.^2 - alpha2.*K.^2.*c0.^m0)) - ((r + alpha0.*K.^2 + D2.*K.^4 - alpha1.*n0.^m.*K.^2).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0) - (chi0.*n0.*beta0.*K.^2)./((1 + gamma.*n0.^2).^2.*(1 + beta.*c0).^2)))./(1 + r + K.^2.*(1 + alpha0 + D2.*K.^2 - alpha1.*n0.^m - alpha2.*c0.^m0))) + sqrt(((1 + K.^2 + r - alpha2.*c0.^m0.*K.^2) + ((r.*(1 + K.^2 - alpha2.*K.^2.*c0.^m0)) - ((r + alpha0.*K.^2 + D2.*K.^4 - alpha1.*n0.^m.*K.^2).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0) - (chi0.*n0.*beta0.*K.^2)./((1 + gamma.*n0.^2).^2.*(1 + beta.*c0).^2)))./(1 + r + K.^2.*(1 + alpha0 + D2.*K.^2 - alpha1.*n0.^m - alpha2.*c0.^m0))).^2 - 4.*(r.*(1 + K.^2 - alpha2.*K.^2.*c0.^m0)).*(1 + K.^2 + r - alpha2.*c0.^m0.*K.^2)./(1 + r + K.^2.*(1 + alpha0 + D2.*K.^2 - alpha1.*n0.^m - alpha2.*c0.^m0))))./2;
omegai = (1 + r + r./tau + K.^2.*((1 - alpha2.*c0.^m0).*(r + 1./tau) + (alpha0 - alpha1.*n0^m + D2.*K.^2)./tau) )
plot(K, 1./taucrit, 'linewidth', 3); hold on
end

% for D2 = 10^4
% for r = 10^3
for alpha1 = 1
n0 = beta0; c0 = beta0.*n0./(1 + gamma.*n0);
taucrit = (-((1 + K.^2 + r - alpha2.*c0.^m0.*K.^2) + ((r.*(1 + K.^2 - alpha2.*K.^2.*c0.^m0)) - ((r + alpha0.*K.^2 + D2.*K.^4 - alpha1.*n0.^m.*K.^2).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0) - (chi0.*n0.*beta0.*K.^2)./((1 + gamma.*n0.^2).^2.*(1 + beta.*c0).^2)))./(1 + r + K.^2.*(1 + alpha0 + D2.*K.^2 - alpha1.*n0.^m - alpha2.*c0.^m0))) + sqrt(((1 + K.^2 + r - alpha2.*c0.^m0.*K.^2) + ((r.*(1 + K.^2 - alpha2.*K.^2.*c0.^m0)) - ((r + alpha0.*K.^2 + D2.*K.^4 - alpha1.*n0.^m.*K.^2).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0) - (chi0.*n0.*beta0.*K.^2)./((1 + gamma.*n0.^2).^2.*(1 + beta.*c0).^2)))./(1 + r + K.^2.*(1 + alpha0 + D2.*K.^2 - alpha1.*n0.^m - alpha2.*c0.^m0))).^2 - 4.*(r.*(1 + K.^2 - alpha2.*K.^2.*c0.^m0)).*(1 + K.^2 + r - alpha2.*c0.^m0.*K.^2)./(1 + r + K.^2.*(1 + alpha0 + D2.*K.^2 - alpha1.*n0.^m - alpha2.*c0.^m0))))./2;
omegai = (1 + r + r./tau + K.^2.*((1 - alpha2.*c0.^m0).*(r + 1./tau) + (alpha0 - alpha1.*n0^m + D2.*K.^2)./tau) )
plot(K, 1./taucrit, 'linewidth', 3)
xlabel 'k'; ylabel '\tau_{crit}'
end
% legend('m = 5', 'm = 7', 'm = 10')
legend('\alpha_1 = 0.05', '\alpha_1 = 0.1', '\alpha_1 = 1')
% legend('D_2 = 10^{-3}', 'D_2 = 1', 'D_2 = 10^3')
% legend('r = 10^{-3}', 'r = 1', 'r = 10^3')
