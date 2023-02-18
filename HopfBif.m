function HopfBif
% This code plots the critical parameters determined by investigating Hopf
% Bifurcations in the model with Inertia in Porous media

clear all; close all; clc
set(0, 'defaultaxesfontsize', 20, 'defaultaxesfontWeight', 'bold', 'defaultaxesLineWidth', 1)
format long g; 

Klim = 50; K0 = -Klim; Kf = Klim; dK = 0.001; K = 0:dK:Kf;
% Parameters of the model
alpha0 = 0.0025; beta = 100; beta0 = 3.5; gamma = 100; 
chi0 = 80; D2 = 1; r = 100; tau = 10; 

% Parameters characterizing the mediu0m porosity
alpha1 = 0.05; alpha2 = 0.02; m = 6; m0 = 5;
% alpha1 = alpha1.^m; alpha2 = alpha2.^m0;

% The steady state of the system reads
n0 = beta0; c0 = beta0.*n0./(1 + gamma.*n0);

taucrit = (-((1 + K.^2 + r - alpha2.*c0.^m0.*K.^2) + ((r.*(1 + K.^2 - alpha2.*K.^2.*c0.^m0)) - ((r + alpha0.*K.^2 + D2.*K.^4 - alpha1.*n0.^m.*K.^2).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0) - (chi0.*n0.*beta0.*K.^2)./((1 + gamma.*n0.^2).^2.*(1 + beta.*c0).^2)))./(1 + r + K.^2.*(1 + alpha0 + D2.*K.^2 - alpha1.*n0.^m - alpha2.*c0.^m0))) + sqrt(((1 + K.^2 + r - alpha2.*c0.^m0.*K.^2) + ((r.*(1 + K.^2 - alpha2.*K.^2.*c0.^m0)) - ((r + alpha0.*K.^2 + D2.*K.^4 - alpha1.*n0.^m.*K.^2).*(1 + K.^2 - alpha2.*K.^2.*c0.^m0) - (chi0.*n0.*beta0.*K.^2)./((1 + gamma.*n0.^2).^2.*(1 + beta.*c0).^2)))./(1 + r + K.^2.*(1 + alpha0 + D2.*K.^2 - alpha1.*n0.^m - alpha2.*c0.^m0))).^2 - 4.*(r.*(1 + K.^2 - alpha2.*K.^2.*c0.^m0)).*(1 + K.^2 + r - alpha2.*c0.^m0.*K.^2)./(1 + r + K.^2.*(1 + alpha0 + D2.*K.^2 - alpha1.*n0.^m - alpha2.*c0.^m0))))./2;
% taucrit2 = (-((1 + r + K.^2.*(1 - alpha2.*c0.^m0)) + ((r.*(1 + K.^2.*(1 - alpha2.*c0.^m0))) - ((r + K.^2.*(alpha0 + D2.*K.^2 - alpha1.*n0.^m)).*(1 + K.^2.*(1 - alpha2.*c0.^m0)) - chi0.*n0.*beta0.*K.^2./((1 + gamma.*n0).^2.*(1 + beta.*c0).^2)))./(1 + r + K.^2.*(1 + alpha0 + D2.*K.^2 - alpha1.*n0.^m - alpha2.*c0.^m0))) + sqrt(((1 + r + K.^2.*(1 - alpha2.*c0.^m0)) + ((r.*(1 + K.^2.*(1 - alpha2.*c0.^m0))) - ((r + K.^2.*(alpha0 + D2.*K.^2 - alpha1.*n0.^m)).*(1 + K.^2.*(1 - alpha2.*c0.^m0)) - chi0.*n0.*beta0.*K.^2./((1 + gamma.*n0).^2.*(1 + beta.*c0).^2)))./(1 + r + K.^2.*(1 + alpha0 + D2.*K.^2 - alpha1.*n0.^m - alpha2.*c0.^m0))).^2 - 4.*(1 + r + K.^2.*(1 - alpha2.*c0.^m0)).*(r.*(1 + K.^2.*(1 - alpha2.*c0.^m0)))./(1 + r + K.^2.*(1 + alpha0 + D2.*K.^2 - alpha1.*n0.^m - alpha2.*c0.^m0))))./2;

% k2 = (1./6).*( - 8.*((D2 + (alpha0 - alpha1.*n0.^m).*((1 - alpha2.*c0.^m0)))./(D2.*(1 - alpha2.*c0.^m0))).^3 + 36.*((alpha0 - alpha1.*n0.^m + r.*(1 - alpha2.*c0.^m0))./(D2.*(1 - alpha2.*c0.^m0)) - chi0.*n0.*beta0./(D2.*(1 - alpha2.*c0.^m0).*(1 + gamma.*n0).^2.*(1 + beta.*c0).^2)).*((D2 + (alpha0 - alpha1.*n0.^m).*((1 - alpha2.*c0.^m0)))./(D2.*(1 - alpha2.*c0.^m0))) - 108.*(r./(D2.*(1 - alpha2.*c0.^m0))) + 12.*sqrt(12.*(r./(D2.*(1 - alpha2.*c0.^m0))).*((D2 + (alpha0 - alpha1.*n0.^m).*((1 - alpha2.*c0.^m0)))./(D2.*(1 - alpha2.*c0.^m0))).^3 - 3.*((alpha0 - alpha1.*n0.^m + r.*(1 - alpha2.*c0.^m0))./(D2.*(1 - alpha2.*c0.^m0)) - chi0.*n0.*beta0./(D2.*(1 - alpha2.*c0.^m0).*(1 + gamma.*n0).^2.*(1 + beta.*c0).^2)).^2.*((D2 + (alpha0 - alpha1.*n0.^m).*((1 - alpha2.*c0.^m0)))./(D2.*(1 - alpha2.*c0.^m0))).^2 - 54.*(r./(D2.*(1 - alpha2.*c0.^m0))).*((alpha0 - alpha1.*n0.^m + r.*(1 - alpha2.*c0.^m0))./(D2.*(1 - alpha2.*c0.^m0)) - chi0.*n0.*beta0./(D2.*(1 - alpha2.*c0.^m0).*(1 + gamma.*n0).^2.*(1 + beta.*c0).^2)).*((D2 + (alpha0 - alpha1.*n0.^m).*((1 - alpha2.*c0.^m0)))./(D2.*(1 - alpha2.*c0.^m0))) + 12.*((alpha0 - alpha1.*n0.^m + r.*(1 - alpha2.*c0.^m0))./(D2.*(1 - alpha2.*c0.^m0)) - chi0.*n0.*beta0./(D2.*(1 - alpha2.*c0.^m0).*(1 + gamma.*n0).^2.*(1 + beta.*c0).^2)).^3 + 81.*(r./(D2.*(1 - alpha2.*c0.^m0))).^2)).^(1./3) - (6.*((1./3).*((alpha0 - alpha1.*n0.^m + r.*(1 - alpha2.*c0.^m0))./(D2.*(1 - alpha2.*c0.^m0)) - chi0.*n0.*beta0./(D2.*(1 - alpha2.*c0.^m0).*(1 + gamma.*n0).^2.*(1 + beta.*c0).^2)) - (1./9).*((D2 + (alpha0 - alpha1.*n0.^m).*((1 - alpha2.*c0.^m0)))./(D2.*(1 - alpha2.*c0.^m0))).^2))./( - 8.*((D2 + (alpha0 - alpha1.*n0.^m).*((1 - alpha2.*c0.^m0)))./(D2.*(1 - alpha2.*c0.^m0))).^3 + 36.*((alpha0 - alpha1.*n0.^m + r.*(1 - alpha2.*c0.^m0))./(D2.*(1 - alpha2.*c0.^m0)) - chi0.*n0.*beta0./(D2.*(1 - alpha2.*c0.^m0).*(1 + gamma.*n0).^2.*(1 + beta.*c0).^2)).*((D2 + (alpha0 - alpha1.*n0.^m).*((1 - alpha2.*c0.^m0)))./(D2.*(1 - alpha2.*c0.^m0))) - 108.*(r./(D2.*(1 - alpha2.*c0.^m0))) + 12.*sqrt(12.*(r./(D2.*(1 - alpha2.*c0.^m0))).*((D2 + (alpha0 - alpha1.*n0.^m).*((1 - alpha2.*c0.^m0)))./(D2.*(1 - alpha2.*c0.^m0))).^3 - 3.*((alpha0 - alpha1.*n0.^m + r.*(1 - alpha2.*c0.^m0))./(D2.*(1 - alpha2.*c0.^m0)) - chi0.*n0.*beta0./(D2.*(1 - alpha2.*c0.^m0).*(1 + gamma.*n0).^2.*(1 + beta.*c0).^2)).^2.*((D2 + (alpha0 - alpha1.*n0.^m).*((1 - alpha2.*c0.^m0)))./(D2.*(1 - alpha2.*c0.^m0))).^2 - 54.*(r./(D2.*(1 - alpha2.*c0.^m0))).*((alpha0 - alpha1.*n0.^m + r.*(1 - alpha2.*c0.^m0))./(D2.*(1 - alpha2.*c0.^m0)) - chi0.*n0.*beta0./(D2.*(1 - alpha2.*c0.^m0).*(1 + gamma.*n0).^2.*(1 + beta.*c0).^2)).*((D2 + (alpha0 - alpha1.*n0.^m).*((1 - alpha2.*c0.^m0)))./(D2.*(1 - alpha2.*c0.^m0))) + 12.*((alpha0 - alpha1.*n0.^m + r.*(1 - alpha2.*c0.^m0))./(D2.*(1 - alpha2.*c0.^m0)) - chi0.*n0.*beta0./(D2.*(1 - alpha2.*c0.^m0).*(1 + gamma.*n0).^2.*(1 + beta.*c0).^2)).^3 + 81.*(r./(D2.*(1 - alpha2.*c0.^m0))).^2)).^(1./3) - (1./3).*((D2 + (alpha0 - alpha1.*n0.^m).*((1 - alpha2.*c0.^m0)))./(D2.*(1 - alpha2.*c0.^m0)))

% tau = taucrit;
% omegaHopf = (r + 1./tau).*(1 + K.^2) + (r + K.^2.*(alpha0 + D2.*K.^2))./tau - K.^2.*(alpha1.*n0.^m + alpha2.*c0.^m0.*(r + 1./tau));

% the phase velocity read
% vph = sqrt(omegaHopf)./K;

plot(K, 1./taucrit, K, 0.*taucrit, 'linewidth', 3)
xlabel 'k'; ylabel '\tau_{crit}'

% figure
% plot(K, taucrit2, K, 0.*taucrit, 'linewidth', 3)
% xlabel 'k'; ylabel '\tau_{crit}'
% figure
% plot(K, omegaHopf, K, 0.*omegaHopf, 'linewidth', 3)
% xlabel 'k'; ylabel '\Omega_{Hopf}^2'

% figure
% plot(K, vph, 'linewidth', 3)
% xlabel 'k'; ylabel 'V_{ph}'