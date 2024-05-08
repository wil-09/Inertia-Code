function GainVg_2D
tic; clear all; close all; clc; format long e;
set(0,'defaultaxesfontsize',20,'defaultaxesfontweight','bold','defaultaxeslinewidth',1);

i = sqrt(-1); alpha0 = 0.25; beta = 100; beta0 = 3.5; gamma = 100; 
chi0 = 80; r = 0.1; D2 = 1000; tau = 1; u = 1; %k = 0.05;
% Parameters characterizing the medium porosity
alpha1 = 0.05; alpha2 = 0.02; m0 = 0; m = 10;

%%%%%%%%%%%%%%%%%%%%% grid space %%%%%%%%%%%%%%%%%%%     
kd = -20; kf = 20; dk = 0.001;
z = []; x = []; zvg = []; %matrix to keep calculated variables
k = kd;
while k < kf
    gaini = []; Vgi = [];
    n0 = beta0; c0 = beta0.*n0./(1 + gamma.*n0);
    X = 1 + r + 1./tau + k.^2.*(1 - alpha2.*c0.^m0) + 1i.*k.*u;
    Y = (r + 1./tau).*(1 + k.^2 - alpha2.*k.^2.*c0.^m0 + 1i.*k.*u) + (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 1i.*k.*u)./tau;
    Z = (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 1i.*k.*u).*(1 + k.^2.*(1 - alpha2.*c0.^m0) + 1i.*k.*u)./tau - chi0.*n0.*beta0.*k.^2./(tau.*(1 + beta.*c0).^2.*(1 + gamma.*n0).^2);
    omega1 = (1./6).*( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) - (6.*((1./3).*Y - (1./9).*X.^2))./( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) - (1./3).*X;
    omega2 = - (1./12).*( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) + (3.*((1./3).*Y - (1./9).*X.^2))./( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) - (1./3).*X + (1./2.*i).*sqrt(3).*((1./6).*( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) + (6.*((1./3).*Y - (1./9).*X.^2))./( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3));
    omega3 = - (1./12).*( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) + (3.*((1./3).*Y - (1./9).*X.^2))./( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) - (1./3).*X - (1./2.*i).*sqrt(3).*((1./6).*( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3) + (6.*((1./3).*Y - (1./9).*X.^2))./( - 8.*X.^3 + 36.*Y.*X - 108.*Z + 12.*sqrt(12.*X.^3.*Z - 3.*X.^2.*Y.^2 - 54.*X.*Y.*Z + 12.*Y.^3 + 81.*Z.^2)).^(1./3));
    Vg1 = ( - (2.*k.*(1 - alpha2.*c0.^m0) + i.*u).*omega1.^2 - (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + i.*u.*k).*(2.*k - 2.*alpha2.*k.*c0.^m0 + i.*u)./tau - ((r + 1./tau).*(2.*k - 2.*alpha2.*k.*c0.^m0 + i.*u) + (2.*k.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 2.*k.^3.*D2 + i.*u)./tau).*omega1 - (2.*k.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 2.*k.^3.*D2 + i.*u).*(1 + k.^2 - alpha2.*k.^2.*c0.^m0 + i.*u.*k)./tau + 2.*chi0.*n0.*beta0.*k./(tau.*(gamma.*n0 + 1).^2.*(beta.*c0 + 1).^2))./(3.*omega1.^2 + (r + 1./tau).*(1 + k.^2 - alpha2.*k.^2.*c0.^m0 + i.*u.*k) + (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + i.*u.*k)./tau + (2.*(1 + r + 1./tau + k.^2.*(1 - alpha2.*c0.^m0) + i.*u.*k)).*omega1);
    Vg2 = ( - (2.*k.*(1 - alpha2.*c0.^m0) + i.*u).*omega2.^2 - (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + i.*u.*k).*(2.*k - 2.*alpha2.*k.*c0.^m0 + i.*u)./tau - ((r + 1./tau).*(2.*k - 2.*alpha2.*k.*c0.^m0 + i.*u) + (2.*k.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 2.*k.^3.*D2 + i.*u)./tau).*omega2 - (2.*k.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 2.*k.^3.*D2 + i.*u).*(1 + k.^2 - alpha2.*k.^2.*c0.^m0 + i.*u.*k)./tau + 2.*chi0.*n0.*beta0.*k./(tau.*(gamma.*n0 + 1).^2.*(beta.*c0 + 1).^2))./(3.*omega2.^2 + (r + 1./tau).*(1 + k.^2 - alpha2.*k.^2.*c0.^m0 + i.*u.*k) + (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + i.*u.*k)./tau + (2.*(1 + r + 1./tau + k.^2.*(1 - alpha2.*c0.^m0) + i.*u.*k)).*omega2);
    Vg3 = ( - (2.*k.*(1 - alpha2.*c0.^m0) + i.*u).*omega3.^2 - (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + i.*u.*k).*(2.*k - 2.*alpha2.*k.*c0.^m0 + i.*u)./tau - ((r + 1./tau).*(2.*k - 2.*alpha2.*k.*c0.^m0 + i.*u) + (2.*k.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 2.*k.^3.*D2 + i.*u)./tau).*omega3 - (2.*k.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + 2.*k.^3.*D2 + i.*u).*(1 + k.^2 - alpha2.*k.^2.*c0.^m0 + i.*u.*k)./tau + 2.*chi0.*n0.*beta0.*k./(tau.*(gamma.*n0 + 1).^2.*(beta.*c0 + 1).^2))./(3.*omega3.^2 + (r + 1./tau).*(1 + k.^2 - alpha2.*k.^2.*c0.^m0 + i.*u.*k) + (r + k.^2.*(alpha0 - alpha1.*n0.^m + D2.*k.^2) + i.*u.*k)./tau + (2.*(1 + r + 1./tau + k.^2.*(1 - alpha2.*c0.^m0) + i.*u.*k)).*omega3);
    
    omega11 = real(omega1); omega21 = real(omega2); omega31 = real(omega3); 
%     omega11 = imag(omega1); omega21 = imag(omega2); omega31 = imag(omega3); 
    Vg11 = imag(Vg1); Vg21 = imag(Vg2); Vg31 = imag(Vg3);
    v = [omega11 omega21 omega31]; V = [Vg11 Vg21 Vg31];
    pp = zeros(); PPVg = zeros();
    for j = 1:3
        if v(j) > 0
            pp(j) = (v(j));
            PPVg(j) = V(j);
        end
    end
    gain = (pp); gaini = ([gaini gain]);
    vg = (PPVg); Vgi = [Vgi vg];
%     if k == kd
%         y = [y m];
%     end
    gaini = gaini'; z = [z gaini]; 
    Vgi = Vgi'; zvg = [zvg Vgi];
    x = [x k]; k = k + dk;
end

figure 
plot(x, z, 'linewidth', 3)
xlabel('k'); ylabel('G_R') % legend('stand k=4')
%  colormap (jet)
figure
plot(x,zvg, 'linewidth', 3)
xlabel('k'); ylabel('V_g') % legend('stand k=4')

 toc 
 