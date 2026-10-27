function [C_rho, C_rho_der] = C_rho_contour_SDDE(a,bar_nu,eps,N_FFT)
% Constants C of (34), computed as in (35)-(36) (maximum over the lines Im(z) = -+ bar rho), for
%   C_rho     : v = bar u(t-1-eps_2 bar u(t))   (F_SDDE_rig.m)
%   C_rho_der : bar u'(t-1-eps_2 bar u(t))      (L_2 in DF_SDDE_rig.m)
% a = Fourier coefficients of bar u, bar_nu = bar nu.

N = (length(a)-1)/2;
bar_rho = sup(log(intval(bar_nu)));   % bar_rho >= log(bar_nu)
n = (-N:N)';

ipi = intval('pi'); delta = infsup(0,1)*2*ipi/N_FFT;
alpha1 = a.*exp(n*(bar_rho+1i*delta));
alpha2 = a.*exp(n*(-bar_rho+1i*delta));

alpha1 = [zeros(N_FFT/2-N,1);alpha1;zeros(N_FFT/2-N-1,1)];
alpha2 = [zeros(N_FFT/2-N,1);alpha2;zeros(N_FFT/2-N-1,1)];

% U1(k), U2(k) enclose bar u(t -+ i bar_rho) for t in [t_{k-1},t_k], t_k = 2*pi*k/N_FFT, as in (35)
% (the code uses delta = [0,2*pi/N_FFT] instead of [-pi/N_FFT,pi/N_FFT])
U1 = verifyfft(ifftshift(alpha1),-1)*N_FFT; U2 = verifyfft(ifftshift(alpha2),-1)*N_FFT;

V1_max = zeros(N_FFT,1); V2_max = zeros(N_FFT,1);
V3_max = zeros(N_FFT,1); V4_max = zeros(N_FFT,1);

for k = 1:N_FFT
    % enclosures V_k of Section 2.5.1 of the two functions on both lines
    beta1 = a.*exp(n*(bar_rho+1i*delta-1i-1i*eps(2)*U1(k))).*exp(1i*n*(k-1)*2*ipi/N_FFT);
    beta2 = a.*exp(n*(-bar_rho+1i*delta-1i-1i*eps(2)*U2(k))).*exp(1i*n*(k-1)*2*ipi/N_FFT);
    beta3 = (1i*n).*beta1; beta4 = (1i*n).*beta2;
    V1_on_mesh = sum(beta1); V2_on_mesh = sum(beta2);
    V3_on_mesh = sum(beta3); V4_on_mesh = sum(beta4);
    V1_max(k) = sup(abs(V1_on_mesh)); V2_max(k) = sup(abs(V2_on_mesh));
    V3_max(k) = sup(abs(V3_on_mesh)); V4_max(k) = sup(abs(V4_on_mesh));
end

C1_rho = sup(sum(intval(V1_max)))/N_FFT;
C2_rho = sup(sum(intval(V2_max)))/N_FFT;
C3_rho = sup(sum(intval(V3_max)))/N_FFT;
C4_rho = sup(sum(intval(V4_max)))/N_FFT;

C_rho = max(C1_rho,C2_rho);
C_rho_der = max(C3_rho,C4_rho);

end
