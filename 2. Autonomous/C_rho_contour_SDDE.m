function [C_rho1, C_rho2, C_rho3, C_rho4] = C_rho_contour_SDDE(x,par,bar_nu,N_FFT)
% Constants C of (34), computed as in (35)-(36) (maximum over the lines Im(z) = -+ bar rho), for
% bar v_1,...,bar v_4 of (79), i.e. bar u(alpha_1), bar u(alpha_2), bar u'(alpha_1), bar u'(alpha_2)
% with alpha_j = t - a_j/bar beta_1 - c_j bar u(t).
% x = (bar beta_1, bar beta_2, bar U), bar_nu = bar nu.

omega = intval(x(1));
a = x(3:end);
gamma = par(1); a1 = par(2); a2 = par(3); c1 = par(4); c2 = par(5); k1 = par(6); k2 = par(7);

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
V5_max = zeros(N_FFT,1); V6_max = zeros(N_FFT,1);
V7_max = zeros(N_FFT,1); V8_max = zeros(N_FFT,1);

for k = 1:N_FFT
    % enclosures V_k of Section 2.5.1 of bar v_1,...,bar v_4 on both lines
    beta1 = a.*exp(n*(bar_rho+1i*delta - 1i*a1/omega - 1i*c1*U1(k))).*exp(1i*n*(k-1)*2*ipi/N_FFT);
    beta2 = a.*exp(n*(-bar_rho+1i*delta - 1i*a1/omega - 1i*c1*U2(k))).*exp(1i*n*(k-1)*2*ipi/N_FFT);
    beta3 = a.*exp(n*(bar_rho+1i*delta - 1i*a2/omega - 1i*c2*U1(k))).*exp(1i*n*(k-1)*2*ipi/N_FFT);
    beta4 = a.*exp(n*(-bar_rho+1i*delta - 1i*a2/omega - 1i*c2*U2(k))).*exp(1i*n*(k-1)*2*ipi/N_FFT);

    beta5 = (1i*n).*a.*exp(n*(bar_rho+1i*delta - 1i*a1/omega - 1i*c1*U1(k))).*exp(1i*n*(k-1)*2*ipi/N_FFT);
    beta6 = (1i*n).*a.*exp(n*(-bar_rho+1i*delta - 1i*a1/omega - 1i*c1*U2(k))).*exp(1i*n*(k-1)*2*ipi/N_FFT);
    beta7 = (1i*n).*a.*exp(n*(bar_rho+1i*delta - 1i*a2/omega - 1i*c2*U1(k))).*exp(1i*n*(k-1)*2*ipi/N_FFT);
    beta8 = (1i*n).*a.*exp(n*(-bar_rho+1i*delta - 1i*a2/omega - 1i*c2*U2(k))).*exp(1i*n*(k-1)*2*ipi/N_FFT);
    V1_on_mesh = sum(beta1); V2_on_mesh = sum(beta2);
    V3_on_mesh = sum(beta3); V4_on_mesh = sum(beta4);

    V5_on_mesh = sum(beta5); V6_on_mesh = sum(beta6);
    V7_on_mesh = sum(beta7); V8_on_mesh = sum(beta8);

    V1_max(k) = sup(abs(V1_on_mesh)); V2_max(k) = sup(abs(V2_on_mesh));
    V3_max(k) = sup(abs(V3_on_mesh)); V4_max(k) = sup(abs(V4_on_mesh));

    V5_max(k) = sup(abs(V5_on_mesh)); V6_max(k) = sup(abs(V6_on_mesh));
    V7_max(k) = sup(abs(V7_on_mesh)); V8_max(k) = sup(abs(V8_on_mesh));
end

C1_rho1 = sup(sum(intval(V1_max)))/N_FFT;
C2_rho1 = sup(sum(intval(V2_max)))/N_FFT;
C1_rho2 = sup(sum(intval(V3_max)))/N_FFT;
C2_rho2 = sup(sum(intval(V4_max)))/N_FFT;
C1_rho3 = sup(sum(intval(V5_max)))/N_FFT;
C2_rho3 = sup(sum(intval(V6_max)))/N_FFT;
C1_rho4 = sup(sum(intval(V7_max)))/N_FFT;
C2_rho4 = sup(sum(intval(V8_max)))/N_FFT;

C_rho1 = max(C1_rho1,C2_rho1);
C_rho2 = max(C1_rho2,C2_rho2);
C_rho3 = max(C1_rho3,C2_rho3);
C_rho4 = max(C1_rho4,C2_rho4);

end
