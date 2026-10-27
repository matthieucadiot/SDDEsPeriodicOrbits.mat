function [C_rho1,C_rho2] = C_rho_contour_operator(x,par,bar_nu,N_FFT)
% C_rho1(k+N+1), C_rho2(k+N+1) = constants C of (34) for the columns e^{ik alpha_1}, e^{ik alpha_2}
% of L31, L32 (DF_SDDE_rig.m), |k| <= N, alpha_j = t - a_j/bar beta_1 - c_j bar u(t) (maximum over
% the lines Im(z) = -+ bar rho). x = (bar beta_1, bar beta_2, bar U), bar_nu = bar nu.

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
U1 = verifyfft(ifftshift(alpha1),-1)*N_FFT;
U2 = verifyfft(ifftshift(alpha2),-1)*N_FFT;

G_max1 = zeros(1,2*N+1);G_max2 = zeros(1,2*N+1);
G_max3 = zeros(1,2*N+1);G_max4 = zeros(1,2*N+1);

j = -N:N;
for k = 1:N_FFT
    beta1 = exp(j*(bar_rho + 1i*delta -1i*a1/omega - 1i*c1*U1(k))).*exp(1i*j*(k-1)*2*ipi/N_FFT);
    beta2 = exp(j*(-bar_rho + 1i*delta -1i*a1/omega - 1i*c1*U2(k))).*exp(1i*j*(k-1)*2*ipi/N_FFT);
    beta3 = exp(j*(bar_rho + 1i*delta -1i*a2/omega - 1i*c2*U1(k))).*exp(1i*j*(k-1)*2*ipi/N_FFT);
    beta4 = exp(j*(-bar_rho + 1i*delta -1i*a2/omega - 1i*c2*U2(k))).*exp(1i*j*(k-1)*2*ipi/N_FFT);
    G_max1 = G_max1 + abs(beta1); G_max2 = G_max2 + abs(beta2);
    G_max3 = G_max3 + abs(beta3); G_max4 = G_max4 + abs(beta4);
end

C_rho1 = max(G_max1/N_FFT,G_max2/N_FFT);
C_rho2 = max(G_max3/N_FFT,G_max4/N_FFT);

end
