function C_rho = C_rho_contour_operator(a,bar_nu,eps,N_FFT)
% C_rho(k+N+1) = constant C of (34) for the column e^{ik(t-1-eps_2 bar u(t))} of L_3 (63), |k| <= N
% (maximum over the lines Im(z) = -+ bar rho). a = Fourier coefficients of bar u, bar_nu = bar nu.

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

j = -N:N;
for k = 1:N_FFT
    beta1 = exp(j*(bar_rho+1i*delta-1i-1i*eps(2)*U1(k))).*exp(1i*j*(k-1)*2*ipi/N_FFT);
    beta2 = exp(j*(-bar_rho+1i*delta-1i-1i*eps(2)*U2(k))).*exp(1i*j*(k-1)*2*ipi/N_FFT);
    G_max1 = G_max1 + abs(beta1); G_max2 = G_max2 + abs(beta2);
end

C_rho = max(G_max1/N_FFT,G_max2/N_FFT);

end
