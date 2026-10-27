function [C_g1, C_g2, C_vw] = C_rho_contour_SDDE(x,par,bar_nu,N_FFT)
% Constants C of (34), computed as in (35)-(36) (maximum over the lines Im(z) = -+ bar rho,
% bar rho = log(bar nu)), for the components of g(bar beta_2, bar u) (Section 4.1.1) and v(bar w):
%   C_g1, C_g2 : bar g_1 = omega b e^{-mu omega tau} v(w)/v(w(t-tau)) q(w(t-tau)) - omega gamma w,
%                bar g_2 = 1 - v(w)/v(w(t-tau)),  with (w,tau) = (bar w,bar tau), omega = bar beta_2
%   C_vw       : v(bar w) (integral in xi_1 of (113))
% x = (bar beta, bar w, bar tau), par as in DF_threshold.m, N_FFT = number of grid points.

beta  = x(1:4);
omega = intval(beta(2));
u = x(5:end); N = (length(u)-2)/4;
w   = u(1:2*N+1);
tau = u(2*N+2:4*N+2);

b       = par(2);
mu      = par(3);
gamma   = par(4);
m       = mid(par(5));
n_q     = mid(par(6));
q_minus = par(7);
q_plus  = par(8);
Theta_q = par(9);
v_minus = par(10);
v_plus  = par(11);
Theta_v = par(12);

q = @(y)((q_minus*Theta_q^n_q + q_plus*y.^n_q)./(Theta_q^n_q + y.^n_q));
v = @(y)((v_minus*Theta_v^m   + v_plus*y.^m  )./(Theta_v^m   + y.^m  ));

bar_rho = sup(log(intval(bar_nu)));   % bar_rho >= log(bar_nu)
n = (-N:N)';

% cells [t_{k-1},t_k], t_k = 2*pi*k/N_FFT (delta = [0,2*pi/N_FFT] instead of [-pi/N_FFT,pi/N_FFT])
ipi = intval('pi'); delta = infsup(0,1)*2*ipi/N_FFT;

% w(z) and tau(z) for z in [t_{k-1},t_k] -+ i bar_rho, as in (35)
aw1 = [zeros(N_FFT/2-N,1); w.*exp(n*( bar_rho+1i*delta)); zeros(N_FFT/2-N-1,1)];
aw2 = [zeros(N_FFT/2-N,1); w.*exp(n*(-bar_rho+1i*delta)); zeros(N_FFT/2-N-1,1)];
at1 = [zeros(N_FFT/2-N,1); tau.*exp(n*( bar_rho+1i*delta)); zeros(N_FFT/2-N-1,1)];
at2 = [zeros(N_FFT/2-N,1); tau.*exp(n*(-bar_rho+1i*delta)); zeros(N_FFT/2-N-1,1)];

Uw1 = verifyfft(ifftshift(aw1),-1)*N_FFT;  Uw2 = verifyfft(ifftshift(aw2),-1)*N_FFT;
Ut1 = verifyfft(ifftshift(at1),-1)*N_FFT;  Ut2 = verifyfft(ifftshift(at2),-1)*N_FFT;

g1_max1 = zeros(N_FFT,1); g1_max2 = zeros(N_FFT,1);
g2_max1 = zeros(N_FFT,1); g2_max2 = zeros(N_FFT,1);
vw_max1 = zeros(N_FFT,1); vw_max2 = zeros(N_FFT,1);

for k = 1:N_FFT
    % line Im(z) = -bar_rho
    wz = Uw1(k);  tz = Ut1(k);
    % w(z - tau(z))
    wd = sum( w.*exp(n*( bar_rho+1i*delta - 1i*tz)).*exp(1i*n*(k-1)*2*ipi/N_FFT) );
    vW  = v(wz);  vWd = v(wd);  qWd = q(wd);
    ex  = exp(-mu*omega*tz);
    g1z = omega*b*ex*vW*qWd/vWd - omega*gamma*wz;
    g2z = 1 - vW/vWd;
    g1_max1(k) = sup(abs(g1z));
    g2_max1(k) = sup(abs(g2z));
    vw_max1(k) = sup(abs(vW));

    % line Im(z) = +bar_rho
    wz = Uw2(k);  tz = Ut2(k);
    wd = sum( w.*exp(n*(-bar_rho+1i*delta - 1i*tz)).*exp(1i*n*(k-1)*2*ipi/N_FFT) );
    vW  = v(wz);  vWd = v(wd);  qWd = q(wd);
    ex  = exp(-mu*omega*tz);
    g1z = omega*b*ex*vW*qWd/vWd - omega*gamma*wz;
    g2z = 1 - vW/vWd;
    g1_max2(k) = sup(abs(g1z));
    g2_max2(k) = sup(abs(g2z));
    vw_max2(k) = sup(abs(vW));
end

C_g1 = max( sup(sum(intval(g1_max1))) , sup(sum(intval(g1_max2))) )/N_FFT;
C_g2 = max( sup(sum(intval(g2_max1))) , sup(sum(intval(g2_max2))) )/N_FFT;
C_vw = max( sup(sum(intval(vw_max1))) , sup(sum(intval(vw_max2))) )/N_FFT;

end
