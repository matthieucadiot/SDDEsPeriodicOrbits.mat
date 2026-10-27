function F = F_threshold_rig(x,par,C_g1,C_g2,C_vw,nu,N_FFT)
% Rigorous enclosure of F(bar x) = (xi(bar x), f(bar x)) of (111)-(113) for the threshold SDDE,
% in Fourier coefficients, |n| <= N, including the tails |k| > N in the zeroth rows of L_1 g and in
% the integral of xi_1.
%   x     = (beta,u), beta = (eta,omega,beta_3,beta_4), u = (w,tau)
%   par   = parameters (see DF_threshold.m)
%   C_g1, C_g2, C_vw = constants C of (34) for bar g_1, bar g_2 and v(bar w) (C_rho_contour_SDDE.m)
%   nu = bar nu,  N_FFT = number of grid points (power of 2, N_FFT >= 4N)

beta = x(1:4);
eta = beta(1); omega = intval(beta(2));
u = x(5:end); N = (length(u)-2)/4;
w   = u(1:2*N+1);
tau = u(2*N+2:4*N+2);

a       = par(1);
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

if abs(log2(N_FFT)-ceil(log2(N_FFT)))>0 || N>N_FFT/4
    disp('Wrong size for N_FFT'), return
end

n = intval((-N:N)');

% w(t_k), tau(t_k) (inverse FFT, Section 2.5.2)
w_pad   = [intval(zeros(N_FFT/2-N,1)); w;   intval(zeros(N_FFT/2-N-1,1))];
tau_pad = [intval(zeros(N_FFT/2-N,1)); tau; intval(zeros(N_FFT/2-N-1,1))];
w_k   = ifftshift(verifyfft(fftshift(w_pad),  -1))*intval(N_FFT);
tau_k = ifftshift(verifyfft(fftshift(tau_pad),-1))*intval(N_FFT);

t = (2*intval('pi')/N_FFT)*(-N_FFT/2:N_FFT/2-1)';

% w(t_k - tau(t_k)), cf. (37)
wd_k = intval(zeros(N_FFT,1));
for k = 0:N_FFT-1
    wd_k(k+1) = real(sum( w.*exp(1i*n*(t(k+1)-tau_k(k+1))) ));
end

% bar g_1, bar g_2 (Section 4.1.1) on the grid
vW_k   = v(w_k);     % v(w(t))
vWd_k  = v(wd_k);    % v(w(t-tau))
qWd_k  = q(wd_k);    % q(w(t-tau))

g1_k = omega*b*exp(-mu*omega*tau_k).*vW_k.*qWd_k./vWd_k - omega*gamma*w_k;
g2_k = 1 - vW_k./vWd_k;

% Fourier coefficients |n| <= N with the aliasing bound (38)
G1bar = fftshift(verifyfft(ifftshift(g1_k))/intval(N_FFT));
G2bar = fftshift(verifyfft(ifftshift(g2_k))/intval(N_FFT));
VWbar = fftshift(verifyfft(ifftshift(vW_k))/intval(N_FFT));

% 1/(nu^N_FFT-1) written as nu^(-N_FFT)/(1-nu^(-N_FFT)) to avoid overflow
nu_N = exp(-N_FFT*log(intval(nu))); nu_N = nu_N/(1-nu_N);
alias_g1 = infsup(-1,1)*C_g1*(nu.^n + nu.^(-n))*nu_N;
alias_g2 = infsup(-1,1)*C_g2*(nu.^n + nu.^(-n))*nu_N;
alias_vw = infsup(-1,1)*C_vw*(nu.^n + nu.^(-n))*nu_N;

g1      = G1bar(N_FFT/2-N+1:N_FFT/2+N+1,1) + alias_g1;
g2      = G2bar(N_FFT/2-N+1:N_FFT/2+N+1,1) + alias_g2;
v_circ_w= VWbar(N_FFT/2-N+1:N_FFT/2+N+1,1) + alias_vw;

% xi of (113)
xi = intval(zeros(4,1));

tau_at_0 = real(sum(intval(tau)));
int_xi_1 = (1 - exp(-1i*n*tau_at_0))./(1i*n);  int_xi_1(N+1) = tau_at_0;   % vartheta_n (defined after (125))

% tail of the integral in xi_1: |vartheta_n| <= 2/|n| and |(v o w)_n| <= C_vw/nu^|n|, hence
% sum_{|n|>N} <= 4*C_vw*nu/((N+1)(nu-1)nu^(N+1))
tail_xi1 = 4*C_vw*nu/((N+1)*(nu-1)*nu^(N+1));

xi(1) = omega*real(sum(v_circ_w.*int_xi_1)) - a + omega*infsup(-1,1)*tail_xi1;
xi(2) = sum(intval(w)) - 1;
xi(3) = g1(N+1);
xi(4) = g2(N+1) + eta;

% L_1 (Section 2.4)
L1 = intval(zeros(2*N+1));
for j = [-N:-1 1:N]
    L1(j+N+1,j+N+1) = 1/(1i*intval(j));
end
for j = [-N:-1 1:N]
    L1(N+1,j+N+1) = (-1)^(j+1)/(1i*intval(j));
end

% f of (111)
one = [intval(zeros(N,1)); intval(1); intval(zeros(N,1))];
f1 = w   - beta(3)*one - L1*g1;
f2 = tau - beta(4)*one - L1*g2;

% tail |k| > N of the zeroth row of L_1 g, bounded with (34)
f1(N+1) = f1(N+1) + infsup(-1,1)*2*C_g1*nu/((N+1)*(nu-1)*nu^(N+1));
f2(N+1) = f2(N+1) + infsup(-1,1)*2*C_g2*nu/((N+1)*(nu-1)*nu^(N+1));

F = [xi; f1; f2];

end
