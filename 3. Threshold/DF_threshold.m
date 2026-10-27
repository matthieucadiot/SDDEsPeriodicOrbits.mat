function DF = DF_threshold(x,par,N_FFT)
% Jacobian of the map F of (112) in Fourier coefficients, |n| <= N, for the threshold SDDE
% (non-rigorous, used to build A^{<=N}), with the decomposition (119)-(125).
%   x     = (beta,u): beta = (eta,omega,beta_3,beta_4) ((beta_3,beta_4) = u(-pi) at a zero of F),
%           u = (u_1,u_2) = (w,tau) (Fourier coefficients (u_i)_n, |n| <= N; tau is the
%           rescaled delay, Lemma 21)
%   par   = (a,b,mu,gamma,m,n,q^-,q^+,theta_q,v^-,v^+,theta_v) of (100)-(101)
%   N_FFT = number of grid points (power of 2, N_FFT >= 4N)

beta = x(1:4);
omega = beta(2);
u = x(5:end); N = (length(u)-2)/4;
w = u(1:2*N+1); tau = u(2*N+2:4*N+2);

b = par(2);
mu = par(3);
gamma = par(4);
m = par(5);
n = par(6);
q_minus = par(7);
q_plus = par(8);
Theta_q = par(9);
v_minus = par(10);
v_plus = par(11);
Theta_v = par(12);

q = @(y)((q_minus*Theta_q^n+q_plus*y.^n)./(Theta_q^n+y.^n));
v = @(y)((v_minus*Theta_v^m+v_plus*y.^m)./(Theta_v^m+y.^m));
dqdy = @(y)((n.*y.^(n-1).*Theta_q.^n.*(q_plus - q_minus))./(Theta_q.^n + y.^n).^2);
dvdy = @(y)((m.*y.^(m-1).* Theta_v.^m.*(v_plus - v_minus))./(Theta_v.^m + y.^m).^2);

if abs(log2(N_FFT)-ceil(log2(N_FFT)))>0 || N>N_FFT/4
    disp('Wrong size for N_FFT'), return
end

% zero padding to size N_FFT (Section 2.5.2)
w_pad = [zeros(N_FFT/2-N,1);w;zeros(N_FFT/2-N-1,1)];
tau_pad = [zeros(N_FFT/2-N,1);tau;zeros(N_FFT/2-N-1,1)];

% w(t_k), tau(t_k) (inverse FFT)
w_k = ifftshift(ifft(fftshift(w_pad)))*N_FFT;
tau_k = ifftshift(ifft(fftshift(tau_pad)))*N_FFT;

% y = w(t_k - tau(t_k)) and ty = w'(t_k - tau(t_k))
t = (2*pi/N_FFT)*(-N_FFT/2:N_FFT/2-1)';
y = zeros(N_FFT,1);
ty = zeros(N_FFT,1);

n1 = (-N:N)';
wp = w.*(1i.*n1);

for k = 0:N_FFT-1
    y(k+1) = real(sum(w.*exp(1i*n1*(t(k+1)-tau_k(k+1)))));
    ty(k+1) = real(sum(wp.*exp(1i*n1*(t(k+1)-tau_k(k+1)))));
end

v_circ_w = v(w_k);
vp_circ_w = dvdy(w_k);

% derivatives of calG at (w(t-tau), u, omega) on the grid, cf. (120)-(124)
d_1G_1 = omega*b*exp(-mu*omega*tau_k).*v_circ_w.*dqdy(y)./v(y)-omega*b*exp(-mu*omega*tau_k).*dvdy(y).*v_circ_w.*q(y)./(v(y).^2);
d_1G_2 = dvdy(y).*v_circ_w./(v(y).^2);

d_1G_1_ty = d_1G_1.*ty;
d_1G_2_ty = d_1G_2.*ty;

d_2G_11 = omega*b*exp(-mu*omega*tau_k).*vp_circ_w.*q(y)./v(y)-omega*gamma*ones(N_FFT,1);
d_2G_12 = -mu*omega^2*b*exp(-mu*omega*tau_k).*v_circ_w.*q(y)./v(y);
d_2G_21 = -vp_circ_w./v(y);

d_3G_1 = b*exp(-mu*omega*tau_k).*v_circ_w.*q(y).*(ones(N_FFT,1)-omega*mu*tau_k)./v(y)-gamma*w_k;

% columns of d_u g(bar x) (122), from the functions d_u g(bar x) e_k on the grid, and zeroth rows
% (L_1 U)_0 = sum(c0.*U) and int_{-tau(0)}^0 h = sum(th0.*U), with th0 = vartheta (defined after
% (125)), which involve all the Fourier coefficients
nF = (-N_FFT/2:N_FFT/2-1)'; tau_0 = real(sum(tau));
c0 = (-1).^(nF+1)./(1i*(nF+(nF==0))); c0([1 N_FFT/2+1]) = 0;
th0 = (1-exp(-1i*nF*tau_0))./(1i*(nF+(nF==0))); th0(1) = 0; th0(N_FFT/2+1) = tau_0;
F = @(h) fftshift(fft(ifftshift(h))/N_FFT);

g11 = zeros(2*N+1); g12 = g11; g21 = g11; g22 = g11;
r11 = zeros(1,2*N+1); r12 = r11; r21 = r11; r22 = r11; xi1 = r11;
for k =-N:N
    col = exp(1i*k*(t - tau_k)); ek = exp(1i*k*t);
    G = F(d_1G_1.*col + d_2G_11.*ek); g11(:,k+N+1) = G(N_FFT/2-N+1:N_FFT/2+N+1); r11(k+N+1) = sum(c0.*G);
    G = F((d_2G_12 - d_1G_1_ty).*ek); g12(:,k+N+1) = G(N_FFT/2-N+1:N_FFT/2+N+1); r12(k+N+1) = sum(c0.*G);
    G = F(d_1G_2.*col + d_2G_21.*ek); g21(:,k+N+1) = G(N_FFT/2-N+1:N_FFT/2+N+1); r21(k+N+1) = sum(c0.*G);
    G = F(-d_1G_2_ty.*ek);            g22(:,k+N+1) = G(N_FFT/2-N+1:N_FFT/2+N+1); r22(k+N+1) = sum(c0.*G);
    xi1(k+N+1) = sum(th0.*F(vp_circ_w.*ek));
end
rb = sum(c0.*F(d_3G_1)); xb = sum(th0.*F(v_circ_w));

% Fourier coefficients of d_3 calG_1 (bar v_4 of (124))
d_3G_1 = F(d_3G_1);
d_3G_1 = d_3G_1(N_FFT/2-N+1:N_FFT/2+N+1,1);

% d_beta xi
xi_beta = zeros(4);

xi_beta(1,2) = real(xb);
xi_beta(3,2) = d_3G_1(N+1);
xi_beta(4,1) = 1;

% L_1 (Section 2.4)
L1 = zeros(2*N+1);
for j = [-N:-1 1:N]
    L1(j+N+1,j+N+1) = 1/(1i*j) ;
end
for j = [-N:-1 1:N]
    L1(N+1,j+N+1) = (-1)^(j+1)/(1i*j) ;
end

% d_beta f and d_u f = I - L_1 d_u g of (111), with the zeroth rows of L_1 computed above
g_beta = zeros(2*(2*N+1),4);
g_beta(1:2*N+1,2) = -L1*d_3G_1; g_beta(N+1,2) = -rb;

f11 = eye(2*N+1) - L1*g11;
f12 = - L1*g12;
f21 = - L1*g21;
f22 = eye(2*N+1) - L1*g22;
f11(N+1,:) = (n1'==0) - r11; f12(N+1,:) = -r12;
f21(N+1,:) = -r21; f22(N+1,:) = (n1'==0) - r22;

df = [[f11 f12];[f21 f22]];

% d_u xi of (125)
xi_u = zeros(4,2*(2*N+1));

xi_u(1, 1:2*N+1) = omega*xi1;
xi_u(1, 2*N+2:end) = omega*ones(1,2*N+1)*v(real(sum(w.*exp(-1i*n1*tau_0))));
xi_u(2, 1:2*N+1) = ones(1,2*N+1) ; xi_u(2, 2*N+2:end) = zeros(1,2*N+1);
xi_u(3, 1:2*N+1) = g11(N+1,:) ; xi_u(3, 2*N+2:end) = g12(N+1,:);
xi_u(4, 1:2*N+1) = g21(N+1,:) ; xi_u(4, 2*N+2:end) = g22(N+1,:);

g_beta(N+1,3) = -1 ; g_beta(2*N+1 + N+1,4) = -1 ;

DF = [[xi_beta xi_u];[g_beta df]];

end
