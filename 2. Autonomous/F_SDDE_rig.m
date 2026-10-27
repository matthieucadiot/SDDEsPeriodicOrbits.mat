function F = F_SDDE_rig(x,par,C_rho,nu,N_FFT)
% Rigorous enclosure of the map F of (77) at bar x (autonomous SDDE (70)), in Fourier
% coefficients, for |n| <= N (the tail |n| > N of the zeroth row of L_1 G(bar U) is added in
% computation_bounds.m).
%   x     = (beta_1,beta_2,a) = (bar beta_1, bar beta_2, bar U), beta_1 = omega, a = (a_n)_{|n|<=N}
%   par   = (gamma,a_1,a_2,c_1,c_2,kappa_1,kappa_2)
%   C_rho = constant C of (34)-(36) for bar v_1 and bar v_2 of (79),  nu = bar nu
%   N_FFT = number of grid points (power of 2, N_FFT >= 4N)

omega = intval(x(1)); % omega = beta_1
u0_scalar = x(2); % bar beta_2 (not the phase-condition value u_0)
a = x(3:end);
gamma = par(1); a1 = par(2); a2 = par(3); c1 = par(4); c2 = par(5); k1 = par(6); k2 = par(7);

N = (length(a)-1)/2; n = intval((-N:N)');

if abs(log2(N_FFT)-ceil(log2(N_FFT)))>0 || N>N_FFT/4
    disp('Wrong size for N_FFT'), return
end

% zero padding to size N_FFT (Section 2.5.2)
ta = [intval(zeros(N_FFT/2-N,1));a;intval(zeros(N_FFT/2-N-1,1))];

% bar u(t_k) (inverse FFT)
u_k = ifftshift(verifyfft(fftshift(ta),-1))*intval(N_FFT);

% bar v_1(t_k), bar v_2(t_k), cf. (37)
t = (2*intval('pi')/N_FFT)*(-N_FFT/2:N_FFT/2-1)';
v1_at_t_k = intval(zeros(N_FFT,1));
v2_at_t_k = intval(zeros(N_FFT,1));

for k = 0:N_FFT-1
    v1_at_t_k(k+1) = real(sum(a.*exp(1i*n*(t(k+1)- a1/omega - c1*u_k(k+1)))));
    v2_at_t_k(k+1) = real(sum(a.*exp(1i*n*(t(k+1)- a2/omega - c2*u_k(k+1)))));
end

Vbar1 = fftshift(verifyfft(ifftshift(v1_at_t_k))/intval(N_FFT));
Vbar2 = fftshift(verifyfft(ifftshift(v2_at_t_k))/intval(N_FFT));

eps_alias = C_rho*(nu.^n + (1/nu).^n)/(nu^N_FFT - 1);   % aliasing bound (38)
eps_alias = infsup(-1,1)*eps_alias;

V1 = Vbar1(N_FFT/2-N+1:N_FFT/2+N+1,1) + eps_alias;
V2 = Vbar2(N_FFT/2-N+1:N_FFT/2+N+1,1) + eps_alias;

% G(bar U): Fourier coefficients of g(bar beta_1, bar u) of (75)
g = -omega*( k1*V1 + k2*V2 + gamma*a );

% L_1 (Section 2.4)
L1 = intval(zeros(2*N+1));
for j = [-N:-1 1:N]
    L1(j+N+1,j+N+1) = 1/(1i*intval(j)) ;
end
for j = [-N:-1 1:N]
    L1(N+1,j+N+1) = (-1)^(j+1)/(1i*intval(j)) ;
end

bar_g = g(N+1);   % M G(bar U)
u0 = [intval(zeros(N,1));u0_scalar;intval(zeros(N,1))];   % calB of (29)
f = a - u0 - L1*g;

F = [sum(intval(a)); bar_g; f];   % (xi_1, xi_2, f) of (77), with u_0 = 0

end
