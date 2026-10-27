function F = F_SDDE_rig(x,eps,C_rho,nu,N_FFT)
% Rigorous enclosure of F(bar X), the map (29) for the forced toy model (4), for |n| <= N
% (the tail |n| > N of the zeroth row of L_1 G(bar U) is added in computation_bounds.m).
%   x     = (beta,a) = (bar beta, bar U), a = (a_n)_{|n|<=N}
%   eps   = (eps_1,eps_2)
%   C_rho = constant C of (34)-(36) for v = bar u(t-1-eps_2 bar u(t)),  nu = bar nu
%   N_FFT = number of grid points (power of 2, N_FFT >= 4N)

beta = x(1);
a = x(2:end);

N = (length(a)-1)/2; n = intval((-N:N)');

if abs(log2(N_FFT)-ceil(log2(N_FFT)))>0 || N>N_FFT/4
    disp('Wrong size for N_FFT'), return
end

sin_t = [intval(zeros(N-1,1));1i/intval(2);0;-1i/intval(2);intval(zeros(N-1,1))];
cos_t = [intval(zeros(N-1,1));1/intval(2);0;1/intval(2);intval(zeros(N-1,1))];

% zero padding to size N_FFT (Section 2.5.2)
ta = [intval(zeros(N_FFT/2-N,1));a;intval(zeros(N_FFT/2-N-1,1))];

% bar u(t_k) (inverse FFT)
u_k = ifftshift(verifyfft(fftshift(ta),-1))*intval(N_FFT);

% v(t_k) = bar u(t_k - 1 - eps_2 bar u(t_k)), cf. (37)
t = (2*intval('pi')/N_FFT)*(-N_FFT/2:N_FFT/2-1)';
v_at_t_k = intval(zeros(N_FFT,1));

for k = 0:N_FFT-1
    v_at_t_k(k+1) = real(sum(a.*exp(1i*n*(t(k+1)-1-eps(2)*u_k(k+1)))));
end

Vbar = fftshift(verifyfft(ifftshift(v_at_t_k))/intval(N_FFT));

eps_alias = C_rho*(nu.^n + (1/nu).^n)/(nu^N_FFT - 1);   % aliasing bound (38)
eps_alias = infsup(-1,1)*eps_alias;

V = Vbar(N_FFT/2-N+1:N_FFT/2+N+1,1) + eps_alias;

% G(bar U) of (28): Fourier coefficients of g(bar u) = phi - bar u - eps_1 v
G = cos_t + sin_t - a - eps(1)*V;

% L_1 (Section 2.4)
L1 = intval(zeros(2*N+1));
for j = [-N:-1 1:N]
    L1(j+N+1,j+N+1) = 1/(1i*intval(j)) ;
end
for j = [-N:-1 1:N]
    L1(N+1,j+N+1) = (-1)^(j+1)/(1i*intval(j)) ;
end

MG = G(N+1);   % M G(bar U)
calB = [intval(zeros(N,1));beta;intval(zeros(N,1))];   % calB of (29)
f = a - calB - L1*G;

F = [MG;f];

end
