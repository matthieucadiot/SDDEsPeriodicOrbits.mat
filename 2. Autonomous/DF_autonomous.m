function DF = DF_autonomous(x,par,N_FFT)
% Jacobian (81) in Fourier coefficients, |n| <= N, for the autonomous SDDE (70) (non-rigorous,
% used to build A^{<=N}).
%   x     = (beta_1,beta_2,a): beta_1 = omega, beta_2 (= u(-pi) at a zero of F),
%           a = (a_n)_{|n|<=N} = U
%   par   = (gamma,a_1,a_2,c_1,c_2,kappa_1,kappa_2)
%   N_FFT = number of grid points (power of 2, N_FFT >= 4N)

omega = x(1);
a = x(3:end);
gamma = par(1); a1 = par(2); a2 = par(3); c1 = par(4); c2 = par(5); k1 = par(6); k2 = par(7);

N = (length(a)-1)/2; n = (-N:N)';

ta = [zeros(N_FFT/2-N,1); a; zeros(N_FFT/2-N-1,1)];
u_k = ifftshift(ifft(fftshift(ta)))*N_FFT;
t = (2*pi/N_FFT)*(-N_FFT/2:N_FFT/2-1)';

% v_1,...,v_4 of (79) at the grid points
v1_grid = zeros(N_FFT,1); v2_grid = zeros(N_FFT,1);
v3_grid = zeros(N_FFT,1); v4_grid = zeros(N_FFT,1);
for k = 0:N_FFT-1
    delay1 = a1/omega + c1*u_k(k+1);
    delay2 = a2/omega + c2*u_k(k+1);
    v1_grid(k+1) = sum(a.*exp(1i*n*(t(k+1)-delay1))); % u(alpha_1)
    v2_grid(k+1) = sum(a.*exp(1i*n*(t(k+1)-delay2)));
    v3_grid(k+1) = sum((1i*n).*a.*exp(1i*n*(t(k+1)-delay1))); % u'(alpha_1)
    v4_grid(k+1) = sum((1i*n).*a.*exp(1i*n*(t(k+1)-delay2)));
end

% Fourier coefficients of v_1,...,v_4
V1 = fftshift(fft(ifftshift(v1_grid))/N_FFT);
V2 = fftshift(fft(ifftshift(v2_grid))/N_FFT);
V3 = fftshift(fft(ifftshift(v3_grid))/N_FFT);
V4 = fftshift(fft(ifftshift(v4_grid))/N_FFT);

% d_{beta_1} g, cf. (80)
v1 = V1(N_FFT/2-N+1:N_FFT/2+N+1);
v2 = V2(N_FFT/2-N+1:N_FFT/2+N+1);
v3 = V3(N_FFT/2-N+1:N_FFT/2+N+1);
v4 = V4(N_FFT/2-N+1:N_FFT/2+N+1);
g_omega = (-gamma*a - k1*v1 - k2*v2) - (k1*a1/omega)*v3 - (k2*a2/omega)*v4;

% d_u g, cf. (80): coefficients |n| <= 2N of v_3, v_4 for the multiplication matrices
v3_ext = V3(N_FFT/2-2*N+1:N_FFT/2+2*N+1);
v4_ext = V4(N_FFT/2-2*N+1:N_FFT/2+2*N+1);
L2_v3 = multiplicationmatrix(v3_ext,"Fou",N);
L2_v4 = multiplicationmatrix(v4_ext,"Fou",N);

% zeroth row of L_1*g_u: involves all the Fourier coefficients of the functions d_u g e_k
% (e_k = e^{ikt}), computed on the grid
nF = (-N_FFT/2:N_FFT/2-1)';
c0 = (-1).^(nF+1)./(1i*(nF+(nF==0))); c0([1 N_FFT/2+1]) = 0;
L1g_u0 = zeros(1,2*N+1);

L3_alpha1 = zeros(2*N+1); L3_alpha2 = zeros(2*N+1);
for k =-N:N
    c1_grid = exp(1i*k*(t - a1/omega - c1*u_k));
    c2_grid = exp(1i*k*(t - a2/omega - c2*u_k));
    C1 = fftshift(fft(ifftshift(c1_grid))/N_FFT);
    C2 = fftshift(fft(ifftshift(c2_grid))/N_FFT);
    L3_alpha1(:,k+N+1) = C1(N_FFT/2-N+1:N_FFT/2+N+1);
    L3_alpha2(:,k+N+1) = C2(N_FFT/2-N+1:N_FFT/2+N+1);
    ek = exp(1i*k*t);
    g = -omega*(gamma*ek + k1*(c1_grid - c1*v3_grid.*ek) + k2*(c2_grid - c2*v4_grid.*ek));
    L1g_u0(k+N+1) = sum(c0.*fftshift(fft(ifftshift(g))/N_FFT));
end

Id = eye(2*N+1);
g_u = -omega*gamma*Id - omega*k1*(L3_alpha1 - c1*L2_v3) - omega*k2*(L3_alpha2 - c2*L2_v4);

% L_1 (Section 2.4) and DF of (81)
L1 = zeros(2*N+1);
for j = [-N:-1 1:N]
    L1(j+N+1,j+N+1) = 1/(1i*j);
    L1(N+1,j+N+1) = (-1)^(j+1)/(1i*j);
end
e0 = zeros(2*N+1,1); e0(N+1) = 1;

L1g_omega = L1*g_omega;
L1g_omega(N+1) = sum(c0.*fftshift(fft(ifftshift(-gamma*u_k - k1*v1_grid - k2*v2_grid - (k1*a1/omega)*v3_grid - (k2*a2/omega)*v4_grid))/N_FFT));
Df_u = Id - L1*g_u; Df_u(N+1,:) = e0' - L1g_u0;

DF = [     0           0      ones(1,2*N+1) ;
      g_omega(N+1)     0       g_u(N+1,:)   ;
      -L1g_omega      -e0      Df_u  ];
end
