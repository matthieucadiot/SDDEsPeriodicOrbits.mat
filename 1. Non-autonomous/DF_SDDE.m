function DF = DF_SDDE(x,eps,N_FFT)
% Jacobian (17) in Fourier coefficients, truncated to |n| <= N, for the forced toy model (4)
% (non-rigorous, used by newton.m and to build A^{<=N}).
%   x     = (beta,a): beta (= u(-pi) at a zero, Lemma 2), a = (a_n)_{|n|<=N} Fourier coefficients
%           of u
%   eps   = (eps_1,eps_2)
%   N_FFT = number of grid points (power of 2, N_FFT >= 4N)

a = x(2:end) ;
N = (length(a)-1)/2; n = (-N:N)';

if abs(log2(N_FFT)-ceil(log2(N_FFT)))>0 || N>N_FFT/4
    disp('Wrong size for N_FFT'), return
end

% zero padding to size N_FFT (Section 2.5.2)
ta = [zeros(N_FFT/2-N,1);a;zeros(N_FFT/2-N-1,1)];
% u(t_k) (inverse FFT)
u_k = ifftshift(ifft(fftshift(ta)))*N_FFT;
% grid points
t = (2*pi/N_FFT)*(-N_FFT/2:N_FFT/2-1)';

% L_1 (Section 2.4)
L1 = zeros(2*N+1);
for j = [-N:-1 1:N]
    L1(j+N+1,j+N+1) = 1/(1i*j) ;
end
for j = [-N:-1 1:N]
    L1(N+1,j+N+1) = (-1)^(j+1)/(1i*j) ;
end

% L_2 (Section 2.6.5): multiplication by v = u'(t-1-eps_2 u(t)), Toeplitz in V_n, |n| <= 2N
V = zeros(N_FFT,1);
for k = 0:N_FFT-1
    V(k+1) = real(sum((1i*n).*a.*exp(1i*n*(t(k+1)-1-eps(2)*u_k(k+1)))));
end
V = fftshift(fft(ifftshift(V))/N_FFT);
V = V(N_FFT/2-2*N+1:N_FFT/2+2*N+1,1) ;

L2 = multiplicationmatrix(V,"Fou",N);

% L_3 of (63), by columns
L3 = zeros(2*N+1);

for k =-N:N
    col = exp(1i*k*(t - 1 - eps(2)*u_k));
    col = fftshift(fft(ifftshift(col))/N_FFT);
    L3(:,k+N+1) = col(N_FFT/2-N+1:N_FFT/2+N+1,1);
end

DGu = -eye(2*N+1) + eps(1)*eps(2)*L2 - eps(1)*L3;   % d_u g(u), cf. (62)
DFu = eye(2*N+1) - L1*DGu;

% zeroth row of L_1*DGu: involves all the Fourier coefficients of the functions d_u g(u) e_k
% (e_k = e^{ikt}), computed on N0 points
N0 = 4*N_FFT; j0 = (-N0/2:N0/2-1)';
c0 = (-1).^(j0+1)./(1i*(j0+(j0==0))); c0([1 N0/2+1]) = 0;
t0 = (2*pi/N0)*j0;
u0 = ifftshift(ifft(fftshift([zeros(N0/2-N,1);a;zeros(N0/2-N-1,1)])))*N0;
V0 = zeros(N0,1);
for m = 1:N0
    V0(m) = real(sum((1i*n).*a.*exp(1i*n*(t0(m)-1-eps(2)*u0(m)))));
end
for k =-N:N
    g = exp(1i*k*t0).*(eps(1)*eps(2)*V0 - 1) - eps(1)*exp(1i*k*(t0 - 1 - eps(2)*u0));
    DFu(N+1,k+N+1) = (k==0) - sum(c0.*fftshift(fft(ifftshift(g))/N0));
end

e0 = zeros(2*N+1,1) ; e0(N+1) = 1;

DF = [0  DGu(N+1,:);
      -e0  DFu   ];

end
