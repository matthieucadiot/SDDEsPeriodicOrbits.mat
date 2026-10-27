function F = F_SDDE(x,eps,N_FFT)
% Map F of (29) (Fourier representation of (16)) for the forced toy model (4), truncated to
% |n| <= N (non-rigorous, used by newton.m).
%   x     = (beta,a): beta (= u(-pi) at a zero, Lemma 2), a = (a_n)_{|n|<=N} Fourier coefficients
%           of u
%   eps   = (eps_1,eps_2)
%   N_FFT = number of grid points (power of 2, N_FFT >= 4N)

u0 = x(1);
a = x(2:end);

N = (length(a)-1)/2; n = (-N:N)';

if abs(log2(N_FFT)-ceil(log2(N_FFT)))>0 || N>N_FFT/4
    disp('Wrong size for N_FFT'), return
end

sin_t = [zeros(N-1,1);1i/2;0;-1i/2;zeros(N-1,1)];
cos_t = [zeros(N-1,1);1/2;0;1/2;zeros(N-1,1)];

% zero padding to size N_FFT (Section 2.5.2)
ta = [zeros(N_FFT/2-N,1);a;zeros(N_FFT/2-N-1,1)];

% u(t_k) (inverse FFT)
u_k = ifftshift(ifft(fftshift(ta)))*N_FFT;

% u(t_k - 1 - eps_2 u(t_k)), cf. (37)
t = (2*pi/N_FFT)*(-N_FFT/2:N_FFT/2-1)';
G_at_t_k = zeros(N_FFT,1);

for k = 0:N_FFT-1
    G_at_t_k(k+1) = real(sum(a.*exp(1i*n*(t(k+1)-1-eps(2)*u_k(k+1)))));
end

delayed_term_ext = fftshift(fft(ifftshift(G_at_t_k))/N_FFT);

% G(U) of (28): Fourier coefficients of g(u) = phi - u - eps_1 u(t-1-eps_2 u(t))
g = cos_t + sin_t - a - eps(1)*delayed_term_ext(N_FFT/2-N+1:N_FFT/2+N+1,1);

% L_1 (Section 2.4)
L1 = zeros(2*N+1);
for j = [-N:-1 1:N]
    L1(j+N+1,j+N+1) = 1/(1i*j) ;
end
for j = [-N:-1 1:N]
    L1(N+1,j+N+1) = (-1)^(j+1)/(1i*j) ;
end

bar_g = g(N+1);                    % M G(U)
u0 = [zeros(N,1);u0;zeros(N,1)];   % calB of (29)
f = a - u0 - L1*g;

F = [bar_g;f];

end
