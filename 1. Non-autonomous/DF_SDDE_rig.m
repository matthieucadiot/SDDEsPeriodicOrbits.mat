function DF = DF_SDDE_rig(x,eps,C_rho_der,C_rho_col,nu,N_FFT)
% Rigorous enclosure of the Galerkin matrix Pi^{<=N} DF(bar X) Pi^{<=N} of the Jacobian (17)
% for the forced toy model (4) (Section 2.6.5).
%   x         = (beta,a), a = (a_n)_{|n|<=N} = bar U
%   eps       = (eps_1,eps_2)
%   C_rho_der = constant C of (34) for bar u'(t-1-eps_2 bar u(t))
%   C_rho_col = constants C of (34) for the columns e^{ik(t-1-eps_2 bar u(t))} of L_3, |k| <= N
%   nu = bar nu,  N_FFT = number of grid points (power of 2, N_FFT >= 4N)

a = x(2:end) ;
N = (length(a)-1)/2; n = (-N:N)';

if abs(log2(N_FFT)-ceil(log2(N_FFT)))>0 || N>N_FFT/4
    disp('Wrong size for N_FFT'), return
end

% zero padding to size N_FFT (Section 2.5.2)
ta = [intval(zeros(N_FFT/2-N,1));a;intval(zeros(N_FFT/2-N-1,1))];
% bar u(t_k) (inverse FFT)
u_k = ifftshift(verifyfft(fftshift(ta),-1))*N_FFT;
% grid points
t = (2*intval('pi')/N_FFT)*(-N_FFT/2:N_FFT/2-1)';

% L_1 (Section 2.4)
L1 = intval(zeros(2*N+1));
for j = [-N:-1 1:N]
    L1(j+N+1,j+N+1) = 1/(1i*intval(j)) ;
end
for j = [-N:-1 1:N]
    L1(N+1,j+N+1) = (-1)^(j+1)/(1i*intval(j)) ;
end

% L_2 (Section 2.6.5): multiplication by v = bar u'(t-1-eps_2 bar u(t)), Toeplitz in V_n,
% |n| <= 2N, enclosed with the aliasing bound (38)
V = intval(zeros(N_FFT,1));
for k = 0:N_FFT-1
    V(k+1) = real(sum((1i*intval(n)).*a.*exp(1i*intval(n)*(t(k+1)-intval(1)-eps(2)*u_k(k+1)))));
end
Vbar = fftshift(verifyfft(ifftshift(V))/intval(N_FFT));
n2 = (-2*N:2*N)';
eps_alias = C_rho_der*(nu.^n2 + (1/nu).^n2)/(nu^N_FFT - 1);
eps_alias = infsup(-1,1)*eps_alias;

V = Vbar(N_FFT/2-2*N+1:N_FFT/2+2*N+1,1) + eps_alias;

L2 = multiplicationmatrix(V,"Fou", N);

% L_3 of (63), by columns, with the aliasing bound (38)
L3 = intval(zeros(2*N+1));

for k =-N:N
    col = exp(1i*intval(k)*(t - 1 - eps(2)*u_k));
    col_bar = fftshift(verifyfft(ifftshift(col))/intval(N_FFT));
    col_alias = C_rho_col(k+N+1)*(nu.^n + (1/nu).^n)/(nu^N_FFT - 1);
    col_alias = infsup(-1,1)*col_alias;

    L3(:,k+N+1) = col_bar(N_FFT/2-N+1:N_FFT/2+N+1,1) + col_alias;
end

DGu = -intval(eye(2*N+1)) + eps(1)*eps(2)*L2 - eps(1)*L3;   % d_u g(bar u), cf. (62)
DFu = intval(eye(2*N+1)) - L1*DGu;

% zeroth row of L_1*DGu: involves all the Fourier coefficients of the functions d_u g(bar u) e_k
% (e_k = e^{ikt}), enclosed on N0 points with the aliasing bound (38), the coefficients
% |n| >= N0/2 being bounded with (34)
N0 = 4*N_FFT; j0 = (-N0/2:N0/2-1)';
c0 = (-1).^(j0+1)./(1i*intval(j0+(j0==0))); c0([1 N0/2+1]) = 0;
alias0 = infsup(-1,1)*(nu.^(j0-N0) + nu.^(-j0-N0))/(1-nu^(-N0));
tail0 = infsup(-1,1)*2*nu^(1-N0/2)/((N0/2)*(nu-1));
t0 = (2*intval('pi')/N0)*j0;
u0 = ifftshift(verifyfft(fftshift([intval(zeros(N0/2-N,1));a;intval(zeros(N0/2-N-1,1))]),-1))*N0;
V0 = intval(zeros(N0,1));
for m = 1:N0
    V0(m) = real(sum((1i*intval(n)).*a.*exp(1i*intval(n)*(t0(m)-1-eps(2)*u0(m)))));
end
for k =-N:N
    % d_u g(bar u) e_k on the N0 points, and its constant C of (34)
    C = nu^abs(k)*(1 + abs(eps(1)*eps(2))*C_rho_der) + abs(eps(1))*C_rho_col(k+N+1);
    g = exp(1i*intval(k)*t0).*(eps(1)*eps(2)*V0 - 1) - eps(1)*exp(1i*intval(k)*(t0 - 1 - eps(2)*u0));
    G = fftshift(verifyfft(ifftshift(g))/intval(N0)) + C*alias0;
    DFu(N+1,k+N+1) = (k==0) - sum(c0.*G) + C*tail0;
end

e0 = zeros(2*N+1,1) ; e0(N+1) = 1; e0 = intval(e0);

DF = [0  DGu(N+1,:);
      -e0  DFu   ];

end
