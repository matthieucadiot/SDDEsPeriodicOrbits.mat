function [DF, V3, V4] = DF_SDDE_rig(x,par,C_rho,C_rho_col,nu,N_FFT)
% Rigorous enclosure of the Galerkin matrix of the Jacobian (81) at bar x (autonomous SDDE (70)),
% in Fourier coefficients, |n| <= N. V3, V4 = Fourier coefficients of bar v_3, bar v_4 of (79).
%   x         = (beta_1,beta_2,a) = (bar beta_1, bar beta_2, bar U), beta_1 = omega,
%               a = (a_n)_{|n|<=N}
%   par       = (gamma,a_1,a_2,c_1,c_2,kappa_1,kappa_2)
%   C_rho     = [C_rho1;...;C_rho4], constants C of (34)-(36) for bar v_1,...,bar v_4
%   C_rho_col = constant C of (34) for the columns e^{ik alpha_j} of L31, L32
%   nu = bar nu,  N_FFT = number of grid points (power of 2, N_FFT >= 4N)

omega = intval(x(1));
a = x(3:end);
gamma = par(1); a1 = par(2); a2 = par(3); c1 = par(4); c2 = par(5); k1 = par(6); k2 = par(7);

N = (length(a)-1)/2; n = intval((-N:N)');

if abs(log2(N_FFT)-ceil(log2(N_FFT)))>0 || N>N_FFT/4
    disp('Wrong size for N_FFT'), return
end

% zero padding to size N_FFT (Section 2.5.2)
ta = [intval(zeros(N_FFT/2-N,1));a;intval(zeros(N_FFT/2-N-1,1))];
% bar u(t_k) (inverse FFT)
u_k = ifftshift(verifyfft(fftshift(ta),-1))*N_FFT;
% bar v_1(t_k),...,bar v_4(t_k), cf. (37)
t = (2*intval('pi')/N_FFT)*(-N_FFT/2:N_FFT/2-1)';

v1_at_t_k = intval(zeros(N_FFT,1));
v2_at_t_k = intval(zeros(N_FFT,1));
v3_at_t_k = intval(zeros(N_FFT,1));
v4_at_t_k = intval(zeros(N_FFT,1));

for k = 0:N_FFT-1
    v1_at_t_k(k+1) = real(sum(a.*exp(1i*n*(t(k+1)- a1/omega - c1*u_k(k+1)))));
    v2_at_t_k(k+1) = real(sum(a.*exp(1i*n*(t(k+1)- a2/omega - c2*u_k(k+1)))));
    v3_at_t_k(k+1) = real(sum((1i*n).*a.*exp(1i*n*(t(k+1)- a1/omega - c1*u_k(k+1)))));
    v4_at_t_k(k+1) = real(sum((1i*n).*a.*exp(1i*n*(t(k+1)- a2/omega - c2*u_k(k+1)))));
end

Vbar1 = fftshift(verifyfft(ifftshift(v1_at_t_k))/intval(N_FFT));
Vbar2 = fftshift(verifyfft(ifftshift(v2_at_t_k))/intval(N_FFT));
Vbar3 = fftshift(verifyfft(ifftshift(v3_at_t_k))/intval(N_FFT));
Vbar4 = fftshift(verifyfft(ifftshift(v4_at_t_k))/intval(N_FFT));

eps_alias = infsup(-1,1)*(nu.^n + (1/nu).^n)/(nu^N_FFT - 1);   % aliasing bound (38) for C = 1

V1 = Vbar1(N_FFT/2-N+1:N_FFT/2+N+1,1) + C_rho(1)*eps_alias;
V2 = Vbar2(N_FFT/2-N+1:N_FFT/2+N+1,1) + C_rho(2)*eps_alias;
V3 = Vbar3(N_FFT/2-N+1:N_FFT/2+N+1,1) + C_rho(3)*eps_alias;
V4 = Vbar4(N_FFT/2-N+1:N_FFT/2+N+1,1) + C_rho(4)*eps_alias;
% coefficients |n| <= 2N of bar v_3, bar v_4, for the multiplication operators
n2 = intval((-2*N:2*N)'); eps_alias2 = infsup(-1,1)*(nu.^n2 + (1/nu).^n2)/(nu^N_FFT - 1);
V3_ext = Vbar3(N_FFT/2-2*N+1:N_FFT/2+2*N+1,1) + C_rho(3)*eps_alias2;
V4_ext = Vbar4(N_FFT/2-2*N+1:N_FFT/2+2*N+1,1) + C_rho(4)*eps_alias2;

% L_1 (Section 2.4)
L1 = intval(zeros(2*N+1));
for j = [-N:-1 1:N]
    L1(j+N+1,j+N+1) = 1/(1i*intval(j)) ;
end
for j = [-N:-1 1:N]
    L1(N+1,j+N+1) = (-1)^(j+1)/(1i*intval(j)) ;
end

% LV3, LV4: multiplication by bar v_3, bar v_4 (as L_2 in Section 2.6.5)
LV3 = multiplicationmatrix(V3_ext,"Fou", N);
LV4 = multiplicationmatrix(V4_ext,"Fou", N);

% zeroth row of L_1*Du_g: involves all the Fourier coefficients of the functions d_u g(bar x) e_k
% (e_k = e^{ikt}), enclosed on the grid with the aliasing bound (38), the coefficients
% |n| >= N_FFT/2 being bounded with (34)
nF = (-N_FFT/2:N_FFT/2-1)';
c0 = (-1).^(nF+1)./(1i*intval(nF+(nF==0))); c0([1 N_FFT/2+1]) = 0;
alias0 = infsup(-1,1)*(nu.^(nF-N_FFT) + nu.^(-nF-N_FFT))/(1-nu^(-N_FFT));
tail0 = infsup(-1,1)*2*nu^(1-N_FFT/2)/((N_FFT/2)*(nu-1));
L1Du_g0 = intval(zeros(1,2*N+1));

% L31, L32: (n,k) -> F(e^{ik alpha_j(bar beta_1,bar u)})_n (as L_3 in (63)), by columns
L31 = intval(zeros(2*N+1));
L32 = intval(zeros(2*N+1));

for k =-N:N
    col1 = exp(1i*intval(k)*(t - a1/omega - c1*u_k));
    col2 = exp(1i*intval(k)*(t - a2/omega - c2*u_k));
    col_bar1 = fftshift(verifyfft(ifftshift(col1))/intval(N_FFT));
    col_bar2 = fftshift(verifyfft(ifftshift(col2))/intval(N_FFT));
    col_alias = C_rho_col*(nu.^n + (1/nu).^n)/(nu^N_FFT - 1);
    col_alias = infsup(-1,1)*col_alias;

    L31(:,k+N+1) = col_bar1(N_FFT/2-N+1:N_FFT/2+N+1,1) + col_alias;
    L32(:,k+N+1) = col_bar2(N_FFT/2-N+1:N_FFT/2+N+1,1) + col_alias;

    % d_u g(bar x) e_k on the grid, and its constant C of (34)
    ek = exp(1i*intval(k)*t);
    g = -omega*(gamma*ek + k1*(col1 - c1*v3_at_t_k.*ek) + k2*(col2 - c2*v4_at_t_k.*ek));
    C = abs(omega)*((abs(k1)+abs(k2))*C_rho_col + (abs(gamma) + abs(k1*c1)*C_rho(3) + abs(k2*c2)*C_rho(4))*nu^abs(k));
    L1Du_g0(k+N+1) = sum(c0.*(fftshift(verifyfft(ifftshift(g))/intval(N_FFT)) + C*alias0)) + C*tail0;
end

% d_u g(bar x) and d_{beta_1} g(bar x) of (80)
Du_g = -omega*gamma*intval(eye(2*N+1)) - omega*k1*(L31 - c1*LV3) - omega*k2*(L32 - c2*LV4);
Dom_g = (-gamma*a - k1*V1 - k2*V2) - (k1*a1/omega)*V3 - (k2*a2/omega)*V4 ;

Du_f = intval(eye(2*N+1)) - L1*Du_g;

e0 = zeros(2*N+1,1) ; e0(N+1) = 1; e0 = intval(e0);

Du_f(N+1,:) = e0.' - L1Du_g0;
% zeroth row of L_1*Dom_g, as above; C_beta = C of (34) for d_{beta_1} g(bar x) + gamma bar u
% (bar u has order N, so the term gamma bar u has no aliasing error and no tail)
C_beta = abs(k1)*C_rho(1) + abs(k2)*C_rho(2) + abs(k1*a1/omega)*C_rho(3) + abs(k2*a2/omega)*C_rho(4);
dom_g = -gamma*u_k - k1*v1_at_t_k - k2*v2_at_t_k - (k1*a1/omega)*v3_at_t_k - (k2*a2/omega)*v4_at_t_k;
L1Dom_g = L1*Dom_g;
L1Dom_g(N+1) = sum(c0.*(fftshift(verifyfft(ifftshift(dom_g))/intval(N_FFT)) + C_beta*alias0)) + C_beta*tail0;

% DF(bar x) of (81)
DF = [     intval(0)           intval(0)      intval(ones(1,2*N+1)) ;
      Dom_g(N+1)     intval(0)       Du_g(N+1,:)   ;
      -L1Dom_g    -e0      Du_f  ];

end
