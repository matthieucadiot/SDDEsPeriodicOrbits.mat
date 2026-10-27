clc, clear, close all

% Computer-assisted proof of Theorem 28: periodic solution of the threshold SDDE (100) for the
% parameters of Figure 3. Requires INTLAB.
% x = (beta,u): beta = (eta,omega,beta_3,beta_4) ((beta_3,beta_4) = u(-pi) at a zero of F) and
% u = (u_1,u_2) = (w,tau) in Fourier coefficients (after the rescalings (106)-(108));
% par = (a,b,mu,gamma,m,n,q^-,q^+,theta_q,v^-,v^+,theta_v) of (100)-(101); N_u = order of bar u,
% nu = bar nu, N2 = N_2 of Lemma 26.

sol = 'pt2';            % approximate solution bar x of Figure 3 (N_u = 70)
nu  = 1.05;             % bar nu = e^{bar rho}, S_{bar rho} = analyticity strip of Section 4.1.1
N_FFT = 2^14;           % number of grid points (power of 2, N_FFT >= 4*N2 = 8N)
r_star = intval(1e-3);  % r_* of Theorem 3

load(sol)               % x and par

% enclosures of the decimal parameters (e.g. mu = 0.3, gamma = 0.13, v^+ = 0.2)
parI = intval(par); for i = 1:numel(par), parI(i) = intval(mat2str(par(i))); end; par = parI;

beta = x(1:4);
u = x(5:end); N_u = (length(u)-2)/4;
w   = u(1:2*N_u+1);
tau = u(2*N_u+2:4*N_u+2);
% real bar x: real beta and Hermitian Fourier coefficients
beta = real(beta); w = (w + conj(flip(w)))/2; tau = (tau + conj(flip(tau)))/2;
x = [beta; w; tau];

% analyticity on S_{bar rho}, bar rho = log(nu) (Lemma 23)
disp('Checking analyticity on the strip ...')
if ~check_analyticity(x,par,nu,N_FFT)
    disp('Analyticity check failed : decrease nu'), return
end

% constants C of (34) for bar g_1, bar g_2 and v(bar w)
disp('Computing the constants C of (34) ...')
[C_g1, C_g2, C_vw] = C_rho_contour_SDDE(x,par,nu,N_FFT);
C_rho = max([C_g1, C_g2, C_vw]);
disp(['  C_g1 = ',num2str(sup(C_g1)),',  C_g2 = ',num2str(sup(C_g2)), ...
      ',  C_vw = ',num2str(sup(C_vw))])

N = 1000;   % truncation size N (Section 2.3)
w2   = [zeros(N-N_u,1); w;   zeros(N-N_u,1)];
tau2 = [zeros(N-N_u,1); tau; zeros(N-N_u,1)];
x = [beta; w2; tau2];

% constants C of (34) for DF(bar x), including the columns |k| <= N_2 of L_3(bar u_2), and
% lambda of (126)
N2 = 2*N;
disp('Computing operator strip bounds ...')
[C_col, C_DF, Cv] = C_rho_contour_operator(x,par,nu,N_FFT,N2);

N_Z1 = N;   % size of the block of DF inverted to build A^{<=N}

[Y, Z1, Z2] = computation_bounds(x,par,r_star,N_FFT,C_g1,C_g2,C_vw,C_col,C_DF,Cv,nu,N_Z1,N2);
Y = intval(sup(Y)); Z1 = intval(sup(Z1)); Z2 = intval(sup(Z2));

disp('-------------------------------------------------')
disp(['Y  <= ',num2str(sup(Y))])
disp(['Z1 <= ',num2str(sup(Z1))])
disp(['Z2(r_*) <= ',num2str(sup(Z2))])

% smallest root of Y - (1-Z1) r + Z2 r^2/2 (Theorem 3)
r = abs((1 - Z1 - sqrt((1-Z1)^2 - 2*Y*Z2))/Z2);
if (1 - Z1 - r*Z2 > 0) && ((1-Z1)^2 - 2*Y*Z2 > 0) && (sup(r) <= inf(r_star))
    disp(['Proof successful (Theorem 3), smallest admissible radius r_0 = ',num2str(sup(r))])
    disp(['   (r_0 <= r_* = ',num2str(sup(r_star)),')'])
else
    disp('Proof failed : radii polynomial has no admissible root r0 <= r_*.')
end

% radius r_0 of Theorem 28
r0 = intval('2e-9');
if sup(Y + Z1*r0 + Z2*r0^2/2 - r0) <= 0 && sup(Z1 + Z2*r0) < 1 && sup(r0) <= inf(r_star)
    disp('Theorem 28 holds with r_0 = 2e-9')
else
    disp('Theorem 28 is NOT verified with r_0 = 2e-9')
end
