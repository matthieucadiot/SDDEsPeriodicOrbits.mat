clc, clear, close all

% Computer-assisted proof of Theorem 14: periodic solution of the forced toy model (4)
% with (eps_1,eps_2) = (0.5,2). Requires INTLAB.
% x = (beta,U): beta (= u(-pi) at a zero of F, Lemma 2), U = Fourier coefficients of u;
% N_u = order of bar u, nu = bar nu.

N_u = 150; eps = [0.5;2]; N_FFT = 2^12; nu = 1.1;

a0 = [zeros(N_u-1,1);1i/2;0;-1i/2;zeros(N_u-1,1)]; % sin(t), the solution for eps_1 = 0

u0 = 0; % initial guess for beta
x = [u0;a0];

x = newton(x,eps,N_FFT);
x(1) = real(x(1)); x(2:end) = (x(2:end) + conj(flip(x(2:end))))/2; % real bar u: U_{-n} = conj(U_n)

a = x(2:end);
plot_periodic(a);

disp('Computing C_rho... ')
[C_rho, C_rho_der] = C_rho_contour_SDDE(a,nu,eps,N_FFT);
disp(['C_rho = ',num2str(C_rho)])

N = 660;    % truncation size N (Section 2.3)
N_Z1 = N;   % size of the block of DF inverted to build A^{<=N}

a = [zeros(N-N_u,1);a;zeros(N-N_u,1)]; x = [x(1);a];

% one constant per column k = -N,...,N of L_3 (63), hence computed after the padding
disp('Computing C_rho_col... ')
C_rho_col = C_rho_contour_operator(a,nu,eps,N_FFT);

[Y, Z1, Z2]  = computation_bounds(x,eps,N_FFT,C_rho,C_rho_der,C_rho_col,nu,N_Z1);
Y = intval(sup(Y)) ; Z1 = intval(sup(Z1)) ; Z2 = intval(sup(Z2));
r = abs((1 - Z1 - sqrt((1-Z1)^2 - 2*Y*Z2))/(Z2));   % smallest root of Y - (1-Z1) r + Z2 r^2/2

if 1 - Z1 - r*Z2 > 0
   if (1-Z1)^2 - 2*Y*Z2 > 0
        disp(['Proof successful (Theorem 3), smallest admissible radius r_0 = ',num2str(sup(r))])
   else
        disp('Proof failed : (1-Z1)^2 - 2*Y*Z2 <= 0')
   end
else
    disp('Proof failed : 1 - Z1 - r*Z2 <= 0')
end

% radius r_0 of Theorem 14 (Z_2 does not depend on r, so r_* is arbitrary)
r0 = intval('2.06e-7');
if sup(Y + Z1*r0 + Z2*r0^2/2 - r0) <= 0 && sup(Z1 + Z2*r0) < 1
    disp('Theorem 14 holds with r_0 = 2.06e-7')
else
    disp('Theorem 14 is NOT verified with r_0 = 2.06e-7')
end
