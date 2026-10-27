clc, clear, close all

% Computer-assisted proof of Theorem 18: periodic solution of the autonomous SDDE (70) with
% gamma = 4.75, a_1 = 1.3, a_2 = 6, c_1 = c_2 = 1, kappa_1 = 3.4, kappa_2 = 3. Requires INTLAB.
% x = (beta_1,beta_2,U): beta_1 = omega (period 2*pi*beta_1), beta_2 (= u(-pi) at a zero of F),
% U = Fourier coefficients of u (after the rescalings (71)-(72)); u_0 = 0 in the phase
% condition (76).
% par = (gamma,a_1,a_2,c_1,c_2,kappa_1,kappa_2), N_u = order of bar u, nu = bar nu.

load sol1   % x = bar x (N_u = 30), par and N_FFT = 2^12

% enclosures of the decimal parameters (e.g. a_1 = 1.3, kappa_1 = 3.4)
parI = intval(par); for i = 1:numel(par), parI(i) = intval(mat2str(par(i))); end; par = parI;
% real bar x: real (beta_1,beta_2) and U_{-n} = conj(U_n)
x(1:2) = real(x(1:2)); x(3:end) = (x(3:end) + conj(flip(x(3:end))))/2;

nu = 1.1;

a = x(3:end); N_u = (length(a)-1)/2;

[C_rho1, C_rho2, C_rho3, C_rho4] = C_rho_contour_SDDE(x,par,nu,N_FFT);
C_rho = [C_rho1; C_rho2; C_rho3; C_rho4];

N = 800;   % truncation size N (Section 2.3)
disp(['N = ',num2str(N)])

a = [zeros(N-N_u,1);a;zeros(N-N_u,1)]; x = [x(1);x(2);a];

% one constant per column k = -N,...,N of the delay operators, hence computed after the padding
[C_rho_col1, C_rho_col2] = C_rho_contour_operator(x,par,nu,N_FFT);
C_rho_col = [C_rho_col1(:); C_rho_col2(:)];

r_star = intval(1e-3);   % r_* of Theorem 3

N_Z1 = N;   % size of the block of DF inverted to build A^{<=N}

[Y, Z1, Z2]  = computation_bounds(x,par,r_star,N_FFT,C_rho,C_rho_col,nu,N_Z1);
Y = intval(sup(Y)) ; Z1 = intval(sup(Z1)) ; Z2 = intval(sup(Z2));
r = abs((1 - Z1 - sqrt((1-Z1)^2 - 2*Y*Z2))/(Z2));   % smallest root of Y - (1-Z1) r + Z2 r^2/2

if 1 - Z1 - r*Z2 > 0
   if (1-Z1)^2 - 2*Y*Z2 > 0
        if sup(r) <= inf(r_star)
            disp(['Proof successful (Theorem 3), smallest admissible radius r_0 = ',num2str(sup(r))])
        else
            disp('Proof failed : r > r_star (Z2 is only valid on the ball of radius r_star)')
        end
   else
        disp('Proof failed : (1-Z1)^2 - 2*Y*Z2 <= 0')
   end
else
    disp('Proof failed : 1 - Z1 - r*Z2 <= 0')
end

% radius r_0 of Theorem 18
r0 = intval('9.7e-11');
if sup(Y + Z1*r0 + Z2*r0^2/2 - r0) <= 0 && sup(Z1 + Z2*r0) < 1 && sup(r0) <= inf(r_star)
    disp('Theorem 18 holds with r_0 = 9.7e-11')
else
    disp('Theorem 18 is NOT verified with r_0 = 9.7e-11')
end
