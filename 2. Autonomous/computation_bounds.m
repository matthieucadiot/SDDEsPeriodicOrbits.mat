function [Y, Z1, Z2] = computation_bounds(x,par,r_star,N_FFT,C_rho,C_rho_col,nu,N_Z1)
% Bounds Y, Z_1 and Z_2(r_*) of Theorem 3 for the autonomous SDDE (70) (Section 3.1):
%   Y (Lemma 4), Z_1 = Z_{1,0} + Z_{1,1} + Z_{1,2} (Lemma 9 with (83), and Lemma 16) and
%   Z_2(r_*) (Lemma 17), the operator norms on H_0 being computed with Lemmas 11 and 13.
%   x = bar X = (bar beta_1, bar beta_2, bar U) (order N),
%   par = (gamma,a_1,a_2,c_1,c_2,kappa_1,kappa_2),
%   r_star = r_*,  nu = bar nu,  C_rho = constants C of (34) for bar v_1,...,bar v_4 of (79),
%   C_rho_col = constants C of (34) for the columns of the delay operators
%   (C_rho_contour_operator.m),
%   N_Z1 = size of the block of DF inverted to build A^{<=N} (default N).
% Notation: om = omega = bar beta_1, k1, k2 = kappa_1, kappa_2, A = matrix of A^{<=N} (in the
%   comments, A = A^{<=N} + pi^{>N} is the operator of (24)), Ddag = D^dagger,
%   e = (0,0,1) of Lemma 11, thetaN = theta_N.

    C_rho1 = C_rho(1) ; C_rho2 = C_rho(2) ; C_rho3 = C_rho(3) ; C_rho4 = C_rho(4) ;
    C_rho = max(C_rho) ; C_rho_col = max(C_rho_col);

    % only the block |n| <= N_Z1 of DF is inverted, A^{<=N} being the identity for N_Z1 < |n| <= N
    % (x and all the tail bounds use the truncation size N)
    Nx = (length(x(3:end))-1)/2;
    if nargin < 8 || isempty(N_Z1), N_Z1 = Nx; end
    N_Z1 = min(N_Z1,Nx);
    ind = Nx-N_Z1+1 : Nx+N_Z1+1;   io = [1, 2, 2+ind];

    % approximate inverse A^{<=N} (Section 2.3), computed in floating point
    DF = DF_autonomous(mid([x(1); x(2); x(2+ind)]),mid(par),N_FFT);
    A = intval(eye(2*Nx+3));
    A(io,io) = intval(inv(DF));

    omega = intval(x(1)); % omega = bar beta_1
    par = intval(par); nu = intval(nu);

    a = x(3:end); % bar U
    gamma = par(1); a1 = par(2); a2 = par(3); c1 = par(4); c2 = par(5); k1 = par(6); k2 = par(7);

    N = (length(a)-1)/2; n = intval((-N:N)');

    %% Y bound (Lemma 4)

    F = F_SDDE_rig(x,par,C_rho,nu,N_FFT);   % F(bar X), |n| <= N

    % tail |k| > N of the zeroth row of L_1 G(bar U), bounded with (34)
    F(N+3) = F(N+3) + infsup(-1,1)*2*C_rho*nu/((N+1)*(nu-1)*nu^(N+1))*omega*(abs(k1) + abs(k2));

    AF = A*F ;

    Y = max([abs(AF(1)) abs(AF(2)) norm(AF(3:end),1) norm(n.*AF(3:end),1)]) ;
    % |n| > N: bar W_n = -G_n/(in), with G_n bounded with (34) (as in Section 2.6.1)
    Y = Y +  2*abs(C_rho)*omega*(abs(k1) + abs(k2))*nu/((nu-1)*nu^(N+1)) ;
    Y = abs(Y);
    disp(['Y <= ',num2str(sup(Y))])

    %% Z2 bound (Lemma 17, with r = r_*)

    L1 = intval(zeros(2*N+1));
    for j = [-N:-1 1:N]
        L1(j+N+1,j+N+1) = 1/intval(1i*j) ;
    end
    for j = [-N:-1 1:N]
        L1(N+1,j+N+1) = (-1)^(j+1)/intval(1i*j) ;
    end
    Ddag = [infsup(1,1)  intval(zeros(1,2*N+2)) ;   % D^dagger = diag(1, 1, L_1)
         infsup(0,0)  infsup(1,1)  intval(zeros(1,2*N+1))
          intval(zeros(2*N+1,2)) L1] ;

    D = diag([infsup(1,1) ; infsup(1,1) ; 1i*n]);

    % ||A^{<=N} D^dagger|| and ||D A^{<=N} D^dagger||: first inequality of (60) (Lemma 11)
    thetaN = kolmogorov_theta_N(N);
    if sup(thetaN) > 1, error('theta_N > 1 : Lemma 9 does not apply'), end
    e = intval(zeros(2*N+3,1)); e(N+3) = 1;
    DADdag = D*A*Ddag;
    normAD  = comp_mat_norm(A*Ddag,N_FFT) + thetaN*norm(A*e,1);
    normDAD = comp_mat_norm(DADdag,N_FFT) + thetaN*norm(D*(A*e),1);

    % max{||A D^dagger||, ||D A D^dagger||} (second part of Lemma 11), PiJ = pi^{<=N} J
    PiJ = diag([intval(0) ; intval(0) ; intval(ones(2*N+1,1))]);
    norm_ADdag = max(normAD + thetaN, comp_mat_norm(DADdag-PiJ,N_FFT) + thetaN*norm(D*(A*e),1) + 1);
    norm_der2_u = norm(n.^2.*a,1) ;   % >= ||bar u''||_{C^0}
    norm_der_u = norm(n.*a,1) ;       % >= ||bar u'||_{C^0}

    om = intval(omega) ; r = r_star ;
    if sup(r) >= inf(abs(om)), error('r_* >= |bar beta_1| : Lemma 17 does not apply'), end

    % delta_j and ell_j of (85)
    delta1 = abs(a1)/(abs(om)*(abs(om) - r)) ; ell1 = delta1 + abs(c1) ;
    delta2 = abs(a2)/(abs(om)*(abs(om) - r)) ; ell2 = delta2 + abs(c2) ;

    % hat Z_2^beta(r_*) and hat Z_2^u(r_*) of (86)
    hatZ2_beta = abs(gamma) + abs(k1)*(1 + ell1*norm_der_u + delta1*(norm_der_u + r) + abs(om)*delta1*(1 + ell1*norm_der2_u)) ...
                         + abs(k2)*(1 + ell2*norm_der_u + delta2*(norm_der_u + r) + abs(om)*delta2*(1 + ell2*norm_der2_u)) ;
    hatZ2_u    = abs(gamma) + abs(k1)*(1 + abs(om)*ell1 + abs(c1)*(norm_der_u + r + abs(om)*(1 + ell1*norm_der2_u))) ...
                         + abs(k2)*(1 + abs(om)*ell2 + abs(c2)*(norm_der_u + r + abs(om)*(1 + ell2*norm_der2_u))) ;

    % Z_2(r_*), with hat Z_2(r) = hat Z_2^beta(r_*) + hat Z_2^u(r_*) for r <= r_*
    % (proof of Theorem 18)
    Z2 = abs(norm_ADdag*(hatZ2_beta + hatZ2_u));

    disp(['Z2(r_*) <= ',num2str(sup(Z2))])

    %% Z1 bound (Lemma 9 with (83), and Lemma 16)

    [DF, V3, V4] =  DF_SDDE_rig(x,par,[C_rho1;C_rho2;C_rho3;C_rho4],C_rho_col,nu,N_FFT);

    % tails of the geometric decay: sum_{|n|>N} nu^{-|n|} and sum_{|n|>N} |n| nu^{-|n|}
    tail_0 = 2*nu^(-N-1)/(1-1/nu);
    tail_1 = 2*(nu^(-N-1)*(N+1-N/nu))/(1-1/nu)^2;

    % upper bounds of ||bar v_j||_{C^0} and ||bar v_j'||_{C^0}, j = 3,4, using (34)
    norm_v3 = norm(V3,1) + C_rho3*tail_0 ; norm_dv3 = norm(n.*V3,1) + C_rho3*tail_1 ;
    norm_v4 = norm(V4,1) + C_rho4*tail_0 ; norm_dv4 = norm(n.*V4,1) + C_rho4*tail_1 ;

    % ||d_u g(bar x)||_{B(C^0)} and sup_{||h||_{C^1}=1} ||d_t[d_u g(bar x)h]||_{C^0} (Lemma 16)
    norm_dug   = abs(om*gamma) + abs(om*k1)*(1 + abs(c1)*norm_v3) + abs(om*k2)*(1 + abs(c2)*norm_v4);
    norm_dtdug = abs(om*gamma) + abs(om*k1)*(1 + abs(c1)*norm_der_u + abs(c1)*(norm_v3 + norm_dv3)) ...
                               + abs(om*k2)*(1 + abs(c2)*norm_der_u + abs(c2)*(norm_v4 + norm_dv4));

    % ||Pi^{>N} F(d_{beta_1} g(bar x))||_{ell^1} <= C_beta*tail_0 in (83), with d_{beta_1} g(bar x)
    % of (80)
    C_beta = abs(k1)*C_rho1 + abs(k2)*C_rho2 + abs(k1*a1/om)*C_rho3 + abs(k2*a2/om)*C_rho4;

    Z11 = max(normAD,normDAD)*norm_dug*thetaN;
    Z12 = C_beta*tail_0 + norm_dtdug*thetaN;

    % Z_{1,0} of (54), with the second inequality of (60) for B = I - A DF(bar x) and
    % B = D(I - A DF(bar x))
    defect = intval(eye(2*N+3)) - A*DF;   % I - A^{<=N} DF(bar x) on the Galerkin space
    Z10 = max(comp_mat_norm(defect*Ddag,N_FFT) + (1+thetaN)*norm(defect*e,1), ...
              comp_mat_norm(D*defect*Ddag,N_FFT) + (1+thetaN)*norm(D*(defect*e),1));

    Z1 = abs(Z10 + Z11 + Z12);
    disp(['Z1 <= ',num2str(sup(Z1))])

end
