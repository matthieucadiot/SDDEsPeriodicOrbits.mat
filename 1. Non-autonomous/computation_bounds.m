function [Y, Z1, Z2] = computation_bounds(x,eps,N_FFT,C_rho,C_rho_der,C_rho_col,nu,N_Z1)
% Bounds Y, Z_1 and Z_2 of Theorem 3 for the forced toy model (4) (Section 2.6):
%   Y (Lemma 4), Z_2 (Lemma 5) and Z_1 = Z_{1,0} + Z_{1,1} + Z_{1,2} (Lemmas 9 and 10), the
%   operator norms on H_0 being computed with Lemmas 11 and 13.
%   x = bar X = (bar beta, bar U) (order N),  eps = (eps_1,eps_2),  nu = bar nu,
%   C_rho, C_rho_der, C_rho_col = constants C of (34) (C_rho_contour_SDDE.m and
%   C_rho_contour_operator.m),
%   N_Z1 = size of the block of DF inverted to build A^{<=N} (default N).
% Notation: A = matrix of A^{<=N} (in the comments, A = A^{<=N} + pi^{>N} is the operator of (24)),
%   Ddag = D^dagger, e = (0,1) of Lemma 11, thetaN = theta_N.

    N = (length(x)-2)/2;

    % only the block |n| <= N_Z1 of DF is inverted, A^{<=N} being the identity for N_Z1 < |n| <= N
    % (x and all the tail bounds use the truncation size N)
    if nargin < 8 || isempty(N_Z1), N_Z1 = N; end
    N_Z1 = min(N_Z1,N);
    ind = N-N_Z1+1 : N+N_Z1+1;   io = [1, 1+ind];

    % approximate inverse A^{<=N} (Section 2.3), computed in floating point (before intval(eps))
    DF = DF_SDDE([x(1); x(1+ind)],eps,N_FFT);
    A = intval(eye(2*N+2));
    A(io,io) = intval(inv(DF));

    eps = intval(eps); nu = intval(nu);
    a = intval(x(2:end));
    n = intval((-N:N)');

    %% Y bound (Lemma 4)

    disp('F_SDDE_rig... ')
    F = F_SDDE_rig(x,eps,C_rho,nu,N_FFT);

    % tail |k| > N of the zeroth row of L_1 G(bar U), bounded with (34)
    F(N+2) = F(N+2) + infsup(-1,1)*2*C_rho*abs(eps(1))*nu/((N+1)*(nu-1)*nu^(N+1));

    AF = A*F ;   % (beta_0, bar W) of Lemma 4, |n| <= N

    Y0 = max([abs(AF(1)) norm(AF(2:end),1) norm(n.*AF(2:end),1)]) ;
    % |n| > N: bar W_n = eps_1 V_n/(in), bounded with (34) (Section 2.6.1)
    Y1 = max(2*C_rho*abs(eps(1))*nu/((N+1)*(nu-1)*nu^(N+1)), 2*C_rho*abs(eps(1))*nu/((nu-1)*nu^(N+1)));
    Y = Y0 + Y1;
    disp(['Y <= ',num2str(sup(Y))])

    %% Z2 bound (Lemma 5)
    disp('Computing Z2... ')
    L1 = intval(zeros(2*N+1));
    for j = [-N:-1 1:N]
        L1(j+N+1,j+N+1) = 1/intval(1i*j) ;
    end
    for j = [-N:-1 1:N]
        L1(N+1,j+N+1) = (-1)^(j+1)/intval(1i*j) ;
    end
    Ddag = [intval(1)  intval(zeros(1,2*N+1)) ;   % D^dagger = diag(1, L_1)
          intval(zeros(2*N+1,1)) L1] ;

    D = diag([1 ; 1i*n]);

    % ||A^{<=N} D^dagger|| and ||D A^{<=N} D^dagger||: first inequality of (60) (Lemma 11)
    thetaN = kolmogorov_theta_N(N);
    if sup(thetaN) > 1, error('theta_N > 1 : Lemma 9 does not apply'), end
    e = intval(zeros(2*N+2,1)); e(N+2) = 1;
    DADdag = D*A*Ddag;
    normAD  = comp_mat_norm(A*Ddag,N_FFT) + thetaN*norm(A*e,1);
    normDAD = comp_mat_norm(DADdag,N_FFT) + thetaN*norm(D*(A*e),1);

    % max{||A D^dagger||, ||D A D^dagger||} (second part of Lemma 11), PiJ = pi^{<=N} J
    PiJ = diag([intval(0) ; intval(ones(2*N+1,1))]);
    norm_ADdag = max(normAD + thetaN, comp_mat_norm(DADdag-PiJ,N_FFT) + thetaN*norm(D*(A*e),1) + 1);
    norm_der2_u = norm(n.^2.*a,1) ;   % >= ||bar u''||_{C^0}

    ae1 = abs(eps(1)) ; ae2 = abs(eps(2));

    Z2 = norm_ADdag*ae1*ae2*(ae2*norm_der2_u + 2);   % Lemma 5, with hat Z_2 of (41)

    disp(['Z2 <= ',num2str(sup(Z2))])

    %% Z1 bound (Lemmas 9 and 10)
    disp('Computing Z1... ')

    norm_der_u = norm(n.*a,1) ;   % >= ||bar u'||_{C^0}

    % ||d_u g(bar u)||_{B(C^0)} and sup_{||h||_{C^1}=1} ||d_t[d_u g(bar u)h]||_{C^0} (Lemma 10)
    norm_dug  = 1 + ae1*(1 + ae2*norm_der_u);
    norm_dtdug = 1 + ae1*(1 + ae2*norm_der_u)*(1 + ae2*norm_der2_u) + ae1*ae2*norm_der_u;

    Z11 = max(normAD,normDAD)*norm_dug*thetaN;
    Z12 = norm_dtdug*thetaN;

    DF =  DF_SDDE_rig(x,eps,C_rho_der,C_rho_col,nu,N_FFT);
    defect = intval(eye(2*N+2)) - A*DF;   % I - A^{<=N} DF(bar x) on the Galerkin space

    % Z_{1,0} of (54), with the second inequality of (60) for B = I - A DF(bar x) and
    % B = D(I - A DF(bar x)): pi^{<=N} B pi^{<=N} D^dagger pi^{<=N} = defect*Ddag
    Z10 = max(comp_mat_norm(defect*Ddag,N_FFT) + (1+thetaN)*norm(defect*e,1), ...
              comp_mat_norm(D*defect*Ddag,N_FFT) + (1+thetaN)*norm(D*(defect*e),1));

    disp(['   Z_{1,0} <= ',num2str(sup(Z10)),' , Z_{1,1} <= ',num2str(sup(Z11)), ...
          ' , Z_{1,2} <= ',num2str(sup(Z12))])

    Z1 = Z10 + Z11 + Z12;
    disp(['Z1 <= ',num2str(sup(Z1))])

end
