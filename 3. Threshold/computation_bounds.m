function [Y, Z1, Z2] = computation_bounds(x,par,r_star,N_FFT,C_g1,C_g2,C_vw,C_col,C_DF,Cv,nu,N_Z1,N2)
% Bounds Y, Z_1 and Z_2(r_*) of Theorem 3 for the threshold SDDE (100) (Section 4.1):
%   Y (Lemma 4), Z_1 = Z_{1,0} + Z_{1,1} + Z_{1,2} (Lemma 26) and Z_2(r_*) (Lemma 27), the
%   operator norms on H_0 being computed with Lemmas 11 and 13.
%   x = (bar beta, bar u_1, bar u_2) (order N),  par as in DF_threshold.m,  r_star = r_*,
%   C_g1, C_g2, C_vw = constants C of (34) for bar g_1, bar g_2, v(bar w) (C_rho_contour_SDDE.m),
%   C_col, C_DF, Cv = output of C_rho_contour_operator.m (only Cv is used),  nu = bar nu,
%   N_Z1 = size of the block of DF inverted to build A^{<=N} (default N),  N2 = N_2 of Lemma 26.
% Notation: (w,tau) = (bar u_1,bar u_2), omega = bar beta_2, A = matrix of A^{<=N} (in the
%   comments, A = A^{<=N} + pi^{>N} is the operator of (115)), Ddag = D^dagger,
%   ew, et = (0,(1,0)), (0,(0,1)), which replace e of Lemma 11 (proof of Theorem 28),
%   thetaN = theta_N, thN2 = theta_{N_2}, lam = lambda of (126), E = calE, X = calX and
%   normA0 = calA of Lemma 26, T = A^{<=N} DF(bar x) pi^{(N,N_2]} D^dagger.

beta  = x(1:4);
omega = intval(beta(2));
nu    = intval(nu);

u = x(5:end); N = (length(u)-2)/4;
w   = u(1:2*N+1);
tau = u(2*N+2:4*N+2);

% only the block |n| <= N_Z1 of DF is inverted, A^{<=N} being the identity for N_Z1 < |n| <= N
% (x and all the tail bounds use the truncation size N)
if nargin < 12 || isempty(N_Z1), N_Z1 = N; end
N_Z1 = min(N_Z1,N);
ind = N-N_Z1+1 : N+N_Z1+1;
io  = [1:4, 4+ind, 4+2*N+1+ind];

n = intval((-N:N)');
Wc = intval(w);  Tc = intval(tau);

% v, v', v'' of (101), for interval evaluations on boxes
m_=mid(par(5)); nq_=mid(par(6)); qm_=par(7); qp_=par(8); Tq_=par(9);
vm_=par(10); vp_=par(11); Tv_=par(12);
vH   = @(z)( (vm_*Tv_^m_ + vp_*z.^m_)./(Tv_^m_ + z.^m_) );
vdH  = @(z)( m_*Tv_^m_*(vp_-vm_).*z.^(m_-1)./(Tv_^m_+z.^m_).^2 );
vddH = @(z)( m_*Tv_^m_*(vp_-vm_).*((m_-1)*z.^(m_-2).*(Tv_^m_+z.^m_)-2*m_*z.^(2*m_-2))./(Tv_^m_+z.^m_).^3 );
mag  = @(z) max(sup(abs(z(:))));

% approximate inverse A^{<=N} (Section 2.3), and Galerkin matrix of DF(bar x)
x_op = [beta; w(ind); tau(ind)];
A = intval(eye(4*N+6));
A(io,io) = intval(inv(DF_threshold(x_op,mid(par),N_FFT)));
[DF_rig, vbar, DFt] = DF_threshold_rig(x,par,Cv,nu,N_FFT,N2);

%% Y bound (Lemma 4)
F = F_threshold_rig(x,par,C_g1,C_g2,C_vw,nu,N_FFT);
AF = A*F;
AF_b = AF(1:4); AF_w = AF(5:2*N+5); AF_t = AF(2*N+6:4*N+6);

Y = max([ abs(AF_b(1)) ; abs(AF_b(2)) ; abs(AF_b(3)) ; abs(AF_b(4)) ; ...
          norm(AF_w,1) ; norm(AF_t,1) ; norm(n.*AF_w,1) ; norm(n.*AF_t,1) ]);
% |n| > N: bar W_n = -(bar g_i)_n/(in), with (bar g_i)_n bounded with (34)
Y = Y + max([ 2*C_g1*nu/((N+1)*(nu-1)*nu^(N+1)) ; 2*C_g1*nu/((nu-1)*nu^(N+1)) ; ...
              2*C_g2*nu/((N+1)*(nu-1)*nu^(N+1)) ; 2*C_g2*nu/((nu-1)*nu^(N+1)) ]);
Y = abs(Y);
disp(['Y  <= ',num2str(sup(Y))])

%% operator norms of Lemma 11
% L_1, Ddag = D^dagger = diag(I_4, L_1, L_1) and D
L1blk = intval(zeros(2*N+1));
for j = [-N:-1 1:N]
    L1blk(j+N+1,j+N+1) = 1/(1i*intval(j));
    L1blk(N+1,j+N+1)   = (-1)^(j+1)/(1i*intval(j));
end
Ddag = intval(zeros(4*N+6));
Ddag(1:4,1:4)               = intval(eye(4));
Ddag(5:2*N+5,5:2*N+5)       = L1blk;
Ddag(2*N+6:4*N+6,2*N+6:4*N+6)= L1blk;
D = diag([intval(ones(4,1)); 1i*n; 1i*n]);

thetaN = kolmogorov_theta_N(N);
if sup(thetaN) > 1, error('theta_N > 1 : Lemma 9 does not apply'), end

% ||A^{<=N} D^dagger|| and ||D A^{<=N} D^dagger||: first inequality of (60), with ew and et
ew = intval(zeros(4*N+6,1)); ew(N+5) = 1;
et = intval(zeros(4*N+6,1)); et(3*N+6) = 1;
DADdag = D*A*Ddag;
normAD  = comp_mat_norm(A*Ddag,N_FFT) + thetaN*(norm(A*ew,1) + norm(A*et,1));
normDAD = comp_mat_norm(DADdag,N_FFT) + thetaN*(norm(D*(A*ew),1) + norm(D*(A*et),1));

% max{||A D^dagger||, ||D A D^dagger||} (second part of Lemma 11), PiJ = pi^{<=N} J
PiJ = diag([intval(zeros(4,1)); intval(ones(4*N+2,1))]);
norm_ADdag = max(normAD + thetaN, comp_mat_norm(DADdag-PiJ,N_FFT) ...
                 + thetaN*(norm(D*(A*ew),1) + norm(D*(A*et),1)) + 1);

% upper bounds of norms of bar u
normWp  = norm(n.*Wc,1);     % ||bar u_1'||_{C^0}
normWpp = norm(n.^2.*Wc,1);  % ||bar u_1''||_{C^0}
normTp  = norm(n.*Tc,1);     % ||bar u_2'||_{C^0}
normT   = norm(Tc,1);        % ||bar u_2||_{C^0}

% sum_{|n|>N} |h_n| and sum_{|n|>N} |n h_n| when |h_n| <= C/nu^|n|
tailC  = @(C) 2*C*nu^(-N-1)/(1-1/nu);
dtailC = @(C) 2*C*nu^(-N-1)*(N+1-N/nu)/(1-1/nu)^2;

%% Z1 bound (Lemma 26)
% Z_{1,0} of (54), with the second inequality of (60) for B = I - A DF(bar x) and
% B = D(I - A DF(bar x)), e being replaced by ew and et
defect = intval(eye(4*N+6)) - A*DF_rig;
Z10 = max( comp_mat_norm(defect*Ddag,N_FFT) + (1+thetaN)*(norm(defect*ew,1) + norm(defect*et,1)) , ...
           comp_mat_norm(D*defect*Ddag,N_FFT) + (1+thetaN)*(norm(D*(defect*ew),1) + norm(D*(defect*et),1)) );

% upper bounds of ||(bar v_j)_i||_{C^0} and ||(bar v_j)_i'||_{C^0}, using (34) for |n| > N
nv1_1=norm(vbar.v1_1,1)+tailC(Cv.v1_1);  nd1_1=norm(n.*vbar.v1_1,1)+dtailC(Cv.v1_1);
nv1_2=norm(vbar.v1_2,1)+tailC(Cv.v1_2);  nd1_2=norm(n.*vbar.v1_2,1)+dtailC(Cv.v1_2);
nv2_1=norm(vbar.v2_1,1)+tailC(Cv.v2_1);  nd2_1=norm(n.*vbar.v2_1,1)+dtailC(Cv.v2_1);
nv2_2=norm(vbar.v2_2,1)+tailC(Cv.v2_2);  nd2_2=norm(n.*vbar.v2_2,1)+dtailC(Cv.v2_2);
nv3_1=norm(vbar.v3_1,1)+tailC(Cv.v3_1);  nd3_1=norm(n.*vbar.v3_1,1)+dtailC(Cv.v3_1);
nv3_2=norm(vbar.v3_2,1)+tailC(Cv.v3_2);  nd3_2=norm(n.*vbar.v3_2,1)+dtailC(Cv.v3_2);

nv1 = max(nv1_1,nv1_2); nv2 = max(nv2_1,nv2_2); nv3 = max(nv3_1,nv3_2);
nd1 = max(nd1_1,nd1_2); nd2 = max(nd2_1,nd2_2); nd3 = max(nd3_1,nd3_2);

% range of bar u_1 (K boxes) and bar u_1(-bar u_2(0)), for calX
W0 = real(Wc(N+1)); sWv = sum(abs(Wc)) - abs(Wc(N+1));
T0 = real(Tc(N+1)); sTv = sum(abs(Tc)) - abs(Tc(N+1));
K = 400;
Yb0 = split(W0 + infsup(-1,1)*sWv,K);
tau_at_0 = real(sum(Tc));
w_mtau0  = sum(Wc.*exp(-1i*n*tau_at_0));

% Z_{1,1} and Z_{1,2} of (127)
T    = A*DFt;
thN2 = kolmogorov_theta_N(N2);
lam  = Cv.lam;
if lam*(N2+1) <= N, error('lam*(N2+1) <= N : increase N2'), end
E = 4*nu^(N+1)/((nu-1)*(N2+1))*( max(Cv.v1_1,Cv.v1_2)*nu^(-intval(lam)*(N2+1))/(1-nu^(-intval(lam))) ...
                                + (max(Cv.v2_1,Cv.v2_2) + max(Cv.v3_1,Cv.v3_2))*nu^(-N2)/(nu-1) );
X = abs(omega)*( (2*mag(vdH(Yb0)) + abs(tau_at_0)*mag(vddH(Yb0))*normWp)*thN2^2 + mag(vH(w_mtau0))*thN2 );
normA0 = 0;    % calA: columns of A and D*A for beta_1, beta_3, beta_4 and the constants c_1, c_2
for i = [1 3 4 N+5 3*N+6]
    normA0 = max([normA0 sup(norm(A(:,i),1)) sup(norm(D*A(:,i),1))]);
end
Z11 = max(comp_tail_norm(T,N_FFT,N2), comp_tail_norm(D*T,N_FFT,N2)) ...
      + normA0*(X + 2*E + 2*thetaN*thN2*(nv1 + nv2 + nv3)) + max(normAD,normDAD)*E;
Z12 = tailC(Cv.v4) + ((1+normTp)*nv1 + nv2 + nv3 + nd1 + nd2 + nd3)*thetaN;

Z1 = abs(Z10 + Z11 + Z12);
disp(['Z1 <= ',num2str(sup(Z1))])

%% Z2 bound (Lemma 27, with r = r_*)
% J_r of (131), with y and z_1 in K sub-boxes
rr = sup(r_star);
Yb  = split(W0 + infsup(-1,1)*(sWv + rr),K);   % y
Z1b = Yb;                                       % z_1
Z2b = T0 + infsup(-1,1)*(sTv + rr);            % z_2
Wb  = omega + infsup(-1,1)*rr;                  % omega

GD = Gderivs_threshold(Yb*ones(1,K),ones(K,1)*Z1b.',Z2b,Wb,par);   % all the pairs (y,z_1)

p1 = 1 + normWp;

% calM_k(r_*) of (132), k = 1,2,3 (calG_2 does not depend on z_2 and omega)
M1 = p1*max(mag(GD.d1d1G1),mag(GD.d1d1G2)) + max(mag(GD.dz1d1G1),mag(GD.dz1d1G2)) ...
        + mag(GD.dz2d1G1) + mag(GD.dwdyG1);
M2 = p1*max([mag(GD.dz1d1G1) mag(GD.dz2d1G1) mag(GD.dz1d1G2)]) ...
        + max([mag(GD.dz1dz1G1) mag(GD.dz1dz2G1) mag(GD.dz1dz1G2)]) ...
        + max([mag(GD.dz1dz2G1) mag(GD.dz2dz2G1)]) ...
        + max([mag(GD.dwdz1G1) mag(GD.dwdz2G1)]);
M3 = p1*mag(GD.dwdyG1) + mag(GD.dwdz1G1) + mag(GD.dwdz2G1) + mag(GD.dwdwG1);

% Z_{21},...,Z_{25} of Lemma 27
Z21 = mag(vH(Yb)) + normT*mag(vdH(Yb));
Z22 = (normT + rr)*mag(vdH(Yb)) + abs(omega)*( mag(vdH(Yb)) + normT*mag(vddH(Yb)) );
Z23 = mag(vH(Yb)) + abs(omega)*mag(vdH(Yb))*p1;
Z24 = M3;
Z25 = M1*(p1 + rr) + 2*M2 + nv1*(2 + normWpp);

Z2 = abs(norm_ADdag*max(Z21 + Z22 + Z23, Z24 + Z25));   % (133)
disp(['Z2(r_*) <= ',num2str(sup(Z2))])

end

function f = comp_tail_norm(T,N_FFT,N2)
% Upper bound for the operator norm on H_0 of the (4N+6) x 4(N2-N) interval matrix T, whose
% columns are the modes N < |k| <= N_2 of u_1 and u_2: as in comp_mat_norm.m, with Lemma 13 for
% N_2 instead of N (as described after Lemma 26).

N  = (size(T,1)-6)/4;
iw = 5:2*N+5;  it = 2*N+6:4*N+6;
nt = size(T,2)/2;  i1 = 1:nt;  i2 = nt+1:2*nt;

[~, Finv] = DFT_matrix(N_FFT,N);
[F, ~, c] = DFT_matrix(N_FFT,N2);
F = F([1:N2-N N2+N+2:2*N2+1],:);      % modes N < |k| <= N2
trick = @(M) c*norm(Finv*M*F,"inf");

f = max([ norm(T(1,:),1) ; norm(T(2,:),1) ; norm(T(3,:),1) ; norm(T(4,:),1) ; ...
          trick(T(iw,i1)) + trick(T(iw,i2)) ; ...
          trick(T(it,i1)) + trick(T(it,i2)) ]);
end

function S = split(I,K)
% K sub-intervals covering the interval I
e = linspace(inf(I),sup(I),K+1)'; e([1 end]) = [inf(I) sup(I)];
S = infsup(e(1:K),e(2:K+1));
end
