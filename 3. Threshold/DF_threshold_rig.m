function [DF, vbar, DFt] = DF_threshold_rig(x,par,Cv,nu,N_FFT,N2)
% Rigorous enclosure of the Galerkin matrix of DF(bar x) (Section 4.1.2) for the threshold SDDE:
% the computations of DF_threshold.m, with interval arithmetic and the aliasing bound (38) (one
% constant C of (34) per function, from Cv) for every Fourier coefficient computed by FFT.
%   vbar = struct of the Fourier coefficients |n| <= N of the components of bar v_1, bar v_2,
%          bar v_3 of (123), of bar v_4 (124), of v(bar w) and of v'(bar w)
%   DFt  = columns N < |k| <= N_2 (u_1 then u_2) of pi^{<=N} DF(bar x) pi^{(N,N_2]} D^dagger,
%          for Z_{1,1} of (127)
% x = (bar beta, bar w, bar tau), par as in DF_threshold.m, nu = bar nu, N2 = N_2 of Lemma 26.

beta = x(1:4);
omega = intval(beta(2));
u = x(5:end); N = (length(u)-2)/4;
w   = u(1:2*N+1);
tau = u(2*N+2:4*N+2);

b       = par(2);
mu      = par(3);
gamma   = par(4);
m       = mid(par(5));
n_q     = mid(par(6));
q_minus = par(7);
q_plus  = par(8);
Theta_q = par(9);
v_minus = par(10);
v_plus  = par(11);
Theta_v = par(12);

q    = @(y)((q_minus*Theta_q^n_q + q_plus*y.^n_q)./(Theta_q^n_q + y.^n_q));
v    = @(y)((v_minus*Theta_v^m   + v_plus*y.^m  )./(Theta_v^m   + y.^m  ));
dqdy = @(y)( n_q.*y.^(n_q-1).*Theta_q.^n_q.*(q_plus-q_minus)./(Theta_q.^n_q+y.^n_q).^2 );
dvdy = @(y)( m  .*y.^(m  -1).*Theta_v.^m  .*(v_plus-v_minus)./(Theta_v.^m  +y.^m  ).^2 );

if abs(log2(N_FFT)-ceil(log2(N_FFT)))>0 || N>N_FFT/4
    disp('Wrong size for N_FFT'), return
end

n = intval((-N:N)');

% w(t_k), tau(t_k) (inverse FFT, Section 2.5.2)
w_pad   = [intval(zeros(N_FFT/2-N,1)); w;   intval(zeros(N_FFT/2-N-1,1))];
tau_pad = [intval(zeros(N_FFT/2-N,1)); tau; intval(zeros(N_FFT/2-N-1,1))];
w_k   = ifftshift(verifyfft(fftshift(w_pad),  -1))*intval(N_FFT);
tau_k = ifftshift(verifyfft(fftshift(tau_pad),-1))*intval(N_FFT);
t = (2*intval('pi')/N_FFT)*(-N_FFT/2:N_FFT/2-1)';

wp = w.*(1i*n);
y_k  = intval(zeros(N_FFT,1));   % w(t-tau)
ty_k = intval(zeros(N_FFT,1));   % w'(t-tau)
for k = 0:N_FFT-1
    y_k(k+1)  = real(sum( w .*exp(1i*n*(t(k+1)-tau_k(k+1))) ));
    ty_k(k+1) = real(sum( wp.*exp(1i*n*(t(k+1)-tau_k(k+1))) ));
end

v_circ_w  = v(w_k);
vp_circ_w = dvdy(w_k);

% derivatives of calG at (w(t-tau), u, omega) on the grid, cf. (120)-(124)
d_1G_1 = omega*b*exp(-mu*omega*tau_k).*v_circ_w.*dqdy(y_k)./v(y_k) ...
       - omega*b*exp(-mu*omega*tau_k).*dvdy(y_k).*v_circ_w.*q(y_k)./(v(y_k).^2);
d_1G_2 = dvdy(y_k).*v_circ_w./(v(y_k).^2);
d_1G_1_ty = d_1G_1.*ty_k;
d_1G_2_ty = d_1G_2.*ty_k;

d_2G_11 = omega*b*exp(-mu*omega*tau_k).*vp_circ_w.*q(y_k)./v(y_k) - omega*gamma;
d_2G_12 = -mu*omega^2*b*exp(-mu*omega*tau_k).*v_circ_w.*q(y_k)./v(y_k);
d_2G_21 = -vp_circ_w./v(y_k);

d_3G_1  = b*exp(-mu*omega*tau_k).*v_circ_w.*q(y_k).*(1-omega*mu*tau_k)./v(y_k) - gamma*w_k;

% Fourier coefficients |n| <= N with the aliasing bound (38), with one constant C per function;
% 1/(nu^N_FFT-1) is written as nu^(-N_FFT)/(1-nu^(-N_FFT)) to avoid overflow
nu_N = exp(-N_FFT*log(intval(nu))); nu_N = nu_N/(1-nu_N);
alias = @(C) infsup(-1,1)*C*(nu.^n + nu.^(-n))*nu_N;
toF = @(h,C) trunc(fftshift(verifyfft(ifftshift(h))/intval(N_FFT)),N_FFT,N) + alias(C);

% columns of d_u g(bar x) (122), from the functions d_u g(bar x) e_k on the grid, and zeroth rows
% (L_1 U)_0 = sum(c0.*U) and int_{-tau(0)}^0 h = sum(th0.*U), with th0 = vartheta (defined after
% (125)), the coefficients |n| >= N_FFT/2 being bounded with (34) and |c0|, |th0| <= 2/|n|
nF = (-N_FFT/2:N_FFT/2-1)'; tau_at_0 = real(sum(intval(tau)));
c0  = (-1).^(nF+1)./(1i*intval(nF+(nF==0))); c0([1 N_FFT/2+1]) = 0;
th0 = (1 - exp(-1i*intval(nF)*tau_at_0))./(1i*intval(nF+(nF==0))); th0(1) = 0; th0(N_FFT/2+1) = tau_at_0;
alias0 = infsup(-1,1)*(nu.^nF + nu.^(-nF))*nu_N;
tail0 = infsup(-1,1)*4*nu^(1-N_FFT/2)/((N_FFT/2)*(nu-1));
F = @(h,C) fftshift(verifyfft(ifftshift(h))/intval(N_FFT)) + C*alias0;

g11 = intval(zeros(2*N+1)); g12 = g11; g21 = g11; g22 = g11;
r11 = intval(zeros(1,2*N+1)); r12 = r11; r21 = r11; r22 = r11; xi1 = r11;
for k = -N:N
    % C = constants C of (34) of the five functions below
    col = exp(1i*intval(k)*(t - tau_k)); ek = exp(1i*intval(k)*t); nk = nu^abs(k);
    C = [Cv.col1 + nk*Cv.v2_1, nk*Cv.v3_1, Cv.col2 + nk*Cv.v2_2, nk*Cv.v3_2, nk*Cv.vpw];
    G = F(d_1G_1.*col + d_2G_11.*ek, C(1)); g11(:,k+N+1) = trunc(G,N_FFT,N); r11(k+N+1) = sum(c0.*G) + C(1)*tail0;
    G = F((d_2G_12 - d_1G_1_ty).*ek, C(2)); g12(:,k+N+1) = trunc(G,N_FFT,N); r12(k+N+1) = sum(c0.*G) + C(2)*tail0;
    G = F(d_1G_2.*col + d_2G_21.*ek, C(3)); g21(:,k+N+1) = trunc(G,N_FFT,N); r21(k+N+1) = sum(c0.*G) + C(3)*tail0;
    G = F(-d_1G_2_ty.*ek, C(4));            g22(:,k+N+1) = trunc(G,N_FFT,N); r22(k+N+1) = sum(c0.*G) + C(4)*tail0;
    xi1(k+N+1) = sum(th0.*F(vp_circ_w.*ek, C(5))) + C(5)*tail0;
end
rb = sum(c0.*F(d_3G_1,Cv.v4)) + Cv.v4*tail0;         % zeroth row of L1*d_3G_1

% DFt: columns N < |k| <= N_2 (u_1 then u_2), divided by ik (D^dagger)
vw0 = v(real(sum(w.*exp(-1i*n*tau_at_0))));         % v(w(-tau(0)))
kt = [-N2:-N-1 N+1:N2]; nt = numel(kt); iL = [1./(1i*intval(-N:-1)'); 0; 1./(1i*intval(1:N)')];
DFt = intval(zeros(4*N+6,2*nt));
for j = 1:nt
    k = kt(j); col = exp(1i*intval(k)*(t - tau_k)); ek = exp(1i*intval(k)*t); nk = nu^abs(k);
    C = [Cv.col1 + nk*Cv.v2_1, nk*Cv.v3_1, Cv.col2 + nk*Cv.v2_2, nk*Cv.v3_2, nk*Cv.vpw];
    G1 = F(d_1G_1.*col + d_2G_11.*ek, C(1)); G2 = F(d_1G_2.*col + d_2G_21.*ek, C(3));
    G3 = F((d_2G_12 - d_1G_1_ty).*ek, C(2)); G4 = F(-d_1G_2_ty.*ek, C(4));
    x1 = omega*(sum(th0.*F(vp_circ_w.*ek, C(5))) + C(5)*tail0);
    DFt(:,j)    = [x1; 0; G1(N_FFT/2+1); G2(N_FFT/2+1); -iL.*trunc(G1,N_FFT,N); -iL.*trunc(G2,N_FFT,N)];
    DFt(:,nt+j) = [omega*vw0; 0; G3(N_FFT/2+1); G4(N_FFT/2+1); -iL.*trunc(G3,N_FFT,N); -iL.*trunc(G4,N_FFT,N)];
    DFt([N+5 3*N+6],[j nt+j]) = -[sum(c0.*G1) + C(1)*tail0, sum(c0.*G3) + C(2)*tail0; sum(c0.*G2) + C(3)*tail0, sum(c0.*G4) + C(4)*tail0];
    DFt(:,[j nt+j]) = DFt(:,[j nt+j])/(1i*k);
end
xb = sum(th0.*F(v_circ_w,Cv.vw)) + Cv.vw*tail0;      % int_{-tau(0)}^0 v(w)

d_1G_1    = toF(d_1G_1,Cv.v1_1);      d_1G_2    = toF(d_1G_2,Cv.v1_2);
d_2G_11   = toF(d_2G_11,Cv.v2_1);     d_2G_12   = toF(d_2G_12,Cv.d2G12);
d_2G_21   = toF(d_2G_21,Cv.v2_2);
d_3G_1    = toF(d_3G_1,Cv.v4);
v_circ_w  = toF(v_circ_w,Cv.vw);      vp_circ_w = toF(vp_circ_w,Cv.vpw);
d_1G_1_ty = toF(d_1G_1_ty,Cv.v1_1ty); d_1G_2_ty = toF(d_1G_2_ty,Cv.v1_2ty);

% bar v_1, bar v_2, bar v_3 of (123)
v1_1 = d_1G_1; v1_2 = d_1G_2;
v2_1 = d_2G_11; v2_2 = d_2G_21;
v3_1 = -d_1G_1_ty + d_2G_12;  v3_2 = -d_1G_2_ty;

vbar = struct('v1_1',v1_1,'v1_2',v1_2,'v2_1',v2_1,'v2_2',v2_2, ...
              'v3_1',v3_1,'v3_2',v3_2,'v4',d_3G_1,'vw',v_circ_w,'vpw',vp_circ_w);

% d_beta xi
xi_beta = intval(zeros(4));
xi_beta(1,2) = real(xb);
xi_beta(3,2) = d_3G_1(N+1);
xi_beta(4,1) = 1;

% L_1 (Section 2.4)
L1 = intval(zeros(2*N+1));
for j = [-N:-1 1:N]
    L1(j+N+1,j+N+1) = 1/(1i*intval(j));
    L1(N+1,j+N+1)   = (-1)^(j+1)/(1i*intval(j));
end

% d_beta f and d_u f = I - L_1 d_u g of (111), with the zeroth rows of L_1 computed above
g_beta = intval(zeros(2*(2*N+1),4));
g_beta(1:2*N+1,2) = -L1*d_3G_1;  g_beta(N+1,2) = -rb;

f11 = intval(eye(2*N+1)) - L1*g11;   f12 = - L1*g12;
f21 = - L1*g21;                      f22 = intval(eye(2*N+1)) - L1*g22;
f11(N+1,:) = ((-N:N)==0) - r11;      f12(N+1,:) = -r12;
f21(N+1,:) = -r21;                   f22(N+1,:) = ((-N:N)==0) - r22;
df  = [[f11 f12];[f21 f22]];

% d_u xi of (125)
xi_u = intval(zeros(4,2*(2*N+1)));
xi_u(1, 1:2*N+1)     = omega*xi1;
xi_u(1, 2*N+2:end)   = omega*intval(ones(1,2*N+1))*vw0;
xi_u(2, 1:2*N+1)     = intval(ones(1,2*N+1));
xi_u(3, 1:2*N+1)     = g11(N+1,:);   xi_u(3, 2*N+2:end) = g12(N+1,:);
xi_u(4, 1:2*N+1)     = g21(N+1,:);   xi_u(4, 2*N+2:end) = g22(N+1,:);

g_beta(N+1,3) = -1;  g_beta(2*N+1 + N+1,4) = -1;

DF = [[xi_beta xi_u];[g_beta df]];

end

function c = trunc(C,N_FFT,N)
% coefficients |n| <= N of a vector of N_FFT Fourier coefficients
    c = C(N_FFT/2-N+1:N_FFT/2+N+1,1);
end
