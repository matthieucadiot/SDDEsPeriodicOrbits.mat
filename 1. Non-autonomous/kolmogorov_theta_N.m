function [ThN,d] = kolmogorov_theta_N(N,opts,d0)
%KOLMOGOROV_THETA_N   Verified enclosure (INTLAB interval) of the constant theta_N of (46),
%
%       theta_N = (1/pi) int_0^pi |kappa_N(t)| dt,  kappa_N(t) = pi - t - 2 sum_{k=1}^N sin(kt)/k,
%
%   the optimal constant in ||pi^{>N} h||_{C^0} <= theta_N ||h'||_{C^0} (Lemma 6), computed with
%   Proposition 8:
%     Step 1: Newton's method gives approximations taubar_j of the zeros tau_j of kappa_N in
%             (t_{j-1},t_j), t_j = 2*pi*j/(2N+1);
%     Step 2: a Newton-Kantorovich test encloses tau_j in [taubar_j-r_j, taubar_j+r_j];
%     Step 3: theta_N is enclosed with (53).
%   Only Steps 2 and 3 enter the proof.
%
%   [ThN,d] = KOLMOGOROV_THETA_N(N,opts,d0) also returns the offsets d = taubar - m of the zeros
%   from the midpoints m_j = (2j-1)*pi/(2N+1); the output d at N-1 can be passed as d0 (optional
%   starting point for Newton's method).
%   opts (optional struct, default values in brackets):
%     eps      radius of the a priori box in the Z bound      [1e-9]
%     tol      Newton stopping tolerance                      [1e-11]
%     maxit    maximum number of Newton steps                 [8]
%     blkelem  block size (matrix entries) of the trig sums   [2e6]
%   Cost: O(N^2) in time and O(opts.blkelem) in memory.

if nargin < 2 || isempty(opts), opts = struct; end
if nargin < 3, d0 = []; end
def = struct('eps',1e-9,'tol',1e-11,'maxit',8,'blkelem',2e6);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opts,f{i}), opts.(f{i}) = def.(f{i}); end
end
if ~(isscalar(N) && N == round(N) && N >= 1)
    error('kolmogorov_theta_N:N','N must be a positive integer');
end

Pi = 4*atan(intval(1));             % verified enclosure of pi

% grid and constants; the indices k are in descending order, which reduces the rounding errors
% of the length-N trigonometric sums (rad(theta_N) then grows like N instead of N^2)
j    = (1:N)';
kk   = N:-1:1;   jr = kk(:);
h    = 2*pi/(2*N+1);
m    = (2*j-1)*pi/(2*N+1);          % midpoints of (t_{j-1},t_j)
tlo  = (j-1)*h;   thi = j*h;        % t_{j-1} , t_j
blk  = max(1,floor(opts.blkelem/N));
q    = 51 - max(0,ceil(log2(N)));   % see QUANT: makes taubar*k exact

% Step 1: starting point and Newton's method (non-rigorous)
tb = m;
if ~isempty(d0)
    n0 = min(numel(d0),N);
    tb(1:n0) = tb(1:n0) + d0(1:n0)*((2*N-1)/(2*N+1));
end
tb = quant(tb,q);

w1 = 1./jr;   cold = inf;
for it = 1:opts.maxit
    kap  = pi - tb - 2*tsum(tb,kk,w1,blk,@sin);
    dkap = -sin(((2*N+1)/2)*tb)./sin(tb/2);
    corr = kap./dkap;
    tb   = quant(min(max(tb-corr,tlo+0.05*h),thi-0.05*h),q);
    cnew = max(abs(corr));
    if cnew <= opts.tol || cnew > 0.5*cold, break, end   % converged / stalled
    cold = cnew;
end
d = tb - m;

% Step 2: verified enclosure of the zeros tau_j
% check that the products taubar_j*k, k = 1,...,N, are exact in floating point (see QUANT)
if max(tb)*N*pow2(1,q) > 2^53
    error('kolmogorov_theta_N:exactness', ...
          'products taubar*k are not exact at N = %d',N);
end
W1 = intval(1)./intval(jr);         % 1/k   (descending order)
W2 = intval(1)./intval(jr.^2);      % 1/k^2 (jr.^2 exact in double)

kapI = Pi - intval(tb) - 2*tsum(tb,kk,W1,blk,@sin);     % kappa_N(taubar)
dkI  = -Dn(N,intval(tb));                               % kappa_N'(taubar)
Ieps = midrad(tb,opts.eps*ones(N,1));
dkE  = -Dn(N,Ieps);                                     % kappa_N'([taubar+-eps])

if any(mig(dkI) == 0)
    error('kolmogorov_theta_N:derivative', ...
          'kappa_N''(taubar) not bounded away from 0 at N = %d',N);
end
Y = mag(kapI./dkI);                      % rigorous upper bound of |kappa_N/kappa_N'|
Z = mag(1 - dkE./dkI);
if ~all(Z < 1)
    error('kolmogorov_theta_N:Zbound','Z >= 1 at N = %d',N);
end
r = sup(intval(Y)./(1-intval(Z)));       % rigorous upper bound of Y/(1-Z)
if ~( all(r <= opts.eps) && all(tb-r > tlo) && all(tb+r < thi) )
    error('kolmogorov_theta_N:verification', ...
          'Newton-Kantorovich test failed at N = %d (max r = %.3e)',N,max(r));
end
% tau_j in [tb_j-r_j, tb_j+r_j], contained in (t_{j-1},t_j); by Proposition 8, these are all
% the zeros of kappa_N in (0,pi)

% Step 3: enclosure of theta_N with (53)
% Kcal_N(tau_j) = Kcal_N(tb) + kappa_N(tb)*(tau-tb) + kappa_N'(xi)*(tau-tb)^2/2
Kbar = Pi.*intval(tb) - intval(tb).^2/2 + 2*tsum(tb,kk,W2,blk,@cos);
Kj   = Kbar + kapI.*infsup(-r,r) + 0.5*dkE.*infsup(0*r,sup(intval(r).^2));

% alternating sum, shifted by c ~ Kcal_N(tau_j) to avoid cancellations
c   = Pi^2/3;
sg  = (-1).^(j-1);
S1  = c*mod(N,2) + sum(sg.*(Kj-c));
S2  = sum(W2);                      % sum_{k=1}^N 1/k^2
S3  = sum(((-1).^jr).*W2);          % sum_{k=1}^N (-1)^k/k^2

ThN = (2/Pi)*S1 - (2/Pi)*S2 + (-1)^N*( Pi/2 + (2/Pi)*S3 );
end
% =====================================================================
function y = Dn(N,t)
%DN   Dirichlet kernel  D_N(t) = sin((N+1/2)t)/sin(t/2)  (= -kappa_N').
M = intval(2*N+1)/2;                    % (2N+1)/2 exact, division by 2 exact
y = sin(M.*t)./sin(t/2);
end
% ---------------------------------------------------------------------
function s = tsum(tb,kk,w,blk,fun)
%TSUM   s(i) = sum_k fun(tb(i)*k)*w(k), computed blockwise.
%   tb is quantised by QUANT, so that every product tb(i)*k is exact. kk is the row vector of
%   the indices k and w the column of weights: double (non-rigorous) or intval (rigorous).
%   kk is passed separately since NUMEL is not overloaded for intval.
n = numel(tb);
rig = isa(w,'intval');
if rig, s = intval(zeros(n,1)); else, s = zeros(n,1); end
for i0 = 1:blk:n
    i1 = min(i0+blk-1,n);
    T  = tb(i0:i1)*kk;                  % exact outer product
    if rig, T = intval(T); end
    s(i0:i1) = fun(T)*w;
end
end
% ---------------------------------------------------------------------
function y = quant(x,q)
%QUANT   Round to a multiple of 2^-q (exact operations).
%   With q = 51 - ceil(log2(N)) and 0 < x < pi, x*k is exactly representable in double
%   precision for every k = 1,...,N.
s = pow2(1,q);
y = round(x*s)/s;
end
