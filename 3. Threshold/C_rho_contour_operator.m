function [C_col, C_DF, Cv] = C_rho_contour_operator(x,par,bar_nu,N_FFT,N2)
% Constants C of (34) for DF(bar x) and the bound Z_1 (Lemma 26), computed as in (35)-(36)
% (maximum over the lines Im(z) = -+ bar rho, bar rho = log(bar nu)), and lambda of (126):
%   C_col : C_col(k+N2+1) for the column e^{ik(t-bar tau(t))} of L_3(bar u_2), |k| <= N_2
%   C_DF  : maximum of max(C_col) and of the constants v1_1,...,v1_2ty of Cv (not used)
%   Cv    : struct with the constants for the components of bar v_1, bar v_2, bar v_3 of (123)
%           (v1_1,...,v3_2), bar v_4 = d_omega calG_1 (v4), v(bar w) (vw), v'(bar w) (vpw),
%           d_{z_2} calG_1 (d2G12) and (bar v_1)_i bar w'(t-bar tau) (v1_1ty, v1_2ty); for
%           (bar v_1)_i e^{ik(t-bar tau(t))}, i = 1,2, maximum over |k| <= N_2 (col1, col2);
%           max(C_col) (col); and lambda of (126) (lam)
% x = (bar beta, bar w, bar tau), par as in DF_threshold.m, N2 = N_2 of Lemma 26.

beta  = x(1:4);
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

q   = @(y)((q_minus*Theta_q^n_q + q_plus*y.^n_q)./(Theta_q^n_q + y.^n_q));
v   = @(y)((v_minus*Theta_v^m   + v_plus*y.^m  )./(Theta_v^m   + y.^m  ));
qd  = @(y)( n_q.*y.^(n_q-1).*Theta_q.^n_q.*(q_plus-q_minus)./(Theta_q.^n_q+y.^n_q).^2 );
vd  = @(y)( m  .*y.^(m  -1).*Theta_v.^m  .*(v_plus-v_minus)./(Theta_v.^m  +y.^m  ).^2 );

bar_rho = sup(log(intval(bar_nu)));   % bar_rho >= log(bar_nu)
n  = (-N:N)';
jc = (-N2:N2);                     % column indices |k| <= N_2

ipi = intval('pi'); delta = infsup(0,1)*2*ipi/N_FFT;

% w(z) and tau(z) for z in [t_{k-1},t_k] -+ i bar_rho, t_k = 2*pi*k/N_FFT, as in (35)
% (delta = [0,2*pi/N_FFT] instead of [-pi/N_FFT,pi/N_FFT])
aw_p=[zeros(N_FFT/2-N,1); w.*exp(n*( bar_rho+1i*delta)); zeros(N_FFT/2-N-1,1)];
aw_m=[zeros(N_FFT/2-N,1); w.*exp(n*(-bar_rho+1i*delta)); zeros(N_FFT/2-N-1,1)];
at_p=[zeros(N_FFT/2-N,1); tau.*exp(n*( bar_rho+1i*delta)); zeros(N_FFT/2-N-1,1)];
at_m=[zeros(N_FFT/2-N,1); tau.*exp(n*(-bar_rho+1i*delta)); zeros(N_FFT/2-N-1,1)];
Uw{1}=verifyfft(ifftshift(aw_p),-1)*N_FFT; Uw{2}=verifyfft(ifftshift(aw_m),-1)*N_FFT;
Ut{1}=verifyfft(ifftshift(at_p),-1)*N_FFT; Ut{2}=verifyfft(ifftshift(at_m),-1)*N_FFT;
sgn = [ bar_rho, -bar_rho ];

names = {'v1_1','v1_2','v2_1','v2_2','v3_1','v3_2','v4','vw','vpw','d2G12','v1_1ty','v1_2ty'};
for s = 1:2, acc{s} = intval(zeros(numel(names),1)); Gcol{s} = intval(zeros(1,2*N2+1)); G1col{s} = Gcol{s}; G2col{s} = Gcol{s}; end

for s = 1:2
    br = sgn(s);
    for k = 1:N_FFT
        wz = Uw{s}(k);  tz = Ut{s}(k);
        ph = exp(n*(br+1i*delta - 1i*tz)).*exp(1i*n*(k-1)*2*ipi/N_FFT);   % e^{in(z-tau(z))}
        yv  = sum( w.*ph );                       % w(z-tau(z))
        ypv = sum( (1i*n).*w.*ph );               % w'(z-tau(z))
        vY=v(yv); vpY=vd(yv); qY=q(yv); qpY=qd(yv);
        vZ=v(wz); vpZ=vd(wz);
        E0 = exp(-mu*omega*tz);

        % derivatives of calG at (w(z-tau(z)), u(z), omega), cf. (120)-(124)
        d1G1 = omega*b*E0*vZ.*( qpY./vY - vpY.*qY./vY.^2 );
        d1G2 = vpY.*vZ./vY.^2;
        d2G11= omega*b*E0*vpZ./vY.*qY - omega*gamma;
        d2G21= -vpZ./vY;
        d2G12= -mu*omega^2*b*E0*vZ./vY.*qY;
        d3G1 = b*E0*vZ./vY.*qY.*(1-omega*mu*tz) - gamma*wz;
        v3_1 = -d1G1.*ypv + d2G12;
        v3_2 = -d1G2.*ypv;

        vals = [ d1G1; d1G2; d2G11; d2G21; v3_1; v3_2; d3G1; vZ; vpZ; ...
                 d2G12; d1G1.*ypv; d1G2.*ypv ];
        acc{s} = acc{s} + sup(abs(vals));

        col = exp(jc*(br+1i*delta - 1i*tz)).*exp(1i*jc*(k-1)*2*ipi/N_FFT);  % e^{ik(z-tau(z))}
        Gcol{s} = Gcol{s} + sup(abs(col));
        G1col{s} = G1col{s} + sup(abs(d1G1*col));   % d_1G_1(z) e^{ik(z-tau(z))}
        G2col{s} = G2col{s} + sup(abs(d1G2*col));
    end
end

C1 = sup(acc{1})/N_FFT;  C2 = sup(acc{2})/N_FFT;  Cvec = max(C1,C2);
C_col = max( sup(Gcol{1})/N_FFT , sup(Gcol{2})/N_FFT );

Cv.v1_1=Cvec(1); Cv.v1_2=Cvec(2); Cv.v2_1=Cvec(3);  Cv.v2_2=Cvec(4);
Cv.v3_1=Cvec(5); Cv.v3_2=Cvec(6); Cv.v4=Cvec(7);    Cv.vw=Cvec(8);
Cv.vpw=Cvec(9);  Cv.d2G12=Cvec(10);
Cv.v1_1ty=Cvec(11); Cv.v1_2ty=Cvec(12);
Cv.col = max(C_col(:));
Cv.col1 = max([sup(G1col{1}) sup(G1col{2})])/N_FFT;   % columns of d_1G_1 L_3(bar u_2)
Cv.col2 = max([sup(G2col{1}) sup(G2col{2})])/N_FFT;   % columns of d_1G_2 L_3(bar u_2)
% lambda of (126): Im(tau(t+i bar_rho)) <= (1-lam)*bar_rho, using both lines
Cv.lam = inf(1 - max([sup(-imag(Ut{1})); sup(imag(Ut{2}))])/intval(bar_rho));

C_DF = max([ Cvec ; Cv.col ]);

end
