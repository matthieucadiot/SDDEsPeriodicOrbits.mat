function s = Gderivs_threshold(Yb,Z1b,Z2b,Wb,par)
% Interval enclosures, for (y,z_1,z_2,omega) in the box (Yb,Z1b,Z2b,Wb), of the first and second
% derivatives of calG(y,z,omega) of (108),
%   calG_1 = omega b e^{-mu omega z_2} v(z_1)/v(y) q(y) - omega gamma z_1,
%   calG_2 = 1 - v(z_1)/v(y),
% for the bounds calM_k of (132) (Lemma 27). In the field names, d1 = d_y, dz1 = d_{z_1},
% dz2 = d_{z_2}, d3 = dw = d_omega, and G1, G2 are the components of calG.

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

% Hill functions (101) and their first two derivatives
v   = @(z)( (v_minus*Theta_v^m + v_plus*z.^m)./(Theta_v^m + z.^m) );
vd  = @(z)( m*Theta_v^m*(v_plus-v_minus).*z.^(m-1)./(Theta_v^m+z.^m).^2 );
vdd = @(z)( m*Theta_v^m*(v_plus-v_minus).*((m-1)*z.^(m-2).*(Theta_v^m+z.^m) ...
              - 2*m*z.^(2*m-2))./(Theta_v^m+z.^m).^3 );
q   = @(z)( (q_minus*Theta_q^n_q + q_plus*z.^n_q)./(Theta_q^n_q + z.^n_q) );
qd  = @(z)( n_q*Theta_q^n_q*(q_plus-q_minus).*z.^(n_q-1)./(Theta_q^n_q+z.^n_q).^2 );
qdd = @(z)( n_q*Theta_q^n_q*(q_plus-q_minus).*((n_q-1)*z.^(n_q-1-1).*(Theta_q^n_q+z.^n_q) ...
              - 2*n_q*z.^(2*n_q-2))./(Theta_q^n_q+z.^n_q).^3 );

y = Yb; z1 = Z1b; z2 = Z2b; om = Wb;
E  = exp(-mu*om*z2);

vY=v(y); vpY=vd(y); vppY=vdd(y);
qY=q(y); qpY=qd(y); qppY=qdd(y);
vZ=v(z1); vpZ=vd(z1); vppZ=vdd(z1);

% first derivatives
s.d1G1  = om*b*E*vZ.*( qpY./vY - vpY.*qY./vY.^2 );
s.d1G2  = vpY.*vZ./vY.^2;
s.dz1G1 = om*b*E*vpZ./vY.*qY - om*gamma;
s.dz1G2 = -vpZ./vY;
s.dz2G1 = -mu*om.^2*b*E.*vZ./vY.*qY;
s.d3G1  = b*E.*vZ./vY.*qY.*(1-om*mu.*z2) - gamma*z1;

% second derivatives in (y,z_1,z_2)
s.d1d1G1 = om*b*E*vZ.*( qppY./vY - 2*qpY.*vpY./vY.^2 - vppY.*qY./vY.^2 + 2*vpY.^2.*qY./vY.^3 );
s.d1d1G2 = vZ.*( vppY./vY.^2 - 2*vpY.^2./vY.^3 );
s.dz1d1G1 = om*b*E*vpZ.*( qpY./vY - vpY.*qY./vY.^2 );
s.dz1d1G2 = vpY.*vpZ./vY.^2;
s.dz2d1G1 = -mu*om.^2*b*E.*vZ.*( qpY./vY - vpY.*qY./vY.^2 );
s.dz2d1G2 = intval(0);
s.dz1dz1G1 = om*b*E*vppZ./vY.*qY;
s.dz1dz1G2 = -vppZ./vY;
s.dz1dz2G1 = -mu*om.^2*b*E.*vpZ./vY.*qY;
s.dz2dz2G1 = mu.^2*om.^3*b*E.*vZ./vY.*qY;

% second derivatives involving omega
s.dwdyG1  = b*E.*vZ.*( qpY./vY - vpY.*qY./vY.^2 ).*(1-om*mu.*z2);
s.dwdz1G1 = b*E*vpZ./vY.*qY.*(1-om*mu.*z2) - gamma;
s.dwdz2G1 = b*E.*vZ./vY.*qY.*( mu.^2*z2.*om.^2 - 2*mu*om );
s.dwdwG1  = b*E.*vZ./vY.*qY.*( mu.^2*z2.^2.*om - 2*mu*z2 );

end
