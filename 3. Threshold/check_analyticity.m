function check = check_analyticity(x,par,bar_nu,N_FFT)
% Verification (Section 4.1.1) that bar g_1, bar g_2 and bar v_1,...,bar v_4 are analytic on
% S_{bar rho}, bar rho = log(bar nu): by Lemma 23, it suffices that the images of the lines
% Im(z) = -+ bar rho by bar w and by bar w(.-bar tau(.)) satisfy the half-plane condition (118)
% (avoidpole.m) for the roots (117) of theta_v^m + y^m, and of v^- theta_v^m + v^+ y^m and
% theta_q^n + y^n, respectively. The roots of theta_v^m + y^m are also excluded for
% bar w(.-bar tau(.)), where v and v' are evaluated separately in the code.
% x = (bar beta, bar w, bar tau), par as in DF_threshold.m.

u = x(5:end); N = (length(u)-2)/4;
w   = u(1:2*N+1);
tau = u(2*N+2:4*N+2);

m = mid(par(5)); n_q = mid(par(6)); Theta_q = par(9);
v_minus = par(10); v_plus = par(11); Theta_v = par(12);

% roots (117)
k = (0:m-1)';
p_v  = Theta_v*exp(1i*intval('pi')*(2*k+1)/m);                              % Theta_v^m + y^m
p_1v = Theta_v*exp(log(intval(v_minus)/v_plus)/m)*exp(1i*intval('pi')*(2*k+1)/m); % v_minus*Theta_v^m + v_plus*y^m
k = (0:n_q-1)';
p_q  = Theta_q*exp(1i*intval('pi')*(2*k+1)/n_q);                            % Theta_q^n + y^n

bar_rho = sup(log(intval(bar_nu)));   % bar_rho >= log(bar_nu)
n = (-N:N)';
ipi = intval('pi'); delta = infsup(0,1)*2*ipi/N_FFT;

% enclosures of bar w and bar w(.-bar tau(.)) on the cells [t_k,t_{k+1}] of both lines, as in (35)
W1 = intval(complex(zeros(N_FFT,1))); W2 = W1; Wd1 = W1; Wd2 = W1;

for k = 1:N_FFT
    e1 = exp(n*( bar_rho+1i*delta)).*exp(1i*n*k*2*ipi/N_FFT);   % e^{inz} on Im z = -bar_rho
    e2 = exp(n*(-bar_rho+1i*delta)).*exp(1i*n*k*2*ipi/N_FFT);   % e^{inz} on Im z = +bar_rho
    W1(k)  = sum(w.*e1);                                       % w(z)
    W2(k)  = sum(w.*e2);
    Wd1(k) = sum(w.*e1.*exp(-1i*n*sum(tau.*e1)));              % w(z-tau(z))
    Wd2(k) = sum(w.*e2.*exp(-1i*n*sum(tau.*e2)));
end

W = [W1;W2]; Wd = [Wd1;Wd2];

check = true;
for j = 1:m
    check = check && avoidpole(p_v(j),W) && avoidpole(p_v(j),Wd) && avoidpole(p_1v(j),Wd);
end
for j = 1:n_q
    check = check && avoidpole(p_q(j),Wd);
end

end
