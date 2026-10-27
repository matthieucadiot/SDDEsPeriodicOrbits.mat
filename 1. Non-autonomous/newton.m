function x = newton(x,eps,N_FFT)
% Newton's method for the truncated map F of (29) (forced toy model (4), non-rigorous).
%   x = (beta,a): initial guess, beta (= u(-pi) at a zero, Lemma 2), a = (a_n)_{|n|<=N} Fourier
%   coefficients of u

tol = 1e-15;

F = F_SDDE(x,eps,N_FFT);
nf = norm(F,1);
it = 0;

while (it<=70) && (nf > tol)
    disp(['||F||_1 = ',num2str(nf)])
    DF = DF_SDDE(x,eps,N_FFT);
    x = x - DF\F_SDDE(x,eps,N_FFT);
    F = F_SDDE(x,eps,N_FFT);
    nf = norm(F,1);
    it = it+1;
end

disp(['After Newton, ||F||_1 = ',num2str(nf)])

end
