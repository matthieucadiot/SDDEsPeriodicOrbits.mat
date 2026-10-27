function [F, Finv, c] = DFT_matrix(N_FFT,N)
% Interval matrices Finv (values on the grid of N_FFT points of the trigonometric polynomials
% with coefficients (a_n)_{|n|<=N}) and F = Finv'/N_FFT, and the constant c of Lemma 13, so that
% (69) reads ||B||_{B(C^0)} <= c*||Finv*B*F||_inf.

if abs(log2(N_FFT)-ceil(log2(N_FFT)))>0 || N>N_FFT/4
    disp('Wrong size for N_FFT'), return
end

    Finv = intval(zeros(N_FFT,2*N+1));
    E = intval(eye(N_FFT));
    N0 = N_FFT/2;
    % values on the grid 2*pi*k/N_FFT, k = -N_FFT/2,...,N_FFT/2-1 (the grid t_k of
    % Section 2.5, up to periodicity)
    for j = -N:N
        ej = E(:,j+N0+1);
        Finv(:,j+N+1) = ifftshift(verifyfft(fftshift(ej),-1))*N_FFT;
    end

    F = Finv'/N_FFT;

    % c >= ||V||_{L^1}/(2*pi*cos(pi*N/N_FFT)) (Lemma 13), with the de la Vallee Poussin kernel V
    % of (68) (degree L), using ||V||_{L^1}/(2*pi) <= P/(P-2L) * (1/P) sum_j |V(s_j)| on P
    % equispaced points s_j (the bound C_1 of Section 2.6.6 with (N,N_FFT) replaced by (L,P))
    L = N_FFT-N-1; P = 2^8*N_FFT; k = (-L:L)';
    Vh = min(1,(L+1-abs(k))/(L+1-N));
    V = verifyfft(fftshift(intval([zeros(P/2-L,1);Vh;zeros(P/2-L-1,1)])),-1)*P;
    c = sup(P/intval(P-2*L)*sum(abs(V))/P/cos(intval('pi')*N/N_FFT));

end
