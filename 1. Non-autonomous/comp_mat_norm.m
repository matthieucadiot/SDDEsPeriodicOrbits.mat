function f = comp_mat_norm(B,N_FFT)
    % Upper bound for the operator norm on H_0 = R x C^0_per of the (2N+2)x(2N+2) interval
    % matrix B (first row and column: R component), Section 2.6.6: Lemma 13 for the Fourier
    % block, ell^1 norms for the R row and column (|U_n| <= ||u||_{C^0}).

    A = B(2:end,2:end);
    N = size(A); N = (N(1)-1)/2;

    [F, Finv, c] = DFT_matrix(N_FFT,N);

    f = c*norm(Finv(:,1:2*N+1)*A*F(1:2*N+1,:),"inf");   % (69)

    f = max(norm(B(1,:),1), f + norm(B(2:end,1),1) );
end
