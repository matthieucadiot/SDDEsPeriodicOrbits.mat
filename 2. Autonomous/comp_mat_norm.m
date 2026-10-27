function f = comp_mat_norm(B,N_FFT)
    % Upper bound for the operator norm on H_0 = R^2 x C^0_per of the (2N+3)x(2N+3) interval
    % matrix B (rows and columns 1,2: R^2 component), Section 2.6.6: Lemma 13 for the Fourier
    % block, ell^1 norms for the scalar rows and columns (|U_n| <= ||u||_{C^0}).

    A = B(3:end,3:end);
    N = size(A,1);  N = (N-1)/2;

    [F, Finv, c] = DFT_matrix(N_FFT,N);

    normA = c*norm(Finv(:,1:2*N+1)*A*F(1:2*N+1,:),"inf");   % (69)

    % maximum over the output blocks of the sum of the norms of the input blocks
    f = max([ norm(B(1,:),1) , ...
              norm(B(2,:),1) , ...
              normA + norm(B(3:end,1),1) + norm(B(3:end,2),1) ]);
end
