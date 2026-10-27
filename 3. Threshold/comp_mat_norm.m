function f = comp_mat_norm(B,N_FFT)
% Upper bound for the operator norm on H_0 = R^4 x C^0_per (u = (u_1,u_2), max norm) of the
% (4N+6)x(4N+6) interval matrix B with blocks [beta ; u_1 ; u_2] of sizes 4, 2N+1 and 2N+1:
% maximum over the output blocks of the sum of the norms of the input blocks, with Lemma 13 for
% the four Fourier blocks and ell^1 norms for the scalar rows and the beta columns
% (|U_n| <= ||u||_{C^0}).

nB = size(B,1);
N  = (nB-6)/4;                       % B is (4N+6) square

iw = 5:2*N+5;                        % u_1 block
it = 2*N+6:4*N+6;                    % u_2 block

[F, Finv, c] = DFT_matrix(N_FFT,N);
trick = @(M) c*norm(Finv(:,1:2*N+1)*M*F(1:2*N+1,:),"inf");   % (69)

% B(C^0) norms of the four Fourier blocks
nWW = trick(B(iw,iw));   nWT = trick(B(iw,it));
nTW = trick(B(it,iw));   nTT = trick(B(it,it));

% contributions of the beta columns
betaW = norm(B(iw,1),1)+norm(B(iw,2),1)+norm(B(iw,3),1)+norm(B(iw,4),1);
betaT = norm(B(it,1),1)+norm(B(it,2),1)+norm(B(it,3),1)+norm(B(it,4),1);

f = max([ norm(B(1,:),1) ; norm(B(2,:),1) ; norm(B(3,:),1) ; norm(B(4,:),1) ; ...
          nWW + nWT + betaW ; ...
          nTW + nTT + betaT ]);
end
