function M = make_hard_matrices(n, seed)
% Returns a struct of high-cond test matrices at size n.
if nargin<2, seed = 1; end, rng(seed);

S = struct();
S.randsvd_hi  = gallery('randsvd', n, 1e12, 3);       % controlled κ2 via decaying σ
S.kahan       = gallery('kahan',   n);                % classic ill-conditioned
S.parter      = gallery('parter',  n);                % near-Toeplitz, big κ
S.lotkin      = gallery('lotkin',  n);                % brutal κ growth
S.minij       = gallery('minij',   n);                
S.lehmer      = gallery('lehmer',  n);
S.dr_rand     = diag(logspace(0,-10,n)) * randn(n);   % diagonal-dominated, near singular
end
