function S = make_hard_matrices(n, seed)
% MAKE_HARD_MATRICES  High-cond test matrices of size n.

n=100;

if nargin < 2, seed = 1; end
rng(seed);

S = struct();

% Safer gallery calls (wrapped in try/catch in case a name isn't available)
S.randsvd_hi = gallery('randsvd', n, 1e12, 3);   % κ2 ≈ 1e12
S.kahan      = gallery('kahan',   n);
try,  S.parter = gallery('parter', n);  catch, S.parter = []; end
try,  S.lotkin = gallery('lotkin', n);  catch, S.lotkin = []; end
S.minij      = gallery('minij',   n);
S.lehmer     = gallery('lehmer',  n);

p='fro'

fprintf( 'randsvd_hi 1e12, 3')
cond(S.randsvd_hi,p)

fprintf( 'kahan')
cond(S.kahan,p )

fprintf( 'parter')
cond(S.parter,p)

fprintf( 'lotkin')
cond(S.lotkin,p)

fprintf( 'minij')
cond(S.minij,p  )

fprintf( 'lehmer')
cond(S.lehmer,p )



% Diagonal dominance collapse — near singular
S.dr_rand    = diag(logspace(0,-10,n)) * randn(n);

fprintf( 'dr_rand Diagonal dominance collapse — near singular')
cond(S.lehmer,p )

% Optional extras (uncomment if you want more pain)
% S.hilbert    = hilb(n);                  % famously awful κ
% S.wilkinson  = gallery('wilkinson', n);  % clustered eigenvalues
end
