function G = cond_growth(A, normType)
% Track cond across principal submatrices A(1:k,1:k).
% normType: '2' (MATLAB cond), 'fro' (Frobenius κ_F = ||A||_F ||A^{-1}||_F)
n = size(A,1); c = zeros(n,1);
for k = 2:n
    Ak = A(1:k,1:k);
    switch normType
        case '2'
            c(k) = cond(Ak);                    % 2-norm
        case 'fro'
            X = Ak \ eye(k);                    % no explicit inv
            c(k) = norm(Ak,'fro') * norm(X,'fro');
        otherwise
            error('normType must be ''2'' or ''fro''.');
    end
end
G.k = (1:n).';  G.k(1)=[];       % start at k=2
G.cond = c(2:end);
end

function plot_growth(G, ttl)
semilogy(G.k, G.cond, '-','LineWidth',1.4);
grid on; xlabel('k (principal submatrix size)'); ylabel('cond'); title(ttl);
end
