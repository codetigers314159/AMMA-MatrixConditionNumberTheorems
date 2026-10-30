function [Vinv, idx] = vortex_map_inverse(A, p)
    if nargin<2, p = 8; end
    [m,n] = size(A); assert(m==n);
    R = sign(randn(n,p));  % Rademacher samples
    [L,U,P] = lu(A,'vector');
    X = U \ (L \ R(P,:)); % solves approximate columns of inv(A)
    g = mean(abs(X),2);    % inverse "leverage" per row
    G = g * g.';          % simple energy outer product

    Vinv = zeros(n-1,n-1);
    for i=1:n-1
        for j=1:n-1
            a=G(i,j); b=G(i+1,j); c=G(i,j+1); d=G(i+1,j+1);
            Vinv(i,j) = abs(a*d - b*c);
        end
    end
    k = max(10, round(0.01*numel(Vinv)));
    [~,ord] = maxk(Vinv(:), k);
    [ii,jj] = ind2sub(size(Vinv), ord);
    idx = [ii, jj];
end
