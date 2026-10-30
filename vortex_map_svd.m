function [Vmap, idx] = vortex_map_svd(A, r)
    if nargin<2, r = max(5, round(min(size(A))/20)); end
    [U,S,V] = svd(A,'econ');
    As = U(:,end-r+1:end) * S(end-r+1:end,end-r+1:end) * V(:,end-r+1:end)';

    [m,n] = size(A); Vmap = zeros(m-1,n-1);
    for i=1:m-1
        for j=1:n-1
            a = As(i,j);    b = As(i+1,j);
            c = As(i,j+1);  d = As(i+1,j+1);
            Vmap(i,j) = abs(a*d - b*c); % 2x2 minor magnitude
        end
    end
    % top-k vortex centers (optional):
    k = max(10, round(0.01*numel(Vmap)));
    [~,ord] = maxk(Vmap(:), k);
    [ii,jj] = ind2sub(size(Vmap), ord);
    idx = [ii, jj];
end
