function T = rank_cuts_by_schur(A, C)
% RANK_CUTS_BY_SCHUR  Evaluate |S| = |a0 - v2'*(A_n \ v1)| for proposed cuts.
% C: struct array from propose_block_cuts (fields i_cut, j_cut).
% Returns table sorted by |S| ascending (smaller = more precarious).

n = size(A,1);
m = numel(C);
cut  = zeros(m,1); Sval = zeros(m,1); ok = true(m,1);

for t=1:m
    ic=C(t).i_cut; jc=C(t).j_cut; cut(t)=ic;
    if ic<2 || jc<2 || ic>n || jc>n || ic~=jc
        ok(t)=false; Sval(t)=NaN; continue;
    end
    k = ic;   % square top-left
    A_n = A(1:k-1, 1:k-1);
    v1  = A(1:k-1, k);
    v2  = A(k, 1:k-1).';
    a0  = A(k,k);
    % Solve A_n x = v1 stably (no explicit inv)
    try
        x = A_n \ v1;
        S = a0 - (v2.' * x);
        Sval(t) = abs(S);
    catch
        ok(t)=false; Sval(t)=Inf;
    end
end

T = table(cut, Sval, ok, 'VariableNames',{'cut','absS','ok'});
T = sortrows(T, 'absS', 'ascend');
end
