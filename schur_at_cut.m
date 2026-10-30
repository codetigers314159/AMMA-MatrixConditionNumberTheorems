function Sabs = schur_at_cut(A, c)
% Compute |S| at diagonal cut c (places hotspot near bottom-right of block).
% A is n×n, c in 2..n
k = c;
A_n = A(1:k-1,1:k-1);
v1  = A(1:k-1,k);
v2  = A(k,1:k-1).';
a0  = A(k,k);
Sabs = abs(a0 - v2.' * (A_n \ v1));
end
