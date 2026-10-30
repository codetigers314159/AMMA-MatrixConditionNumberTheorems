function [Anew, info] = repair_at_cut_edge_nudge(A, cut, tau, epsMax)
% Nudge v1 and v2 in the direction that increases S, bounded by epsMax (e.g., 1e-3*||A||_F/n)
k = cut; n = size(A,1);
A_n = A(1:k-1,1:k-1); v1 = A(1:k-1,k); v2 = A(k,1:k-1).'; a0 = A(k,k);
x  = A_n \ v1;                % S = a0 - v2' * x
grad_v2 = -x;                 % dS/dv2 = -x
grad_v1 = -(A_n' \ v2);       % dS/dv1 = -(A_n'\v2)

alpha = epsMax / (norm([grad_v1;grad_v2],2) + eps);
dv2   = alpha * grad_v2;
dv1   = alpha * grad_v1;

Anew = A;
Anew(k,1:k-1) = (v2 + dv2).';
Anew(1:k-1,k) = (v1 + dv1);

S0 = a0 - v2.'       * x;
S1 = a0 - (v2+dv2).' * (A_n \ (v1+dv1));

info = struct('cut',k, 'S_before',S0, 'S_after',S1, ...
              'dv1_norm',norm(dv1), 'dv2_norm',norm(dv2));
end
