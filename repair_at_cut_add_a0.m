function [Anew, info] = repair_at_cut_add_a0(A, cut, tau)
% Minimalistic repair: tweak a0 at (cut,cut) to target |S|>=tau.
% tau is a small threshold (e.g., 1e-6 relative). Returns Anew and |S| before/after.
k = cut;
A_n = A(1:k-1,1:k-1); v1 = A(1:k-1,k); v2 = A(k,1:k-1).'; a0 = A(k,k);
S  = a0 - v2.' * (A_n \ v1);

if abs(S) < tau
    delta = sign(real(S)) * (tau - abs(S));  % move S to ±tau with minimal |Δa0|
    if isnan(delta) || isinf(delta), delta = tau; end
else
    delta = 0;
end
Anew = A; Anew(k,k) = a0 + delta;

info = struct('cut',k, 'S_before',S, 'S_after', S+delta, ...
              'delta_a0',delta, 'hit', abs(S)<tau);
end
