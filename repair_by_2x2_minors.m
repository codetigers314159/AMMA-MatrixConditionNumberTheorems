function [Anew, edits] = repair_by_2x2_minors(A, steps, eta)
% Greedy: repeatedly shrink the largest 2x2 minor by soft rescaling.
% steps: how many cells to fix; eta: damping (e.g., 0.1)
Anew = A; edits = [];
for t=1:steps
    M = minor_map(Anew);                   % from detect_vortices helpers
    [~,ord] = maxk(M(:),1);
    [i,j] = ind2sub(size(M), ord);
    % 2x2 block
    a = Anew(i,j); b = Anew(i+1,j); c = Anew(i,j+1); d = Anew(i+1,j+1);
    % target: reduce |ad - bc|
    det2 = a*d - b*c;
    if det2==0, continue; end
    % small balanced correction on (a,d) vs (b,c)
    scale = 1 - eta * sign(det2);    % pull towards bc/ad balance
    Anew(i,j)     = a * scale;
    Anew(i+1,j+1) = d * scale;
    Anew(i+1,j)   = b / scale;
    Anew(i,j+1)   = c / scale;
    edits = [edits; i,j,det2]; %#ok<AGROW>
end
end
