function C = propose_block_cuts(peaks, n, varargin)
% PROPOSE_BLOCK_CUTS  Convert vortex cell anchors -> (i_cut=j_cut) candidates.
% We use a single diagonal cut (square A_n) so Schur is scalar a0 at (c,c).
% peaks: kx2 [row,col] from any map (e.g., R.peaks_joint)
% n    : matrix size
% opts.margin (default 2): push the cut a couple steps beyond the hotspot

p = inputParser; p.addParameter('margin',2); p.parse(varargin{:});
marg = p.Results.margin;

k=size(peaks,1);
C = repmat(struct('cut',[],'i_cut',[],'j_cut',[],'weight',1.0),k,1);
for t=1:k
    i=peaks(t,1); j=peaks(t,2);
    c = min(max(max(i,j)+marg, 2), n);  % keep in 2..n
    C(t).cut   = c;
    C(t).i_cut = c; C(t).j_cut = c;
end
end
