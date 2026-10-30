function out = fp_pipeline_tools(cmd, varargin)
% FP_PIPELINE_TOOLS  Floating-point hygiene helpers for linear algebra.
%
%   out = fp_pipeline_tools('stable_solve', A, b, opts)
%   out = fp_pipeline_tools('scale_equilibrate', A, opts)
%   out = fp_pipeline_tools('residual_stats', A, x, b)
%   out = fp_pipeline_tools('kahan_sum', x)
%   out = fp_pipeline_tools('sigma_min_tail', As)      % As: cell array or 3-D stack
%   out = fp_pipeline_tools('detect_subnormals', x)
%
% Options (struct, name-value also OK):
%   opts.scale        : true/false (default: true)  % row/col equilibration
%   opts.max_scale_it : integer iterations (default: 2)
%   opts.structure    : 'auto'|'spd'|'unsym' (default: 'auto')
%   opts.rank_tol     : tolerance for rank (default: max(size(A))*eps(norm(A,2)))
%   opts.refine       : true/false iterative refinement (default: false)
%   opts.refine_steps : iterations (default: 2)
%
% Returns:
%   For 'stable_solve': struct with fields x, stats (kappa2, kappaF, rcond1, digits, sminInv, scaled)
%   For 'scale_equilibrate': struct with fields Ahat, Dr, Dc
%   For 'residual_stats': struct with residual norms (2, inf), Kahan sums, and ulp deltas
%   For 'kahan_sum': scalar sum with compensation
%   For 'sigma_min_tail': struct with inv_smin vector and percentiles
%   For 'detect_subnormals': struct with mask, count, fraction
%
% Author: (you)

% ---- parse
if nargin==0, help fp_pipeline_tools; return; end
if ~ischar(cmd) && ~isstring(cmd), error('First arg must be command string.'); end
cmd = lower(string(cmd));
opts = struct;
if ~isempty(varargin)
    last = varargin{end};
    if isstruct(last)
        opts = last;
        varargin(end) = [];
    else
        % allow name-value pairs
        if mod(numel(varargin),2)==0
            opts = struct(varargin{:});
            varargin = {};
        end
    end
end

switch cmd
    case "stable_solve"
        [A,b] = varargin{:};
        out = stable_solve(A,b,opts);

    case "scale_equilibrate"
        A = varargin{1};
        out = scale_equilibrate(A,opts);

    case "residual_stats"
        [A,x,b] = varargin{:};
        out = residual_stats(A,x,b);

    case "kahan_sum"
        x = varargin{1};
        out = kahan_sum(x);

    case "sigma_min_tail"
        As = varargin{1};
        out = sigma_min_tail(As);

    case "detect_subnormals"
        x = varargin{1};
        out = detect_subnormals(x);

    otherwise
        error('Unknown command: %s', cmd);
end
end

% ======================== IMPLEMENTATIONS ===============================

function out = stable_solve(A,b,opts)
arguments
    A double
    b double
    opts.scale (1,1) logical = true
    opts.max_scale_it (1,1) double = 2
    opts.structure (1,1) string = "auto"   % 'auto'|'spd'|'unsym'
    opts.rank_tol (1,1) double = []
    opts.refine (1,1) logical = false
    opts.refine_steps (1,1) double = 2
end

[m,n] = size(A);
if m~=n, error('A must be square.'); end

scaled = false; Dr = []; Dc = [];
Ahat = A; bhat = b;

if opts.scale
    S = scale_equilibrate(A, struct('max_scale_it',opts.max_scale_it));
    Ahat = S.Ahat; Dr = S.Dr; Dc = S.Dc;
    if ~isempty(b)
        bhat = Dr \ b;  % because Ahat = Dr \ (A * Dc)
    end
    scaled = true;
end

% decide structure
struc = opts.structure;
if struc=="auto"
    if issymmetric(Ahat) && all(diag(Ahat)~=0)
        % cheap SPD probe: try chol with 'lower' and catch
        [~,p] = chol((Ahat+Ahat')*0.5);
        struc = (p==0) * "spd" + (p~=0) * "unsym";
    else
        struc = "unsym";
    end
end

% rank tolerance
if isempty(opts.rank_tol)
    opts.rank_tol = max(size(Ahat)) * eps(norm(Ahat,2));
end

% solve by preference: chol (SPD) → QR with pivoting → SVD
xhat = [];
method = "";

try
    if struc=="spd"
        R = chol(Ahat);                           % Ahat = R' * R
        xhat = R \ (R' \ bhat);
        method = "chol";
    end
catch
    % fall back
end

if isempty(xhat)
    % QR with column pivoting
    [Q,R,E] = qr(Ahat,0);
    r = rank(R, opts.rank_tol);
    if r == n && all(abs(diag(R))>0)
        y = Q' * bhat;
        xhat = E * (R \ y);
        method = "qr";
    else
        % rank-deficient → minimum-norm via SVD
        [U,S,V] = svd(Ahat,'econ');
        s = diag(S);
        tol = max(size(Ahat)) * eps(max(s));
        sInv = 1./s;  sInv(s < tol) = 0;          % Tikhonov-like cutoff
        xhat = V * ( (sInv .* (U' * bhat)) );
        method = "svd";
    end
end

% iterative refinement (optional)
if opts.refine
    for k = 1:opts.refine_steps
        r = bhat - Ahat*xhat;
        dx = linsolve_stable(Ahat, r, struc);
        xhat = xhat + dx;
    end
end

% unscale x if we scaled columns
x = xhat;
if scaled && ~isempty(Dc)
    x = Dc \ xhat;
end

% stats
k2 = cond(Ahat);
s  = svd(Ahat);
kF = hypot(norm(s),0) * hypot(norm(1./s),0);
r1 = rcond(Ahat);
digits = max(0, 16 - log10(k2));
sminInv = 1/s(end);

out.x = x;
out.stats = struct('method',method,'kappa2',k2,'kappaF',kF, ...
    'rcond1',r1,'digits',digits,'sminInv',sminInv, ...
    'scaled',scaled,'Dr',Dr,'Dc',Dc);
end

function dx = linsolve_stable(A, r, struc)
switch struc
    case "spd"
        R = chol(A);
        dx = R \ (R' \ r);
    otherwise
        [Q,R,E] = qr(A,0);
        dx = E * (R \ (Q' * r));
end
end

function out = scale_equilibrate(A,opts)
arguments
    A double
    opts.max_scale_it (1,1) double = 2
end
% Simple iterative max-abs row/col equilibration:
% Ahat = Dr \ (A * Dc), with Dr,Dc diagonal, aiming for O(1) row/col norms.

n = size(A,1);
Dr = speye(n); Dc = speye(n);
Ahat = A;

for it = 1:opts.max_scale_it
    % Row scaling
    rn = max(abs(Ahat), [], 2);
    rn(rn==0) = 1;
    Sr = spdiags(1./rn, 0, n, n);
    Ahat = Sr * Ahat;  Dr = Dr * (Sr \ speye(n));  % maintain Ahat = Dr^{-1} A Dc

    % Column scaling
    cn = max(abs(Ahat), [], 1);
    cn(cn==0) = 1;
    Sc = spdiags(1./cn(:), 0, n, n);
    Ahat = Ahat * Sc;  Dc = Dc * Sc;
end

out = struct('Ahat',Ahat,'Dr',Dr,'Dc',Dc);
end

function out = residual_stats(A,x,b)
r = b - A*x;

% naive sums
r2  = norm(r,2);
rInf= norm(r,inf);

% compensated 2-norm squared (via Kahan on squares)
rsq = kahan_sum(r.^2);
r2_kahan = sqrt(rsq);

% order sensitivity demo for dot product A*x (first column only, illustrative)
if ~isempty(A)
    y1 = A(:,1).' * x(1);
    y2 = sum(A(:,1) .* x(1)); % equivalent here, placeholder
    ulp = ulpdist(y1,y2);
else
    ulp = NaN;
end

out = struct('r2',r2,'r2_kahan',r2_kahan,'rInf',rInf,'ulp_demo',ulp);
end

function s = kahan_sum(x)
c = 0; s = 0;
for i=1:numel(x)
    y = x(i) - c;
    t = s + y;
    c = (t - s) - y;
    s = t;
end
end

function out = sigma_min_tail(As)
% Accepts cell array of matrices OR 3-D stack A(:,:,k)
inv_smin = [];

if iscell(As)
    for k=1:numel(As)
        s = svd(As{k});
        inv_smin(end+1,1) = 1/s(end); %#ok<AGROW>
    end
elseif ndims(As)==3
    K = size(As,3);
    inv_smin = zeros(K,1);
    for k=1:K
        s = svd(As(:,:,k));
        inv_smin(k) = 1/s(end);
    end
else
    % single matrix
    s = svd(As);
    inv_smin = 1/s(end);
end

prcs = [50 90 99];
pct  = prctile(inv_smin, prcs);

out = struct('inv_smin',inv_smin,'percentiles',table(prcs',pct','VariableNames',{'p','value'}));
end

function out = detect_subnormals(x)
mask = (abs(x) < realmin) & (x ~= 0);
cnt  = nnz(mask);
frac = cnt / numel(x);
out  = struct('mask',mask,'count',cnt,'fraction',frac);
end

function d = ulpdist(x,y)
xi = typecast(x,'uint64'); yi = typecast(y,'uint64');
d  = double(abs(double(xi) - double(yi)));
end
