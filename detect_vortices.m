function R = detect_vortices(A, opts)
% DETECT_VORTICES  Flag "vortices" (local swirl/instability) from A- and B-views.
%   R = DETECT_VORTICES(A)                  % defaults, makes a 2-panel figure
%   R = DETECT_VORTICES(A, struct(...))     % override any defaults
%
% Outputs: V_A, V_inv, Z_A, Z_inv, Joint, peaks_* (A,inv,joint), scores_*, meta

% DOCUMENTATION - 10-16-25

% PROPOSE_BLOCK_CUTS
% ------------------
% Converts vortex coordinates into suggested (row,col) split points.
% The goal is to make the vortex fall near the bottom-right of the
% submatrix inspected by Theorem 6.2, where the Schur complement S lives.
%
% Each "cut" defines A_n, v1, v2, a0 so you can compute S = a0 - v2'*inv(A_n)*v1.
% Small |S| means the local block is nearly singular — a prime target for
% repair or condition-number reduction.


% ---------- sanity ----------
[n,m] = size(A); assert(n==m,'A must be square.');
if nargin<2 || ~isstruct(opts), opts=struct(); end

% ---------- defaults ----------
r_small   = def(opts,'r_small',   max(5, round(n/20)));   % # smallest-σ modes
p_rhs     = def(opts,'p_rhs',     8);                     % RHS for inverse-view
k_peaks   = def(opts,'k_peaks',   max(10, round(0.01*(n-1)^2)));
seed      = def(opts,'seed',      []);
normalize = def(opts,'normalize', true);
doplot    = def(opts,'plot',      true);
saveDir   = def(opts,'saveDir',   '');
titleStr  = def(opts,'title',     'Vortex Maps');

% ---------- preproc ----------
if normalize, nf=norm(A,'fro'); if nf>0, A=A/nf; end, end
if ~isempty(seed), rng(seed); end

% ---------- A-view: smallest-σ reconstruction (no inv) ----------
As = smallest_sigma_recon(A, r_small);
V_A = minor_map(abs(As));                 % 2x2 minor magnitude map

% ---------- B-view: solve-based inverse proxy (no inv formed) ----------
[L,U,P] = lu(A,'vector');
Rhs = sign(randn(n,p_rhs));               % Rademacher RHS
X   = U \ (L \ Rhs(P,:));                 % columns ≈ inv(A) columns
g   = mean(abs(X),2);                     % per-row inverse leverage
G   = g * g.';                            % simple energy proxy
V_inv = minor_map(G);

% ---------- normalize & fuse (joint agreement) ----------
Z_A   = zmap(V_A);
Z_inv = zmap(V_inv);
Joint = Z_A .* Z_inv; Joint(Joint<0)=0;   % keep agreement only

% ---------- peaks ----------
[peaks_A,   scores_A]   = topk_peaks(V_A,   k_peaks);
[peaks_inv, scores_inv] = topk_peaks(V_inv, k_peaks);
[peaks_joint, scores_joint] = topk_peaks(Joint, k_peaks);

% ---------- kappa proxy (best-effort) ----------
[svals_headtail,kappa_est] = kappa_proxy(A);

% ---------- plot ----------
if doplot
    figure('Color','w','Position',[100 100 1300 560]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile; imagesc(V_A); axis image; colorbar;
    title(sprintf('%s — A-view (small-σ modes)',titleStr),'Interpreter','none');
    hold on; plot_peaks(peaks_A,'k.'); plot_peaks(peaks_joint,'wo'); hold off;

    nexttile; imagesc(V_inv); axis image; colorbar;
    title(sprintf('%s — B-view (solve-based)',titleStr),'Interpreter','none');
    hold on; plot_peaks(peaks_inv,'k.'); plot_peaks(peaks_joint,'wo'); hold off;

    if ~isempty(saveDir)
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        exportgraphics(gcf, fullfile(saveDir,'vortex_maps.png'),'Resolution',220);
    end
end

% ---------- save CSVs ----------
if ~isempty(saveDir)
    writematrix(V_A,           fullfile(saveDir,'V_A.csv'));
    writematrix(V_inv,         fullfile(saveDir,'V_inv.csv'));
    writematrix(Z_A,           fullfile(saveDir,'Z_A.csv'));
    writematrix(Z_inv,         fullfile(saveDir,'Z_inv.csv'));
    writematrix(Joint,         fullfile(saveDir,'Joint.csv'));
    writematrix(peaks_A,       fullfile(saveDir,'peaks_A.csv'));
    writematrix(peaks_inv,     fullfile(saveDir,'peaks_inv.csv'));
    writematrix(peaks_joint,   fullfile(saveDir,'peaks_joint.csv'));
end

% ---------- pack ----------
R = struct();
R.V_A=V_A; R.V_inv=V_inv; R.Z_A=Z_A; R.Z_inv=Z_inv; R.Joint=Joint;
R.peaks_A=peaks_A; R.peaks_inv=peaks_inv; R.peaks_joint=peaks_joint;
R.scores_A=scores_A; R.scores_inv=scores_inv; R.scores_joint=scores_joint;
R.meta=struct('n',n,'r_small',r_small,'p_rhs',p_rhs,'k_peaks',k_peaks,...
              'normalized',normalize,'kappa_est',kappa_est,'svals_headtail',svals_headtail);
end

% ===== helpers =====
function x = def(s, f, d), if isfield(s,f) && ~isempty(s.(f)), x=s.(f); else, x=d; end, end
function Z = zmap(M), mu=mean(M(:)); sg=std(M(:)); Z=(M-mu)/(sg+eps); end
function As = smallest_sigma_recon(A, r_small)
    n = size(A,1);
    if n<=400
        [U,S,V]=svd(A,'econ'); r=min(r_small,size(S,1));
        As = U(:,end-r+1:end)*S(end-r+1:end,end-r+1:end)*V(:,end-r+1:end)';
    else
        r = min(r_small,max(6,round(n/50)));
        try
            [Vsm,Ssm,Usm]=svds(A,r,'smallest'); As=Usm*Ssm*Vsm';
        catch
            Omega=randn(n,2*r); Y=A*(A'*Omega); [Q,~]=qr(Y,0);
            B=Q'*A; [Ub,Sb,Vb]=svd(B,'econ'); r=min(r,size(Sb,1));
            As=(Q*Ub(:,end-r+1:end))*Sb(end-r+1:end,end-r+1:end)*Vb(:,end-r+1:end)';
        end
    end
end
function M = minor_map(X)  % 2x2 minor magnitudes |ad-bc|
    [m,n]=size(X); M=zeros(m-1,n-1);
    for i=1:m-1
        a=X(i,1:n-1); b=X(i+1,1:n-1); c=X(i,2:n); d=X(i+1,2:n);
        M(i,:)=abs(a.*d - b.*c);
    end
end
function [idx,sc]=topk_peaks(M,k)
    v=M(:); [sc,ord]=maxk(v, min(k,numel(v))); [ii,jj]=ind2sub(size(M),ord); idx=[ii,jj];
end
function plot_peaks(P,style), if isempty(P),return; end, plot(P(:,2),P(:,1),style,'MarkerSize',7,'LineWidth',1.1); end
function [ht,kappa]=kappa_proxy(A)
    n=size(A,1); ht=[]; kappa=NaN;
    try
        if n<=300, s=svd(A); kappa=s(1)/s(end); ht=[s(1:min(5,end)); NaN; s(max(1,end-4):end)];
        else, smax=svds(A,1); smin=svds(A,1,'smallest'); kappa=smax/smin; ht=[smax;NaN;smin];
        end
    catch, end
end
