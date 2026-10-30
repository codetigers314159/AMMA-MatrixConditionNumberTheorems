function S = run_vortex_vs_schur_demo(A, tag)
% Runs growth charts and applies three strategies, producing bar plots.

if nargin<2, tag = 'Demo'; end
n = size(A,1);

% Baseline growth
G2  = cond_growth(A,'2');
GFr = cond_growth(A,'fro');

figure('Color','w','Position',[60 60 1160 460]);
subplot(1,2,1); plot_growth(G2,  sprintf('%s — cond_2 growth',tag));
subplot(1,2,2); plot_growth(GFr, sprintf('%s — cond_F growth',tag));

% Global cond before
cond2_before  = cond(A);
condF_before  = norm(A,'fro') * norm(A\eye(n),'fro');

% --- Strategy A: Schur-first (scan best |S| among proposed cuts from joint peaks)
R  = detect_vortices(A, struct('k_peaks', max(40, round(0.01*(n-1)^2)), 'plot', true, 'title',[tag ' — Vortex maps']));
C  = propose_block_cuts(R.peaks_joint, n, 'margin',2);
T  = rank_cuts_by_schur(A, C);              % sorted by |S|
cuts_A = T.cut(1:min(5,height(T)));         % try up to 5 top cuts

A_S = A;
S_logs = [];
for t = 1:numel(cuts_A)
    [A_S, info] = repair_at_cut_edge_nudge(A_S, cuts_A(t), 1e-8, 1e-3*norm(A,'fro')/n);
    S_logs = [S_logs; info]; %#ok<AGROW>
end

% --- Strategy B: Vortex-only (top joint peaks → cuts → same repair)
cuts_B = min(max(max(R.peaks_joint,[],2)+2,2), n);
cuts_B = unique(cuts_B,'stable'); cuts_B = cuts_B(1:min(5,numel(cuts_B)));
A_V = A;
for t = 1:numel(cuts_B)
    [A_V, ~] = repair_at_cut_edge_nudge(A_V, cuts_B(t), 1e-8, 1e-3*norm(A,'fro')/n);
end

% --- Strategy C: 2×2 minor shrink (greedy local)
[A_M, ~] = repair_by_2x2_minors(A, 20, 0.1);

% After cond
cond2_after_S = cond(A_S); condF_after_S = norm(A_S,'fro')*norm(A_S\eye(n),'fro');
cond2_after_V = cond(A_V); condF_after_V = norm(A_V,'fro')*norm(A_V\eye(n),'fro');
cond2_after_M = cond(A_M); condF_after_M = norm(A_M,'fro')*norm(A_M\eye(n),'fro');

figure('Color','w','Position',[80 80 980 380]);
subplot(1,2,1);
bar([cond2_before cond2_after_S cond2_after_V cond2_after_M]);
set(gca,'YScale','log','XTickLabel',{'before','Schur-guided','Vortex-guided','2x2-greedy'});
ylabel('cond_2'); title([tag ' — 2-norm cond before/after']); grid on;
subplot(1,2,2);
bar([condF_before condF_after_S condF_after_V condF_after_M]);
set(gca,'YScale','log','XTickLabel',{'before','Schur-guided','Vortex-guided','2x2-greedy'});
ylabel('cond_F'); title([tag ' — Frobenius cond before/after']); grid on;

S = struct('cond2_before',cond2_before,'condF_before',condF_before, ...
           'cond2_after_S',cond2_after_S,'condF_after_S',condF_after_S, ...
           'cond2_after_V',cond2_after_V,'condF_after_V',condF_after_V, ...
           'cond2_after_M',cond2_after_M,'condF_after_M',condF_after_M, ...
           'schur_logs',S_logs);
end
