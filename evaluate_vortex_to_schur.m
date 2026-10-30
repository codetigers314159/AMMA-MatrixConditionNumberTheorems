function M = evaluate_vortex_to_schur(A, R, k, q_percent)
% Correlate joint-vortex peaks with Schur smallness.
% A: matrix, R: result from detect_vortices, k: top-k peaks to test,
% q_percent: define "small S" threshold (e.g., 5 = bottom 5%).
n = size(A,1);
P = R.peaks_joint(1:min(k,size(R.peaks_joint,1)), :);
% map each [i,j] cell to a diagonal cut c slightly beyond the hotspot
marg = 2; cands = min(max(max(P,[],2)+marg, 2), n);
Svals = arrayfun(@(c) schur_at_cut(A, c), cands);

% baseline: sample same number of random cuts
rng(1); baseC = randi([2 n], numel(cands), 1);
Svals_base = arrayfun(@(c) schur_at_cut(A,c), baseC);

% compute thresholds and metrics
AllCuts = unique([cands; baseC]);
AllS = arrayfun(@(c) schur_at_cut(A,c), AllCuts);
thr = quantile(AllS, q_percent/100);

hit_rate    = mean(Svals <= thr);
hit_rate_b  = mean(Svals_base <= thr);

% rank correlation
[~,ord_joint] = sort(Svals,'ascend');
joint_rank = 1:numel(Svals);
[~,ord_all] = sort(AllS,'ascend');
all_rank_map = containers.Map(num2cell(AllCuts), num2cell(1:numel(AllCuts)));
joint_ranks_in_all = arrayfun(@(c) all_rank_map(c), cands);
rho = corr(joint_ranks_in_all(:), (1:numel(cands))', 'type','Spearman');

M = struct('k',k,'q_percent',q_percent,'hit_rate',hit_rate,'hit_rate_baseline',hit_rate_b,...
           'spearman_on_all_ranks',rho,'Svals_joint',Svals,'Svals_baseline',Svals_base,...
           'S_threshold',thr);
end
