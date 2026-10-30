% Compare_norms_1_2_Frobenius_percentiles.m
% Percentiles of condition numbers vs n for 1-, 2-, and Frobenius norms.

rng(1);
nset = 5:5:100;           % sizes to sample
T    = 500;               % draws per size

P1 = zeros(numel(nset),3);   % [median, p90, p99] for cond_1
P2 = zeros(numel(nset),3);   % for cond_2
PF = zeros(numel(nset),3);   % for kappa_F

for i = 1:numel(nset)
    n  = nset(i);
    K1 = zeros(T,1); K2 = zeros(T,1); KF = zeros(T,1);

    for t = 1:T
        A    = randn(n,n);
        K1(t)= cond(A,1);
        K2(t)= cond(A);                     % 2-norm (SVD)
        KF(t)= norm(A,'fro') * norm(inv(A),'fro');  % Frobenius
    end

    % keep only finite values (avoid NaN/Inf in percentiles)
    f1 = K1(isfinite(K1));  if isempty(f1), f1 = NaN; end
    f2 = K2(isfinite(K2));  if isempty(f2), f2 = NaN; end
    fF = KF(isfinite(KF));  if isempty(fF), fF = NaN; end

    % IMPORTANT: keep as ROW vectors (1x3), no transpose
    P1(i,:) = prctile(f1,[50 90 99]);
    P2(i,:) = prctile(f2,[50 90 99]);
    PF(i,:) = prctile(fF,[50 90 99]);
end

figure; clf; hold on; box on
set(gca,'YScale','log','XScale','log')
axis([min(nset) max(nset) 1 1e7])

% Heuristic guides
plot(nset, sqrt(nset), '-', 'Color',[0 .66 0], 'DisplayName','n^{1/2}');
plot(nset, 500*(nset.^3), '-', 'Color',[0 .66 0], 'HandleVisibility','off');

% 2-norm
loglog(nset, P2(:,1), 'k-',  'DisplayName','cond_2 median');
loglog(nset, P2(:,3), 'k--', 'DisplayName','cond_2 p99');

% Frobenius
loglog(nset, PF(:,1), 'm-',  'DisplayName','k_F median');
loglog(nset, PF(:,3), 'm--', 'DisplayName','k_F p99');

% 1-norm (optional to show)
loglog(nset, P1(:,1), 'b-',  'DisplayName','cond_1 median');
loglog(nset, P1(:,3), 'b--', 'DisplayName','cond_1 p99');

xlabel('n'); ylabel('condition number (log scale)');
legend('Location','northwest'); title('Condition number percentiles vs matrix size');
grid on
