% === campaign_demo.m ===
n = 100; M = make_hard_matrices(n, 7);

names = fieldnames(M);
for t = 1:numel(names)
    A = M.(names{t});
    fprintf('\n=== %s ===\n', names{t});
    S = run_vortex_vs_schur_demo(A, names{t}); %#ok<NASGU>
end
