A = rand(100);
R = detect_vortices(A, struct('k_peaks',80,'plot',true));
M = evaluate_vortex_to_schur(A, R, 40, 5);   % top-40 peaks vs bottom-5% |S|
disp(M.hit_rate), disp(M.hit_rate_baseline), disp(M.spearman_on_all_ranks)
