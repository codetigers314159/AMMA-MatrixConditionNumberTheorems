A = rand(100);
saveDir = 'Z:\MATRIX VORTICES SVD'
fprintf('DETECT VORTICES RUNNING - ON MATRIX RAND A(100)'); 

R = detect_vortices(A, struct('k_peaks',60,'plot',true,'saveDir',saveDir))

%'Z:\APP_COND_NUMBER_6_2\TechInfo\vortex_demo'));

C = propose_block_cuts(R.peaks_joint, size(A,1), 'margin',2);

T = rank_cuts_by_schur(A, C);    % look at T(1:10,:)

fprintf('plots saved in: %',saveDir); 