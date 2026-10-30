function report = analyzeMatrixBasics(M, outDir)
% analyzeMatrixBasics  Quick, publication-ready matrix diagnostics.
% Inputs:
%   M      : numeric matrix (real or complex, square or rectangular)
%   outDir : folder to save figures and summary.txt (default 'figs/matrix_report')
% Outputs:
%   report : struct with key stats
%
% Figures saved:
%   values_hist.png         - histogram of entries (real; or magnitudes if complex)
%   counts_bar.png          - bars: pos/neg/zero, integer/single-like/double-like, complex count
%   eig_scatter.png         - eigenvalues Re vs Im (square case)
%   edge_block_metrics.png  - if min(size)>=2: ||A_n||_F, ||V1||, ||V2||, |b0|
%
% Notes:
% - “Symmetric” means symmetric if real, Hermitian if complex.
% - Integer test is by value (fractional part ≈ 0), not by storage type.
% - “single-like” means value rounds to single precision within relative tol.

    if nargin < 2 || isempty(outDir), outDir = 'figs/matrix_report'; end
    if ~exist(outDir,'dir'), mkdir(outDir); end

    % ---- tolerances ----
    tol_sym     = 1e-9;     % relative Frobenius tolerance for symmetry
    tol_zero    = 1e-12;    % zero test on value/magnitude
    tol_integer = 1e-12;    % integer-by-value test
    rtol_single = 1e-7;     % single-like relative rounding tolerance

    [m,n]   = size(M);
    isSq    = (m == n);
    issp    = issparse(M);
    isReal  = isreal(M);

    % ---- symmetry (Hermitian if complex) ----
    symRel = NaN;
    if isSq
        if isReal
            symRel = norm(M - M.', 'fro') / max(1, norm(M,'fro'));
        else
            symRel = norm(M - M',  'fro') / max(1, norm(M,'fro')); % Hermitian check
        end
    end
    isSym = isSq && (symRel <= tol_sym);

    % ---- trace / diagonal stats ----
    dn = min(m,n);
    d  = diag(M(1:dn,1:dn));
    tr = sum(d);
    if isempty(d)
        diagStats = struct('mean',NaN,'std',NaN,'min',NaN,'max',NaN);
    else
        diagStats = struct('mean',mean(d), 'std',std(d), 'min',min(d), 'max',max(d));
    end

    % ---- values & basic extrema ----
    vals = M(:);
    % plotting measure: real values if real, magnitude if complex
    if isReal
        valsPlot = real(vals);
    else
        valsPlot = abs(vals);
    end

    % min/max by real part (for a simple scalar view)
    realVals = real(vals);
    if isempty(realVals)
        minVal = NaN; maxVal = NaN; rangeVal = NaN;
        min3 = []; max3 = [];
    else
        minVal = min(realVals);
        maxVal = max(realVals);
        rangeVal = maxVal - minVal;
        [sorted,~] = sort(realVals);
        k3 = min(3,numel(sorted));
        min3 = sorted(1:k3);
        max3 = sorted(end-k3+1:end);
    end

    % ---- sign counts (by real part) & zero count (by magnitude) ----
    posCnt  = sum(realVals >  0);
    negCnt  = sum(realVals <  0);
    zeroCnt = sum(abs(vals)   <= tol_zero);  % true zero by magnitude

    % ---- integer / single-like / double-like counts (by value) ----
    % integer test: |x - round(x)| <= max(tol, rel tol)
    relscale = max(1, abs(vals));
    isInt    = abs(vals - round(vals)) <= max(tol_integer, 1e-12 .* relscale);
    isInt    = reshape(isInt, size(vals));
    intCnt   = sum(isInt);

    % single-like: rounds to single within relative tolerance
    isSingleLike = ~isInt & ( abs(vals - single(vals)) <= rtol_single .* relscale );
    singleCnt    = sum(isSingleLike);
    doubleCnt    = numel(vals) - intCnt - singleCnt;
    if doubleCnt < 0, doubleCnt = 0; end  % guard

    % ---- complex count ----
    complexCnt = sum(abs(imag(vals)) > 0);

    % ---- eigenvalues & determinant (square only) ----
    eigvals    = [];
    detVal     = NaN;
    logabsdet  = NaN;
    detSign    = NaN;
    if isSq
        % eigenvalues (may be complex)
        eigvals = eig(full(M));  % full for sparse safety

        % determinant via LU: det(P)*prod(diag(U)); accumulate in log to avoid under/overflow
        try
            [L,U,P] = lu(full(M)); %#ok<ASGLU>
            diagU   = diag(U);
            logabsdet = sum(log(abs(diagU)));
            % sign(det) for real matrices; for complex, this is not meaningful in the same way
            if isReal
                detSign = det(P) * prod(sign(diagU));
                detVal  = detSign * exp(logabsdet);
            else
                detSign = NaN;
                detVal  = exp(logabsdet); % magnitude only
            end
        catch
            % LU may fail for NaNs/Inf; leave det fields as NaN
        end
    end

    % ---- Edge block (Natural Matrix Partition) if size >= 2 ----
    edge = struct();
    if m >= 2 && n >= 2
        n0 = min(m,n) - 1;
        An = M(1:n0,1:n0);
        V1 = M(1:n0,n0+1);
        V2 = M(n0+1,1:n0);
        b0 = M(n0+1,n0+1);
        edge.AnF = norm(An,'fro');
        edge.V1  = norm(V1);
        edge.V2  = norm(V2);
        edge.b0  = abs(b0);
    end

    % ===================== FIGURES =====================

    % (1) Histogram of values (or magnitudes if complex)
    f1 = figure('Color','w');
    histogram(valsPlot, 60); grid on
    if isReal
        xlabel('value');
    else
        xlabel('|value|');
    end
    ylabel('count');
    title(sprintf('Histogram of %d entries', numel(vals)));
    exportgraphics(f1, fullfile(outDir,'values_hist.png'), 'Resolution',600);
    close(f1);

    % (2) Counts bar
    cats   = {'positive','negative','zero','integer','single-like','double-like','complex'};
    counts = [posCnt, negCnt, zeroCnt, intCnt, singleCnt, doubleCnt, complexCnt];
    f2 = figure('Color','w');
    bar(counts); grid on
    set(gca,'XTick',1:numel(cats),'XTickLabel',cats);
    xtickangle(20);
    ylabel('count');
    title('Entry counts by sign/type');
    exportgraphics(f2, fullfile(outDir,'counts_bar.png'), 'Resolution',600);
    close(f2);

    % (3) Eigenvalue scatter (square only)
    if isSq
        f3 = figure('Color','w');
        plot(real(eigvals), imag(eigvals), 'o'); grid on
        xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
        title(sprintf('Eigenvalues (%dx%d)', m,n));
        exportgraphics(f3, fullfile(outDir,'eig_scatter.png'), 'Resolution',600);
        close(f3);
    end

    % (4) Edge metrics bar (if available)
    if ~isempty(fieldnames(edge))
        f4 = figure('Color','w');
        bar([edge.AnF, edge.V1, edge.V2, edge.b0]); grid on
        set(gca,'XTick',1:4,'XTickLabel',{'||A_n||_F','||V_1||_2','||V_2||_2','|b_0|'});
        ylabel('value');
        title('Natural Matrix Partition metrics (last row/col)');
        exportgraphics(f4, fullfile(outDir,'edge_block_metrics.png'), 'Resolution',600);
        close(f4);
    end

    % ===================== REPORT =====================
    report = struct();
    report.size        = [m n];
    report.issparse    = issp;
    report.isreal      = isReal;
    report.isSquare    = isSq;
    report.isSymmetric = isSym;
    report.symRel      = symRel;
    report.trace       = tr;
    report.diagStats   = diagStats;
    report.min3        = min3;
    report.max3        = max3;
    report.minVal      = minVal;
    report.maxVal      = maxVal;
    report.range       = rangeVal;
    report.counts      = struct('positive',posCnt,'negative',negCnt,'zero',zeroCnt, ...
                                'integer',intCnt,'singleLike',singleCnt,'doubleLike',doubleCnt, ...
                                'complex',complexCnt);
    if isSq
        report.det        = detVal;
        report.logabsdet  = logabsdet;
        report.detSign    = detSign;
        report.eigs       = eigvals;
    end
    if ~isempty(fieldnames(edge)), report.edge = edge; end

    % Text summary
    fid = fopen(fullfile(outDir,'summary.txt'),'w');
    fprintf(fid,'size=%dx%d\nsparse=%d real=%d square=%d symmetric=%d (rel=%.3e)\n', ...
            m,n,issp,isReal,isSq,isSym,symRel);
    fprintf(fid,'trace=%.16g\nrange=[%.16g, %.16g] (Δ=%.16g)\n', tr, minVal, maxVal, rangeVal);
    fprintf(fid,'diag: mean=%.6g std=%.6g min=%.6g max=%.6g\n', ...
            diagStats.mean, diagStats.std, diagStats.min, diagStats.max);
    fprintf(fid,'counts: pos=%d neg=%d zero=%d  integer=%d singleLike=%d doubleLike=%d complex=%d\n', ...
            posCnt,negCnt,zeroCnt,intCnt,singleCnt,doubleCnt,complexCnt);
    if isSq
        if isReal
            fprintf(fid,'det≈%.6e (sign=%+g, log|det|=%.6e)\n', report.det, report.detSign, report.logabsdet);
        else
            fprintf(fid,'|det|≈%.6e (log|det|=%.6e, complex matrix)\n', abs(report.det), report.logabsdet);
        end
    end
    if ~isempty(fieldnames(edge))
        fprintf(fid,'edge metrics: ||A_n||F=%.6e  ||V1||=%.6e  ||V2||=%.6e  |b0|=%.6e\n', ...
                edge.AnF, edge.V1, edge.V2, edge.b0);
    end
    fclose(fid);
end
