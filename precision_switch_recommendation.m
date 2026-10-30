function rec = precision_switch_recommendation(AorK, opts)
% PRECISION_SWITCH_RECOMMENDATION  Rule-of-thumb guidance for precision choice.
%
%   rec = precision_switch_recommendation(A)           % uses cond(A) (2-norm)
%   rec = precision_switch_recommendation(kappa)       % you pass κ directly
%   rec = precision_switch_recommendation(..., opts)
%
% opts fields (all optional):
%   .normType   : '2'|'fro'  (default '2')  % which κ to compute if A is given
%   .kappaFrob  : []         (if you already have κ_F; saves recompute)
%   .kappa2     : []         (if you already have κ_2; saves recompute)
%   .thresholds : [t1 t2 t3] (default [1e8 1e12 1e14])
%
% Output struct rec:
%   .kappa_used       : scalar κ used for decision
%   .digits_reliable  : ~max(0, 16 - log10(κ_2))  (informative)
%   .bucket           : 'safe'|'watch'|'risky'
%   .advice           : short text recommendation
%   .detail           : cellstr of bullet-point actions
%
% Example:
%   A = randn(100); rec = precision_switch_recommendation(A)
%   fprintf('%s\n', rec.advice); disp(string(rec.detail))
%
% Note: if you pass κ directly, it’s treated as κ_2 for the digits estimate.

arguments
    AorK
    opts.normType   (1,1) string = "2"
    opts.kappaFrob  double = []
    opts.kappa2     double = []
    opts.thresholds (1,3) double = [1e8 1e12 1e14]
end

% --- decide kappa
if isscalar(AorK) && isnumeric(AorK)
    kappa2 = AorK;
    kappa_used = kappa2;
else
    A = AorK;
    switch lower(opts.normType)
        case "2"
            kappa2 = i_get_kappa2(A, opts.kappa2);
            kappa_used = kappa2;
        case "fro"
            kappaF = i_get_kappaF(A, opts.kappaFrob);
            % still compute kappa2 cheaply for digits estimate
            kappa2 = i_get_kappa2(A, opts.kappa2);
            kappa_used = kappaF;
        otherwise
            error('normType must be "2" or "fro".');
    end
end

% --- digits heuristic based on κ2
digits = max(0, 16 - log10(max(realmin, kappa2)));

% --- thresholds
t1 = opts.thresholds(1); t2 = opts.thresholds(2); t3 = opts.thresholds(3);

if kappa_used < t1
    bucket = "safe";
    advice = "Double is fine (~8+ digits expected).";
    detail = { ...
        "Use stable factorizations (chol/qr/svd) as appropriate.", ...
        "Basic residual check is sufficient." ...
    };
elseif kappa_used < t2
    bucket = "watch";
    advice = "Borderline: monitor residuals; consider iterative refinement.";
    detail = { ...
        "Solve in double with QR/SVD; compute residual r = b - A*x.", ...
        "If r stagnates, try one or two refinement steps.", ...
        "Consider scaling (row/col equilibrate) before solve." ...
    };
elseif kappa_used < t3
    bucket = "risky";
    advice = "High κ: expect accuracy loss; recondition or refine aggressively.";
    detail = { ...
        "Apply your structure-preserving edits / scaling first.", ...
        "Prefer SVD-based solve; use iterative refinement.", ...
        "If accuracy still poor, try mixed-precision with higher-precision accumulation." ...
    };
else
    bucket = "very_risky";
    advice = "Very high κ: double likely insufficient; recondition or use higher precision.";
    detail = { ...
        "Run your reconditioning (Theo 6.2 edits) or rescale; re-evaluate κ.", ...
        "If κ remains large, switch to VPA/quad or mixed-precision solvers.", ...
        "Document loss of reliable digits to set expectations." ...
    };
end

rec = struct('kappa_used',kappa_used, ...
             'digits_reliable',digits, ...
             'bucket',bucket, ...
             'advice',advice, ...
             'detail',{detail});
end

% ---------------- helpers ----------------

function k2 = i_get_kappa2(A, k2_opt)
    if ~isempty(k2_opt), k2 = k2_opt; return; end
    % exact 2-norm condition via SVD (robust)
    s = svd(A);
    k2 = s(1)/s(end);
end

function kF = i_get_kappaF(A, kF_opt)
    if ~isempty(kF_opt), kF = kF_opt; return; end
    s = svd(A);
    kF = sqrt(sum(s.^2)) * sqrt(sum((1./s).^2));  % ||A||_F * ||A^{-1}||_F
end
