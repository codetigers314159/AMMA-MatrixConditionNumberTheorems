
function test_analyzeMatrixBasics()
% Quick smoke tests for analyzeMatrixBasics (generic matrices only).

% ---- PATH SETUP (adjust if needed) ----
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);                  % current folder
addpath(genpath(thisDir));         % and any subfolders

% Sanity: verify the analyzer is visible
if exist('analyzeMatrixBasics','file') ~= 2
    error(['Could not find analyzeMatrixBasics.m on the path.\n' ...
           'Make sure it is in the same folder as this test, or addpath() to the folder.']);
end

rng(42);

CASES = {
  struct('name','real_rect_50x40',    'M', randn(50,40)), 
  struct('name','complex_square_60',  'M', randn(60) + 1i*randn(60)), 
  struct('name','ill_conditioned_80', 'M', makeIllCond(80, 1e10)), 
  struct('name','integer_like_40',    'M', randi([-9,9],40) + 1e-10*randn(40)), 
  struct('name','mixed_scale_70x30',  'M', 1e6*randn(70,30)), 
  struct('name','zeros_heavy_90x20',  'M', sprinkleZeros(randn(90,20), 0.15))
};

OUTROOT = fullfile('figs','matrix_smoke_tests');
if ~exist(OUTROOT,'dir'), mkdir(OUTROOT); end

fprintf('Running %d analyzeMatrixBasics tests...\n', numel(CASES));

for k = 1:numel(CASES)
    name = CASES{k}.name;
    M    = CASES{k}.M;
    out  = fullfile(OUTROOT, name);
    try
        rep = analyzeMatrixBasics(M, out);
        sz  = rep.size;
        fprintf('  ✓ %-24s  size=%dx%d  real=%d  square=%d  symmetric=%d\n', ...
                name, sz(1), sz(2), rep.isreal, rep.isSquare, rep.isSymmetric);
    catch ME
        fprintf('  ✗ %-24s  ERROR: %s\n', name, ME.message);
    end
    pause(0.1); close all;
end

fprintf('\nDone. See output folders under %s\n', OUTROOT);
end

% ---------- helpers (local-only) ----------
function M = makeIllCond(n, targetCond)
[U,~] = qr(randn(n));
[V,~] = qr(randn(n));
s = logspace(0, -log10(targetCond), n);  % singular values 1 ... 1/cond
M = U*diag(s)*V';
end

function A = sprinkleZeros(A, frac)
if nargin<2, frac = 0.1; end
mask = rand(size(A)) < frac;
A(mask) = 0;
end


