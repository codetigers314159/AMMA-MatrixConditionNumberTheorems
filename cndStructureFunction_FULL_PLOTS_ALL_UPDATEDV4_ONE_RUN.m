
% ---------------------------------------------------------------
% Author: Raul E. Garcia

% Script: cndStructureFunction_FULL_PLOTS_ALL_UPDATEDV4_ONE_RUN.m
% UPDATED: 9/21/25 REG
% Purpose: Full matrix decomposition + all plots in one file
% Includes: Run prompt, titles with run#, cond#, norm, dim - ONE RUN 
% Random: Allows random number in a range - user defined - def is 0 to 1
% Verified: Matrix & submatricesa contents
% ---------------------------------------------------------------
clear; clc;

% === User Input Parameters ===
d=0;
% 
% SAMPLE
% prompt = {'Enter matrix size:','Enter colormap name:'};
% dlgtitle = 'Input';
% fieldsize = [1 45; 1 45];
% definput = {'20','hsv'};
% answer = inputdlg(prompt,dlgtitle,fieldsize,definput)

%using inputdlg(prompt,dlgtitle,fieldsize,definput)
prompt = {'Random Square Matrix size (dim): ','Max# of Trials (to find large matrices)','Min Condition Number: ','Max Condition Number: ','Min Matrix#','Max Matrix#','User Id or Run Id (5 chars):'};
dlgtitle = 'Input Theo 6.2 (Theorem 1 on SIAM Paper) Cond Numbers';
fieldsize = [ 1 70; 1 70; 1 70; 1 70;1 70;1 70;1 70;];
definput = {'10','10000000','10000','1000000','0','1','SIAM'};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput)

d=str2double(answer{1})
%this to allow we have enough runs to find 1 valid matrix cond#
maxTrials=str2double(answer{2})
cnMin1= str2double(answer{3})
cnMax2= str2double(answer{4})
%this allows a range of random values
MinRandMatrixNum= str2double(answer{5})
MaxRandMatrixNum= str2double(answer{6})
userId = answer{7}

normType = 'fro' % Frobenius norm only

% === Unique Run ID and Save Folder ===
% timestamp = datestr(now,'yyyymmdd_HHMMSS');
% datestamp = datestr(now,'yyyymmdd');
% chartFolder = ['Theo62_Charts_' runID];
% if ~exist(chartFolder, 'dir')
%     mkdir(chartFolder);
% end

timestamp = datestr(now,'yyyymmdd_HHMMSS');
datestamp = datestr(now,'yymmdd');
runID = ['' datestamp]
chartFolder = ['Theo6_2_Charts_' userId '_' runID ];
if ~exist(chartFolder, 'dir')
    mkdir(chartFolder);
end
disp(['Charts will be saved in: ' chartFolder]);

%rn = a + (b-a).*rand(n,1)

foundLargeCn=false

r1=MinRandMatrixNum
r2=MaxRandMatrixNum

for x = 1:maxTrials

    md = rand(d);
    m = r1 + (r2-r1)*md

    cn = cond(m, 'fro');
  
    % Print real-time progress
    fprintf('Trial #%d: cond = %.3e\n', x, cn);
 
    if cn > cnMin1 && cn < cnMax2 
        foundLargeCn = true
        condMatrix = m;
        condNumber = cn;
        fprintf('\n>>> Found matrix on trial #%d\n', x);
        fprintf('Condition number: %.3e\n', condNumber);
        disp('Matrix preview (top-left 4x4):')
        format short
        if d>4 
          disp(condMatrix(1:4,1:4));
        end
        break;
    end
 end
 % Summary
 if ~foundLargeCn
    fprintf('\nNo matrix found with condition number > %.0f in %d trials.\n', large, trials);
 else
    fprintf('\nMatrix condMatrix found with condNumber = %.3e\n', condNumber);
 end


runs=1

    normType = 'fro' 
    % Frobenius norm only
   
    %  n+1=d = matrix dimension, below we generate the matrices
    M1 = condMatrix;   % An+1
    %M2 = condMatrix;   % An

    % next we load the submatrix M1 also into the 2nd larger matrix M2
    M2(1:d-1,1:d-1) = M1(1:d-1,1:d-1); % that is M1 is contained in M2 

    % compute each matrix inverses for A
    M1i = inv(M1);
    M2i = inv(M2);
   

    % the edge vecntors for M2  
     % the edge vecntors for M2  
    V1 = M1(1:d-1,d);
    V2 = M1(d,1:d-1);
   % the edge vecntors for M2  
    % V1 = M2(1:d,d+1);
    % V2 = M2(d+1,1:d);

    % the prod of V1 and V2 a nxn matrix
    VX = V1 * V2;

    % the all important constants: k0, b0 (last diagonal in An+1, and m0) 
    k0 = V2 * M2i * V1
    b0 = M1(d,d)
    m0 = 1 / (b0 - k0)

    % Identity nxn matrix for Ln
    I = eye(d-1);
    % The inverse of An times Ln matrix into L1
    L1 = M2i * (I + m0 * VX * M2i);

    % Interemediate matrices related to edge vectors and the inverse 
    F1 = M2i * V1;
    F2 = V2 * M2i;

    format shortG %MAKE TITLES SHORTER & CLEAR

    % === Scalars & Chart Title Metadata ===
    condM1 = cond(M1,normType) %Frobenius norm
    normM1 = norm(M1,normType) %Frobenius norm

    titleInfo = ['UID: ' num2str(userId) ...
                 ' | DIM: ' num2str(d) ...
                 ' | CN: ' num2str(condM1,'%.0f') ...
                 ' | NORM: ' num2str(normM1,'%.0f') ...
                 ' | R1: ' num2str(r1,'%.0f') ...    
                 ' | R2: ' num2str(r2,'%.0f') ... 
                 ];

    % titleInfo = ['Run#' num2str(ix) ' | ' runID ...
    %     ' | dim: ' num2str(d) ...
    %     ' | cond(M1): ' num2str(condM1, '%.2e') ...
    %     ' | ||M1||: ' num2str(normM1, '%.2e')];
    % === Scalar Bar Plot ===
    scalars = [condM1, normM1, norm(M1i,normType), ...
               cond(M2,normType), norm(M2,normType), norm(M2i,normType)];
    scalars2 = [k0, b0, m0 ];
    
    labels = {'cond(M1)','||M1||','||M1⁻¹||','cond(M2)','||M2||','||M2⁻¹||'};
    
    ix=1;

    figure;

    bar(scalars); title(['CN & Norms | ' titleInfo]); ylabel('Value');
    grid on
    set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'XTickLabelRotation', 45);
    text(1:numel(scalars), scalars, compose('%.0f', scalars), 'Rotation',90, ...
        'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',8);
    saveas(gcf, fullfile(chartFolder, ['Scalars_' num2str(ix) '.png'])); 
    
    % main scalars 
    labels = {'k0','b0','m0'};
    figure;
    bar(scalars2); title(['Matrix Scalars | ' titleInfo]); ylabel('Value');
    grid on
    set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'XTickLabelRotation', 45);
    text(1:numel(scalars), scalars, compose('%.2e', scalars), 'Rotation',90, ...
        'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',8);
    saveas(gcf, fullfile(chartFolder, ['Scalars_' num2str(ix) '.png'])); 

    % === Matrix Line Plots ===
    figure; plot(M1'); title(['M1 Rows | ' titleInfo]);
    grid on
    saveas(gcf, fullfile(chartFolder, ['M1_' num2str(ix) '.png'])); 

    figure; plot(M1i'); title(['M1 Inverse | ' titleInfo]);
    grid on
    saveas(gcf, fullfile(chartFolder, ['M1i_' num2str(ix) '.png']));

    figure; plot(M2'); title(['M2 Rows | ' titleInfo]);
    grid on
    saveas(gcf, fullfile(chartFolder, ['M2_' num2str(ix) '.png'])); 

    figure; plot(M2i'); title(['M2 Inverse | ' titleInfo]);
    grid on
    saveas(gcf, fullfile(chartFolder, ['M2i_' num2str(ix) '.png']));

    % === Edge Vectors ===
    figure; plot(V1, 'b-o'); hold on; plot(V2, 'r-*');
    grid on
    title(['V1 and V2 | ' titleInfo]); legend('V1','V2');
    saveas(gcf, fullfile(chartFolder, ['V1V2_' num2str(ix) '.png']));

    % === VX Matrix ===
    figure; imagesc(VX); colorbar;
    title(['V1xV2^T | ' titleInfo]);
    saveas(gcf, fullfile(chartFolder, ['VX_' num2str(ix) '.png'])); 

    % === L1 ===
    figure; plot(L1'); title(['L1 Rows | ' titleInfo]);
    grid on
    saveas(gcf, fullfile(chartFolder, ['L1_' num2str(ix) '.png'])); 

    % === F1, F2 Vectors ===
    figure; plot(F1, 'b'); hold on; plot(F2', 'r');
    grid on
    legend('F1-M1i*V1;','F2-V2*M1i'); title(['F1 and F2 | ' titleInfo]);
    saveas(gcf, fullfile(chartFolder, ['F1F2_' num2str(ix) '.png']));

    disp(['Matrix ' num2str(ix) ' completed.']);


disp(['All charts saved to folder: ' chartFolder]);
