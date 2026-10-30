 % MATLAB Script: Full Matrix Analysis (Eigenvalues, SVD)
 % RAUL GARCIA, Full Analysis: Log Plots, Eigenvalues, SVD
 clear 
 clc
 large = 1e6;  % Threshold for condition number
 d = 100;
 trials = 10000;
 foundLargeCn = false;
 fprintf('Searching for a %dx%d matrix with condition number > %.0f\n', d, d, large);
 for x = 1:trials
    m = rand(d);
    cn = cond(m, 'fro');
    matStoreCn(x) = cn;

   
    if cn > large && ~foundLargeCn
        foundLargeCn = true;
        condMatrix = m;
        condNumber = cn;
        fprintf('\n>>> Found matrix on trial #%d\n', x);
        fprintf('Condition number: %.3e\n', condNumber);
        % === Display top-left preview ===
        disp('Top-left 5x5 preview of condMatrix:');
        disp(condMatrix(1:5,1:5));
        
        % === Matrix statistics ===
        fprintf('Min: %.3e, Max: %.3e\n', min(condMatrix(:)), max(condMatrix(:)));
        fprintf('Mean: %.3e, Std Dev: %.3e\n', mean(condMatrix(:)), std(condMatrix(:)));
     
        % === Heatmap visualization ===
        figure;
        imagesc(condMatrix);
        colorbar;
        title('Heatmap of condMatrix');
        xlabel('Columns'); ylabel('Rows');

        % === Heatmap visualization log10 ===
        figure;
        imagesc(log10(condMatrix));
        colorbar;
        title('Heatmap of Log 10 of condMatrix');
        xlabel('Columns'); ylabel('Rows');

        % === Eigenvalue plot ===
        eigvals = eig(condMatrix);
        figure;
        axPlot=axes; 
        plot(real(eigvals), imag(eigvals), 'o');
        axPlot.XGrid = 'on';
        axPlot.YGrid = 'on';
        title('Eigenvalues of condMatrix (Real vs Imag)');
        xlabel('Real Part'); ylabel('Imaginary Part');
        
        % === Eigenvalue norms plot ===
        eigNorms = log10(abs(eigvals)); % Compute norm (magnitude)
        figure;
        histogram(eigNorms);
        title('Ln of Distribution of condMatrix - Eigenvalue Norms');
        xlabel('Ln Eigenvalues'); ylabel('Count');

        % === Eigenvalue angles plot ===
        eigAngles = angle(eigvals); % Angle in radians; % Compute norm (magnitude)
        figure;
        axPlot = axes;
        plot(eigAngles, 'o'); % Plot norm vs index
        title('Angles of Eigenvalues (Radians) of condMatrix');
        xlabel('Index');
        ylabel('Angles of Eigenvalues');
        grid on;
      
        eigAngles = angle(eigvals)*(180/pi); % Angle in radians; % Compute norm (magnitude)
        figure;
        axPlot = axes;
        plot(eigAngles, 'o'); % Plot norm vs index
        title('Angles of Eigenvalues (Degrees) of condMatrix');
        xlabel('Index');
        ylabel('Angles of Eigenvalues');
        grid on;

        figure;
        scatter(real(eigvals), imag(eigvals), 36, angle(eigvals), 'filled');
        colorbar;
        title('Eigenvalues Colored by Phase of condMatrix');
        xlabel('Real Part');
        ylabel('Imaginary Part');

        % === SVD singular value plot ===
        s = svd(condMatrix);
        figure;
        axPlot=axes; 
        semilogy(s, 'o-');
        axPlot.XGrid = 'on';
        axPlot.YGrid = 'on';
        title('Singular Values of condMatrix');
        xlabel('Index'); ylabel('Singular Value (log scale)');

        % === histogram ===
        figure;
        axPlot=axes;
        histogram(condMatrix(:), 'NumBins',d , 'Normalization', 'pdf');
        axPlot.YGrid = 'on';
        title('Normalized Bins Histogram of condMatrix');

        figure;
        histogram(condMatrix(:));
        title('Distribution of Matrix Values - Histogram of condMatrix');
        xlabel('Value'); ylabel('Frequency');

        figure;
        surf(condMatrix);
        title('Surface Plot of Matrix of of condMatrix');
        xlabel('Column'); ylabel('Row'); zlabel('Value');

        figure;
        surf(log10(condMatrix));
        title('Surface Log 10 Plot of condMatrix');
        xlabel('Column'); ylabel('Row'); zlabel('Value');

        figure;
        vals = condMatrix(:);
        [sortedVals, idx] = sort(vals, 'ascend');
        bar(sortedVals(1:10));
        title('Bottom 10 Smallest Matrix Elements of condMatrix');

        figure;
        boxplot(condMatrix(:));
        title('Boxplot of Matrix Elements of condMatrix');

        figure;
        boxplot(condMatrix, 'Notch', 'on');
        title('Boxplot of Matrix Elements Notch compare Medians of condMatrix');

        % === Export matrix to file ===
        save('LargeCondMatrix.mat', 'condMatrix', 'condNumber', 'eigvals', 's');
        writematrix(condMatrix, 'LargeCondMatrix.csv');
        fprintf('Matrix exported to LargeCondMatrix.mat and .csv\n');

        break; %there may be more & larger?
    end
 end
 if ~foundLargeCn
    fprintf('\nNo matrix found with condition number > %.0f in %d trials.\n', large,  trials);
 end
 % === Log10 Histogram of Condition Numbers ===
 figure;
 axPlot=axes;
 histogram(log10(matStoreCn), 50);
 axPlot.YGrid = 'on';
 title('Log10 Histogram of Condition Numbers of condMatrix');
 xlabel('log10(Condition Number)'); ylabel('Frequency');
 % === Semilogy Plot ===
 figure;
 axPlot=axes;
 semilogy(matStoreCn, '.');
 axPlot.XGrid = 'on';
 axPlot.YGrid = 'on';
 title('Condition Numbers Over Trials (log scale) of condMatrix');
 xlabel('Trial #'); ylabel('Condition Number');
% --- Find the largest value and its index ---
[maxVal, maxIdx] = max(matStoreCn);

% --- Add text annotation at that point ---
hold on;
semilogy(maxIdx, maxVal, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); % Highlight point
text(maxIdx, maxVal, sprintf('Max = %.2e  ', maxVal), ...
    'VerticalAlignment', 'middle', ...
    'HorizontalAlignment', 'right', ...  % ← shifts text left of point
    'FontWeight', 'bold', ...
    'Color', 'red');
