function plotLogBarChart()
% Auto-generating log-scale bar chart with LaTeX labels and title
% based on matrix size, condition number, and norm type.

    % === MATRIX SETUP ===
    d = 100;                  % Dimension
    M1 = rand(d);
    M2 = rand(d);
    condM1 = cond(M1);
    condM2 = cond(M2);
    normType = 'fro';

    % === MEASUREMENTS ===
    normM1     = norm(M1, normType);
    normM1_inv = norm(inv(M1), normType);
    normM2     = norm(M2, normType);
    normM2_inv = norm(inv(M2), normType);

    scalars = [condM1, normM1, normM1_inv, condM2, normM2, normM2_inv];

    % === LABELS ===
    labels = { ...
        'cond(M_1)', ...
        '$$\|M_1\|$$', ...
        '$$\|M_1^{-1}\|$$', ...
        'cond(M_2)', ...
        '$$\|M_2\|$$', ...
        '$$\|M_2^{-1}\|$$' ...
    };

    groupColors = [1 1 1 2 2 2];  % M1 = 1, M2 = 2
    colors = [0.2 0.4 0.8; 0.85 0.2 0.2];  % Blue / Red

    % === TITLE INFO (auto-generated)
    timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    titleLine1 = '\bf\fontsize{16}Log-Scale Condition Numbers and Norms';
    titleLine2 = sprintf('\\fontsize{14}Matrix Size: %dx%d   |   cond(M_1) = %.2e   |   cond(M_2) = %.2e', ...
                         d, d, condM1, condM2);
    titleLine3 = sprintf('\\fontsize{12}Norm Type: %s   |   Generated: %s', ...
                         normType, timestamp);

    % === FILENAME (auto)
    filenameBase = sprintf('M%d_CN1_%.1e_CN2_%.1e', d, condM1, condM2);
    filenameBase = regexprep(filenameBase, '[^\w]', '_');  % Sanitize filename

    chartFolder = 'Charts';
    if ~exist(chartFolder, 'dir'); mkdir(chartFolder); end

    % === CREATE FIGURE
    tiledlayout(1,1);
    ax = nexttile;
    b = barh(scalars, 'FaceColor', 'flat');

    % === LOG SCALE X-AXIS
    safeScalars = scalars(scalars > 0);
    xMin = 10^floor(log10(min(safeScalars))) / 10;
    xMax = 10^ceil(log10(max(safeScalars))) * 10;
    set(gca, 'XScale', 'log');
    xlim([xMin, xMax]);

    % === REVERSE Y AXIS (top-down order)
    set(gca, 'YDir', 'reverse');

    % === COLOR CODE
    for i = 1:numel(scalars)
        b.CData(i,:) = colors(groupColors(i), :);
    end

    % === GRID
    grid on;
    ax.XMinorGrid = 'on';
    ax.MinorGridAlpha = 0.3;
    ax.MinorGridLineStyle = ':';

    % === AXIS LABELS
    xlabel('Value (log scale)', 'FontSize', 14, 'FontName', 'Arial');
    set(gca, 'YTick', 1:numel(labels), ...
             'YTickLabel', labels, ...
             'TickLabelInterpreter', 'latex', ...
             'FontSize', 12, ...
             'FontName', 'Arial');

    % === MARGIN FIX (prevent cutoff)
    ax.Position = [0.2 0.1 0.7 0.8];

    % === TITLE (3 lines)
    title({titleLine1; titleLine2; titleLine3}, 'Interpreter', 'tex');

    % === VALUE LABELS
    for i = 1:numel(scalars)
        text(scalars(i)*1.3, i, ...
             sprintf('%.1e', scalars(i)), ...
             'HorizontalAlignment','left', ...
             'VerticalAlignment','middle', ...
             'FontSize', 10, ...
             'FontName', 'Arial');
    end

    % === LEGEND (manual dummy bars)
    hold on;
    hM1 = barh(nan, nan, 'FaceColor', colors(1,:));
    hM2 = barh(nan, nan, 'FaceColor', colors(2,:));
    legend([hM1, hM2], {'M1 group', 'M2 group'}, ...
           'Location', 'southoutside', 'Orientation', 'horizontal');
    hold off;

    % === EXPORT
    pngFile = fullfile(chartFolder, [filenameBase '.png']);
    pdfFile = fullfile(chartFolder, [filenameBase '.pdf']);

    saveas(gcf, pngFile);
    print(gcf, pdfFile, '-dpdf', '-bestfit');

    fprintf('✅ Chart saved:\n- %s\n- %s\n', pngFile, pdfFile);
end
