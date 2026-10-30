function plotLogBarChart()
% Log-scale bar chart with integer labels and clean text formatting
% Based on matrix condition numbers and norms (Theorem 6.2 format)

    % === MATRIX SETUP ===
    d = 100;                  
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

    % === CLEAN LABELS (no LaTeX, just plain || ||)
    labels = {
        'cond(M1)', ...
        '||M1||', ...
        '||M1^{-1}||', ...
        'cond(M2)', ...
        '||M2||', ...
        '||M2^{-1}||'
    };

    groupColors = [1 1 1 2 2 2];  
    colors = [0.2 0.4 0.8; 0.85 0.2 0.2];  

    % === TITLE INFO
    timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    titleLine1 = '\bf\fontsize{16}Theorem 6.2 – Log-Scale Condition Numbers and Norms';
    titleLine2 = sprintf('\\fontsize{14}Matrix Size: %dx%d   |   cond(M1) = %d   |   cond(M2) = %d', ...
                         d, d, round(condM1), round(condM2));
    titleLine3 = sprintf('\\fontsize{12}Norm Type: %s   |   Generated: %s', ...
                         normType, timestamp);

    % === FILE NAME
    filenameBase = sprintf('M%d_CN1_%d_CN2_%d', d, round(condM1), round(condM2));
    filenameBase = regexprep(filenameBase, '[^\w]', '_');

    chartFolder = 'Charts';
    if ~exist(chartFolder, 'dir'); mkdir(chartFolder); end

    % === CREATE FIGURE & AXES
    figure;
    ax = axes('Position', [0.2 0.1 0.7 0.8]);
    b = barh(ax, scalars, 'FaceColor', 'flat');

    % === LOG SCALE
    safeScalars = scalars(scalars > 0);
    xMin = 10^floor(log10(min(safeScalars))) / 10;
    xMax = 10^ceil(log10(max(safeScalars))) * 10;
    set(ax, 'XScale', 'log');
    xlim(ax, [xMin, xMax]);
    set(ax, 'YDir', 'reverse');

    % === COLOR ASSIGNMENT
    for i = 1:numel(scalars)
        b.CData(i,:) = colors(groupColors(i), :);
    end

    % === GRID
    grid(ax, 'on');
    ax.XMinorGrid = 'on';
    ax.MinorGridAlpha = 0.3;
    ax.MinorGridLineStyle = ':';

    % === AXIS LABELS
    xlabel(ax, 'Value (log scale)', 'FontSize', 14, 'FontName', 'Arial');
    set(ax, 'YTick', 1:numel(labels), ...
            'YTickLabel', labels, ...
            'TickLabelInterpreter', 'tex', ...
            'FontSize', 12, ...
            'FontName', 'Arial');

    % === TITLES (3 lines)
    title(ax, {titleLine1; titleLine2; titleLine3}, 'Interpreter', 'tex');

    % === VALUE LABELS (integers only)
    for i = 1:numel(scalars)
        labelStr = sprintf('%d', round(scalars(i)));
        text(ax, scalars(i)*1.3, i, ...
             labelStr, ...
             'HorizontalAlignment','left', ...
             'VerticalAlignment','middle', ...
             'FontSize', 10, ...
             'FontName', 'Arial');
    end

    % === LEGEND
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
