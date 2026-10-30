% === Input Data ===
scalars = [1e2, 10, 0.1, 1e5, 50, 0.0001];
labels = { ...
    'cond(M_1)', ...
    '$$\|M_1\|$$', ...
    '$$\|M_1^{-1}\|$$', ...
    'cond(M_2)', ...
    '$$\|M_2\|$$', ...
    '$$\|M_2^{-1}\|$$' ...
    };

groupColors = [1 1 1 2 2 2];  % M1 = 1, M2 = 2
colors = [0.2 0.4 0.8; 0.85 0.2 0.2];  % Blue, Red

% === Create figure ===
figure('Color', 'w');

% Horizontal bar plot
b = barh(scalars, 'FaceColor', 'flat');

% Log scale on X-axis
set(gca, 'XScale', 'log');

% Title with large font
title('Log-Scale Condition Numbers and Norms | ID: 250805 | UserID: REG | DIM: 10 | CN: 46', ...
      'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');

% X label
xlabel('Value (log scale)', 'FontSize', 14, 'FontName', 'Arial');

% Tick labels using LaTeX
set(gca, 'YTick', 1:numel(labels), ...
         'YTickLabel', labels, ...
         'TickLabelInterpreter', 'latex', ...
         'FontSize', 12, ...
         'FontName', 'Arial');

% Gridlines
grid on;
ax = gca;
ax.XMinorGrid = 'on';
ax.YGrid = 'on';

% Apply color coding
for i = 1:numel(scalars)
    b.CData(i,:) = colors(groupColors(i), :);
end

% Add value labels
for i = 1:numel(scalars)
    text(scalars(i)*1.3, i, ...
         sprintf('%.1e', scalars(i)), ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','middle', ...
         'FontSize', 10, ...
         'FontName', 'Arial');
end

% Legend (manual dummy bars)
hold on;
hM1 = barh(nan, nan, 'FaceColor', colors(1,:));
hM2 = barh(nan, nan, 'FaceColor', colors(2,:));
legend([hM1 hM2], {'M1 group', 'M2 group'}, ...
       'Location', 'southoutside', 'Orientation', 'horizontal');
hold off;

% === Save to disk ===
chartFolder = 'Charts';
if ~exist(chartFolder, 'dir'); mkdir(chartFolder); end
saveas(gcf, fullfile(chartFolder, 'Scalars_LogOnly_Horizontal.png'));
print(gcf, fullfile(chartFolder, 'Scalars_LogOnly_Horizontal.pdf'), '-dpdf', '-bestfit');
