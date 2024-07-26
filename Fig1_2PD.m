load(fullfile(DataPath(), '2PD_processed'))
breast_image = imread(fullfile(DataPath(), 'BreastDrawings.png'));

%% Make figure
sf = GetSigmoid(2);
SetFont('Arial', 9)
clf;
set(gcf, 'Units', 'Inches', 'Position', [30 1 6.45 5])
xl = [-30 30];
xq = linspace(xl(1), xl(2));

% Examples - Hand
yp1 = 0.625; yp2 = 0.45; yp3 = 0.1;
axes('Position', [0.1 yp1 .14 .5]);
    imsize = 1150; x = 1450; y = 1850;
    imshow(breast_image([y:y+imsize], [x:x+imsize], :))

axes('Position', [0.32 yp1 .14 .5]);
    imsize = 1000; x = 150; y = 1900;
    imshow(breast_image([y:y+imsize], [x:x+imsize], :))

axes('Position', [0.55 yp1 .14 .5]);
    imsize = 1000; x = 1700; y = 500;
    imshow(breast_image([y:y+imsize], [x:x+imsize], :))

axes('Position', [0.76 yp1 .15 .5]);
    imsize = 1000; x = 300; y = 550;
    imshow(breast_image([y:y+imsize], [x:x+imsize+100], :))

xp = [.1, .3167, .5334, .75];
i = 6; % Example subject idx
t = {'Hand', 'Back', 'L.Breast', 'M.Breast'};
for j = 1:4
    axes('Position', [xp(j) yp2 .15 .25]); hold on
    tab = ResponseTable_ConditionMean(subjectData(i).(lf{j}).data);
    plot(xq, sf(subjectData(i).(lf{j}).SigCoeffs, xq), 'Color', [.6 .6 .6], 'LineWidth', 2);
    scatter(tab.distance, tab.response, 40, lf_colors(j,:), 'MarkerFaceColor', lf_colors(j,:))
    % Formatting
    set(gca, 'YLim', [0 1], 'YTick', [0 .5 1], 'XLim', xl)
    [tx,ty] = GetAxisPosition(gca, 5,95);
    text(tx,ty, sprintf('JND =\n%0.1f mm', jnd_table{i, lf{j}}), ...
        'VerticalAlignment', 'top', 'Color', 'k')
    if j == 1
        ylabel('p(above)')
    end
    text(0,1.1, t{j}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
end

% JND distributions
axes('Position', [0.1 yp3 .2 .225]); hold on
    for l = 1:length(lf)
        Swarm(l, jnd_table{:,l}, lf_colors(l,:), 'DS', 'Box', 'SWR', 0.5, 'DW', 0.3);
    end
    set(gca, 'XLim', [0.5, 4.5], ...
             'YLim', [0 20], ...
             'XTick', [1:4], ...
             'XTickLabel', {'Hand', 'Back', 'L.Breast', 'M.Breast'}, ...
             'YTick', [0:10:20])
    ylabel('JND (mm)')

% Breast scaling
axes('Position', [0.4 yp3 .2 .225]); hold on
    % Best fit line
    x = meas_table.delta_bust;
    y = jnd_table.breast;
    nan_idx = isnan(x) | isnan(y);
    x = x(~nan_idx); y = y(~nan_idx);
    p1 = polyfit(x, y, 1);
    [r,p] = corr(x,y);
    xl = [0, 10];
    plot(xl, polyval(p1, xl), 'Color', [.6 .6 .6], 'LineStyle', '--')
    scatter(x, y, 30, lf_colors(3,:), 'MarkerFaceColor', lf_colors(3,:))
    % Formatting
    ylabel('Lat Breast JND (mm)')
    xlabel(sprintf('%s Bust', GetUnicodeChar('Delta')))
    set(gca, 'XLim', xl, 'XTick', [0:5:10], 'YLim', [0 20], 'YTick', [0:10:20])
    [tx,ty] = GetAxisPosition(gca, 5,95);
    text(tx,ty, sprintf('r = %0.3f\n%s', r, pStr(p)), 'VerticalAlignment', 'top')

% NAC scaling
axes('Position', [0.7 yp3 .2 .225]); hold on
    % Best fit line
    x = meas_table.bust;
    y = jnd_table.NAC;
    nan_idx = isnan(x) | isnan(y);
    x = x(~nan_idx); y = y(~nan_idx);
    p1 = polyfit(x, y, 1);
    [r,p] = corr(x,y);
    xl = [30, 45];
    plot(xl, polyval(p1, xl), 'Color', [.6 .6 .6], 'LineStyle', '--')
    scatter(x, y, 30, lf_colors(4,:), 'MarkerFaceColor', lf_colors(4,:))
    % Formatting
    ylabel('Med Breast JND (mm)')
    xlabel('Bust (inches)')
    % set(gca, 'XLim', xl, 'XTick', [0:5:10], 'YLim', [0 20], 'YTick', [0:10:20])
    [tx,ty] = GetAxisPosition(gca, 95,5);
    text(tx,ty, sprintf('r = %0.3f\n%s', r, pStr(p)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')

% Annotation
annotation("textbox", [0.45 yp2-0.14, 0.1 0.1], 'String', 'Distance (mm)', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'EdgeColor', 'none')
char_offset = 64;
annotation("textbox", [0.05 0.925 .05 .05], 'String', char(char_offset+1), ...
'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.05 0.35 .05 .05], 'String', char(char_offset+2), ...
'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.35 0.35 .05 .05], 'String', char(char_offset+3), ...
'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.65 0.35 .05 .05], 'String', char(char_offset+4), ...
'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
shg
print(gcf, "C:\Users\somlab\Downloads\Fig1_2PD.png", '-dpng', '-r300')
return

%% Plot all combinations of JND & measurement
clf; i = 1;
for x = 1:width(meas_table)
    for y = 1:width(jnd_table)
        subplot(width(meas_table), width(jnd_table), i); hold on
        nan_idx = isnan(meas_table{:,x}) | isnan(jnd_table{:,y});
        p1 = polyfit(meas_table{~nan_idx,x}, jnd_table{~nan_idx,y}, 1);
        [r,p] = corr(meas_table{~nan_idx,x}, jnd_table{~nan_idx,y});
        lims = [min(meas_table{~nan_idx,x}), max(meas_table{~nan_idx,x})];
        lims = [lims(1) - range(lims)*0.05, lims(2) + range(lims)*0.05];
        plot(lims, polyval(p1, lims), 'Color', [.6 .6 .6], 'LineStyle', '--')
        scatter(meas_table{~nan_idx,x}, jnd_table{~nan_idx,y}, 50, [.6 .6 .6], 'filled')
        if y == 1
            ylabel(strrep(meas_table.Properties.VariableNames{x}, '_', ' '), 'FontWeight', 'bold')
        end
        if x == 1
            title(jnd_table.Properties.VariableNames{y})
        end
        [tx,ty] = GetAxisPosition(gca, 5,95);
        text(tx,ty, sprintf('r = %0.3f\n%s', r, pStr(p*10)), 'VerticalAlignment', 'top')
        i = i + 1;
    end
end

shg

%% Plot single subject
i = 11;
sf = GetSigmoid(2);
clf;
t = tiledlayout('flow');
for j = 1:length(lf)
    nexttile(); hold on
    % Summarize
    tab = ResponseTable_ConditionMean(subjectData(i).(lf{j}).data);
    % Plot
    scatter(tab.distance, tab.response, 40, [.6 .6 .6])
    xq = linspace(min(tab.distance), max(tab.distance));
    plot(xq, sf(subjectData(i).(lf{j}).SigCoeffs, xq), 'Color', [.6 .6 .6]);
    % Formatting
    title(lf{j})
    set(gca, 'YLim', [0 1], 'YTick', [0 1])
end

shg