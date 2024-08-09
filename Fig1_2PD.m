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
        Swarm(l, jnd_table{:,l}, lf_colors(l,:), 'DS', 'Box', 'SWR', 0.5, 'DW', 0.3, 'ShowStats', true);
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
    x = meas_table.delta_bust .* 25.4;
    y = jnd_table.breast;
    nan_idx = isnan(x) | isnan(y);
    x = x(~nan_idx); y = y(~nan_idx);
    p1 = polyfit(x, y, 1);
    [r,p] = corr(x,y, 'Type', 'Pearson'); p = p * 4; % Bonferonni correction
    xl = [0, 250];
    plot(xl, polyval(p1, xl), 'Color', [.6 .6 .6], 'LineStyle', '--')
    scatter(x, y, 30, lf_colors(3,:), 'MarkerFaceColor', lf_colors(3,:))
    % Formatting
    ylabel('Lat Breast JND (mm)')
    xlabel(sprintf('%s Bust (mm)', GetUnicodeChar('Delta')))
    set(gca, 'XLim', xl, 'XTick', [0:100:200], 'YLim', [0 20], 'YTick', [0:10:20])
    [tx,ty] = GetAxisPosition(gca, 5,95);
    text(tx,ty, sprintf('r = %0.3f\n%s', r, pStr(p)), 'VerticalAlignment', 'top')

% NAC scaling
axes('Position', [0.7 yp3 .2 .225]); hold on
   % Best fit line and scatter
   x = jnd_table.hand;
   y = jnd_table.NAC;
    nan_idx = isnan(x) | isnan(y);
    x = x(~nan_idx); y = y(~nan_idx);
    p1 = polyfit(x, y, 1);
    [r,p] = corr(x,y, 'Type', 'Pearson'); p = p * 4;
    xl = [0, 5];
    plot(xl, polyval(p1, xl), 'Color', [.6 .6 .6], 'LineStyle', '--')
    scatter(x, y, 30, [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6])
    % Formatting
    ylabel('Med Breast JND (mm)')
    xlabel('Hand JND (mm)')
    set(gca, 'XLim', xl, 'XTick', [0:2:6], 'YLim', [0 20], 'YTick', [0:5:20])
    [tx,ty] = GetAxisPosition(gca, 5,95);
    text(tx,ty, sprintf('r = %0.3f\n%s', r, pStr(p)), 'VerticalAlignment', 'top')

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
print(gcf, fullfile(FigurePath, "Fig1_2PD.png"), '-dpng', '-r300')


%% JND Stats
[p, tab, st] = anova1(table2array(jnd_table), [], "off");
% ANOVA effect
fprintf('1-way ANOVA %s\n', pStr(p))
[c,m,~,~] = multcompare(st,'display','off');

% Relative effect size
temp_jnd = table2array(jnd_table);
nan_idx = any(isnan(temp_jnd), 2);
temp_jnd = temp_jnd(~nan_idx, :);
[p,tab,s] = anova2(temp_jnd, 1, "off");
fprintf('Location explained variance: %0.2f\n', tab{2,2} / tab{end,2})
fprintf('Subject explained variance: %0.2f\n', tab{3,2} / tab{end,2})


%% Subject rank similarity
temp_jnd = table2array(jnd_table);
nan_idx = any(isnan(temp_jnd), 2);
temp_jnd = temp_jnd(~nan_idx, :);
jnd_ranks = NaN(size(temp_jnd));
for i = 1:size(jnd_ranks, 2)
    jnd_ranks(:,i) = tiedrank(temp_jnd(:,i));
end
mean_rank = mean(jnd_ranks, 2, 'omitnan');
[~, sort_idx] = sort(mean_rank);
% Sort both
mean_rank = mean_rank(sort_idx);
jnd_ranks = jnd_ranks(sort_idx, :);
[r,p] = corr(mean_rank, jnd_ranks, 'Type', 'Kendall');
p = p .* size(jnd_ranks, 2);

return

%% Supplementary figure 1
clf; 
set(gcf, 'Units', 'Inches', 'Position', [30 1 6.45 3.5])

% JND vs Delta Bust
x = 5; yi = [1,2,4];
lims = [40, 250];
xstart = [.1, .42, .75];
ypos = .55;
x_w = .2;
y_h = .35;
for i = 1:3
    axes('Position', [xstart(i), ypos, x_w, y_h]); hold on
    y = yi(i);
    % Best fit line and scatter
    tx = meas_table{:,x} .* 25.4;
    ty = jnd_table{:,y};
    nan_idx = isnan(tx) | isnan(ty);
    tx = tx(~nan_idx,:); ty = ty(~nan_idx,:);
    p1 = polyfit(tx, ty, 1);
    plot(lims, polyval(p1, lims), 'Color', [.6 .6 .6], 'LineStyle', '--')
    scatter(tx, ty, 50, lf_colors(yi(i),:), 'filled')
    % Formatting
    if i == 2
        xlabel(sprintf('%s Bust (mm)', GetUnicodeChar('Delta')))
    end
    ylabel(sprintf('%s JND (mm)', t{yi(i)}))
    % Stats text
    [r,p] = corr(tx, ty, 'Type', 'Pearson');
    p = p * 4; % Bonferroni correction
    [tx,ty] = GetAxisPosition(gca,100,100);
    text(tx,ty, sprintf('r = %0.3f\n%s', r, pStr(p)), ...
        'VerticalAlignment', 'top','HorizontalAlignment', 'right')
end

xs = [.1];
ypos = .125;
x_w = .1;
y_h = .25;
x_m = 0.085;
% Region x Region JNDs
xi = [1,1,2,2,3];
yi = [2,3,3,4,4];
for i = 1:5
    axes('Position', [xs, ypos, x_w, y_h]); hold on
    x = xi(i);
    y = yi(i);
    % Best fit line and scatter
    tx = meas_table{:,x};
    ty = jnd_table{:,y};
    nan_idx = isnan(tx) | isnan(ty);
    tx = tx(~nan_idx,:); ty = ty(~nan_idx,:);
    p1 = polyfit(tx, ty, 1);
    lims = [min(tx), max(tx)];
    lims = [lims(1) - range(lims)*0.05, lims(2) + range(lims)*0.05];
    plot(lims, polyval(p1, lims), 'Color', [.6 .6 .6], 'LineStyle', '--')
    scatter(tx, ty, 50, [.6 .6 .6], 'filled')
    % Formatting
    xlabel(sprintf('%s JND (mm)', t{x}))
    ylabel(sprintf('%s JND (mm)', t{y}))
    [r,p] = corr(tx, ty, 'Type', 'Pearson');
    p = p * 6; % Bonferroni correction
    [tx,ty] = GetAxisPosition(gca,100,100);
    text(tx,ty, sprintf('r = %0.3f\n%s', r, pStr(p)), ...
        'VerticalAlignment', 'top','HorizontalAlignment', 'right')

    xs = xs + x_m + x_w;
end
AddFigureLabels(gcf, [0.05 -0.03])
shg
print(gcf, fullfile(FigurePath, "SuppFig1_2PD.png"), '-dpng', '-r300')

%% Plot single subject psychometric functions
i = 11; % Which subject index to plot
sf = GetSigmoid(2);
clf;
t = tiledlayout('flow');
for j = 1:length(lf)
    nexttile(); hold on
    % Summarize
    tab = ResponseTable_ConditionMean(subjectData(i).(lf{j}).data);
    % Plot
    scatter(tab.distance, tab.response, 40, lf_colors(j,:), 'filled')
    xq = linspace(min(tab.distance), max(tab.distance));
    plot(xq, sf(subjectData(i).(lf{j}).SigCoeffs, xq), 'Color', [.6 .6 .6]);
    % Formatting
    title(lf{j})
    set(gca, 'YLim', [0 1], 'YTick', [0 1])
end
shg

%% Plot all combinations of JND & measurement
clf; i = 1;
for x = 1:width(meas_table)
    for y = 1:width(jnd_table)
        subplot(width(meas_table), width(jnd_table), i); hold on
        % Best fit line and scatter
        nan_idx = isnan(meas_table{:,x}) | isnan(jnd_table{:,y});
        p1 = polyfit(meas_table{~nan_idx,x}, jnd_table{~nan_idx,y}, 1);
        lims = [min(meas_table{~nan_idx,x}), max(meas_table{~nan_idx,x})];
        lims = [lims(1) - range(lims)*0.05, lims(2) + range(lims)*0.05];
        plot(lims, polyval(p1, lims), 'Color', [.6 .6 .6], 'LineStyle', '--')
        scatter(meas_table{~nan_idx,x}, jnd_table{~nan_idx,y}, 50, [.6 .6 .6], 'filled')
        % Formatting
        if y == 1
            ylabel(strrep(meas_table.Properties.VariableNames{x}, '_', ' '), 'FontWeight', 'bold')
        end
        if x == 1
            title(jnd_table.Properties.VariableNames{y})
        end
        % Stats text
        [r,p] = corr(meas_table{~nan_idx,x}, jnd_table{~nan_idx,y}, 'Type', 'Pearson');
        if p < 0.05 % p-value is not corrected for multiple comparisons here
            fw = 'bold';
        else
            fw = 'normal';
        end
        [tx,ty] = GetAxisPosition(gca, 5,95);
        text(tx,ty, sprintf('r = %0.3f\n%s', r, pStr(p)), 'VerticalAlignment', 'top', 'FontWeight', fw)
        i = i + 1;
    end
end
shg

%% Check subject performance consistency
clf; i = 1;
for x = 1:width(jnd_table)
    for y = 1:width(jnd_table)
        if x <= y
            i = i + 1;
            continue
        end
        subplot(width(jnd_table), width(jnd_table), i); hold on
        
        % Best fit line and scatter
        nan_idx = isnan(jnd_table{:,x}) | isnan(jnd_table{:,y});
        p1 = polyfit(jnd_table{~nan_idx,x}, jnd_table{~nan_idx,y}, 1);
        lims = [min(jnd_table{~nan_idx,x}), max(jnd_table{~nan_idx,x})];
        lims = [lims(1) - range(lims)*0.05, lims(2) + range(lims)*0.05];
        plot(lims, polyval(p1, lims), 'Color', [.6 .6 .6], 'LineStyle', '--')
        scatter(jnd_table{~nan_idx,x}, jnd_table{~nan_idx,y}, 50, [.6 .6 .6], 'filled')
        % Formatting
        ylabel(strrep(jnd_table.Properties.VariableNames{y}, '_', ' '))
        xlabel(jnd_table.Properties.VariableNames{x})
        
        % Stats text
        [r,p] = corr(jnd_table{~nan_idx,x}, jnd_table{~nan_idx,y}, 'Type', 'Pearson');
        p = p * 6; % Bonferroni correction
        if p < 0.05
            fw = 'bold';
        else
            fw = 'normal';
        end
        [tx,ty] = GetAxisPosition(gca, 5,95);
        text(tx,ty, sprintf('r = %0.3f\n%s', r, pStr(p)), 'VerticalAlignment', 'top', 'FontWeight', fw)
        i = i + 1;
    end
end
shg
