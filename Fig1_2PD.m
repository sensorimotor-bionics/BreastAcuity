%% Load 2PD data
tld = 'U:\UserFolders\CharlesGreenspon\BB_Acuity';
xlsx_path = fullfile(tld, 'raw_data', '2PD');
lf = {'hand', 'back', 'breast', 'NAC'};
lf_colors = [rgb(255, 160, 0); rgb(194, 24, 91); rgb(25, 118, 210); rgb(123, 31, 162)];
num_blocks = 3;
breast_image = imread(fullfile(DataPath(), '..', 'BreastDrawings.png'));

data = struct(); di = 1;
for s = 1:length(lf)
    flist = dir(fullfile(xlsx_path, lf{s}, '*xlsx'));
    for f = 1:length(flist)
        % Parse fname
        fsplit = split(flist(f).name(1:end-5), '-');
        % Load data from each block
        block_data = cell(num_blocks, 1);
        for b = 1:num_blocks
            try % Some subjects only have 2 blocks so this will error if the sheet doesn't exist
                block_data{b} = readtable(fullfile(xlsx_path, lf{s}, flist(f).name),...
                    'Sheet', sprintf('Block %d', b));                
            end
        end
        % Assign to data struct
        data(di).Subject = join(fsplit(1:2), '-');
        data(di).Location = lf{s};
        data(di).responses = cat(1, block_data{:});
        di = di + 1;
    end
end

%% Make subject data structure & load measurements
subjectMeta = readtable(fullfile(tld, 'raw_data', 'SubjectMeta.xlsx'));
sl = [data.Subject];
ll = {data.Location};
us = unique(sl);

subjectData = struct();
for s = 1:length(us)
    subjectData(s).Subject = us{s};
    % Add responses to appropriat substructure
    for l = 1:length(lf)
        idx = strcmp(sl, us(s)) & strcmp(ll, lf{l});
        subjectData(s).(lf{l}).data = cat(1, data(idx).responses);
    end
    % Add measurements field
    meta_idx = strcmp(subjectMeta.Subject, us{s});
    if sum(meta_idx) ~= 1
        [subjectData(s).measurements.areola, subjectData(s).measurements.nipple,...
         subjectData(s).measurements.bust, subjectData(s).measurements.underbust] = deal(nan);
        continue
    end
    % Areola and nipple are ranges instead of scalars and thus must be handled differently
    if strcmp(subjectMeta.areola(meta_idx), '-')
        subjectData(s).measurements.areola = NaN;
    else
        temp = char(subjectMeta.areola(meta_idx));
        subjectData(s).measurements.areola = range(str2double(split(temp(2:end-1), ',')));
    end
    if strcmp(subjectMeta.nipple(meta_idx), '-')
        subjectData(s).measurements.nipple = NaN;
    else
        temp = char(subjectMeta.nipple(meta_idx));
        subjectData(s).measurements.nipple = range(str2double(split(temp(2:end-1), ',')));
    end
    subjectData(s).measurements.bust = subjectMeta.bust(meta_idx);
    subjectData(s).measurements.underbust = subjectMeta.underbust(meta_idx);
end

%% Analyze each region
p_threshold = 0.6;
[jnd_table, meas_table, p_correct_table] = deal(nan(length(subjectData), length(lf)));
for s = 1:length(subjectData)
    for l = 1:length(lf)
        % Get mean response for each distance
        summary_table = ResponseTable_ConditionMean(subjectData(s).(lf{l}).data);
        if any(isnan(summary_table.distance))
            summary_table = summary_table(~isnan(summary_table.distance), :);
        end
        [~, coeffs, ~, ~, jnd, warn] = FitSigmoid(summary_table.distance, summary_table.response,...
            'NumCoeffs', 2, 'EnableBackup', false);
        subjectData(s).(lf{l}).SigCoeffs = coeffs;
        subjectData(s).(lf{l}).JND = jnd;
        if warn
            fprintf('Subject %d, location %s\n', s, lf{l})
        end

        % Compute percent correct
        p = (sum(subjectData(s).(lf{l}).data.distance > 0 & subjectData(s).(lf{l}).data.response == 1) + ...
             sum(subjectData(s).(lf{l}).data.distance < 0 & subjectData(s).(lf{l}).data.response == 0)) / ...
             height(subjectData(s).(lf{l}).data);
        [subjectData(s).(lf{l}).p_correct, p_correct_table(s,l)] = deal(p);
        % Assign to tables for easy comparison
        if p > p_threshold
            jnd_table(s,l) = jnd;
            meas_table(s,:) = [subjectData(s).measurements.areola, subjectData(s).measurements.nipple,...
                           subjectData(s).measurements.bust, subjectData(s).measurements.underbust];
        end
    end
end

jnd_table = array2table(jnd_table, 'VariableNames', lf);
meas_table = array2table(meas_table, 'VariableNames', {'areola', 'nipple', 'bust', 'underbust'});
meas_table.delta_bust = meas_table.bust - meas_table.underbust;

%% Make figure
sf = GetSigmoid(2);
SetFont('Arial', 9)
clf;
set(gcf, 'Units', 'Inches', 'Position', [30 1 6.45 4.5])
xl = [-30 30];
xq = linspace(xl(1), xl(2));

% Examples - Hand
yp1 = 0.625; yp2 = 0.475; yp3 = 0.09;
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
for j = 1:4
    axes('Position', [xp(j) yp2 .15 .25]); hold on
    tab = ResponseTable_ConditionMean(subjectData(i).(lf{j}).data);
    plot(xq, sf(subjectData(i).(lf{j}).SigCoeffs, xq), 'Color', [.6 .6 .6], 'LineWidth', 2);
    scatter(tab.distance, tab.response, 40, lf_colors(j,:), 'MarkerFaceColor', lf_colors(j,:))
    % Formatting
    set(gca, 'YLim', [0 1], 'YTick', [0 1], 'XLim', xl)
    [tx,ty] = GetAxisPosition(gca, 5,95);
    text(tx,ty, sprintf('JND =\n%0.1f mm', jnd_table{i, lf{j}}), ...
        'VerticalAlignment', 'top', 'Color', 'k')
    if j == 1
        ylabel('p(above)')
    end
end

% JND distributions
axes('Position', [0.1 yp3 .2 .25]); hold on
    for l = 1:length(lf)
        Swarm(l, jnd_table{:,l}, lf_colors(l,:), 'DS', 'Box', 'SWR', 0.5, 'DW', 0.3);
    end
    set(gca, 'XLim', [0.5, 4.5], ...
             'YLim', [0 20], ...
             'XTick', [1:4], ...
             'XTickLabel', {'Hand', 'Back', 'Breast', 'NAC'}, ...
             'YTick', [0:10:20])
    ylabel('JND (mm)')

% Breast scaling
axes('Position', [0.4 yp3 .2 .25]); hold on
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
axes('Position', [0.7 yp3 .2 .25]); hold on
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