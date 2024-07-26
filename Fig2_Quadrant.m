%% Load QuadrantData from each participant
data_path = fullfile(DataPath(), 'raw_data', 'quadrant');
flist = dir(fullfile(data_path, '*.mat'));
pooledData = struct();
for f = 1:length(flist)
    temp = load(fullfile(data_path, flist(f).name));
    pooledData(f).Subject = temp.data.subjID;
    pooledData(f).areola = temp.data.areola;
    pooledData(f).nipple = temp.data.nipple;
end
clearvars -except pooledData
nSubjects = length(pooledData);

%% For each participant what is the odds of getting x% correct
% Compute percent correct for each subject
[ar_pc, nip_pc] = deal(zeros(nSubjects, 1));
for i = 1:nSubjects
    ar_pc(i) = sum(pooledData(i).areola.location == pooledData(i).areola.stylusresponse) / ...
        length(pooledData(i).areola.location);
    nip_pc(i) = sum(pooledData(i).nipple.location == pooledData(i).nipple.stylusresponse) / ...
        length(pooledData(i).nipple.location);
end

% Simulate odds assuming random guessing for each participant
num_perms = 1e4; % 10k guesses
meta_null = zeros(nSubjects, num_perms, 2); % Keep track of null for each participant
% Technically the nulls would be the same for areola and nipple but we'll keep them separate for added rigour
Classes = unique(pooledData(1).areola.location); % Classes are the same across participants and were equally delivered
[ar_pc_p, nip_pc_p] = deal(zeros(nSubjects, 1)); % pValue for each participant:location
for i = 1:nSubjects
    % Get number of trials for each participant (not sure why it varies)
    naTrials = length(pooledData(i).areola.location);
    nnTrials = length(pooledData(i).nipple.location);
    % Guess for each trial p times (sample with replacement)
    null = zeros(num_perms, 2); % Do areola and nipple separately
    for p = 1:num_perms
        guesses = datasample(Classes, naTrials);
        null(p,1) = sum(guesses == pooledData(i).areola.location) / naTrials;
        guesses = datasample(Classes, nnTrials);
        null(p,2) = sum(guesses == pooledData(i).nipple.location) / nnTrials;
    end
    % Determine the proportion of times our performance exceeded the null and perform one sided test
    ar_pc_p(i) = 1 - (sum(ar_pc(i) > null(:,1)) / num_perms);
    nip_pc_p(i) = 1 - (sum(nip_pc(i) > null(:,2)) / num_perms);
    % Allocate to meta null for cross-participant analysis
    meta_null(i,:,:) = null;
end
% Bonferroni correction
ar_pc_p = ar_pc_p .* nSubjects;
ar_pc_h = ar_pc_p < 0.05;
nip_pc_p = nip_pc_p .* nSubjects;
nip_pc_h = nip_pc_p < 0.05;

% Average across participants in meta_null to get cross-participant null
meta_null_mean = squeeze(mean(meta_null, 1));
meta_null_std = std(meta_null_mean, 1, 1);
% Compute p-value and standardized effect size (z-score) for nipple and areola
ar_meta_pc_p = 1 - (sum(mean(ar_pc) > meta_null_mean(:,1)) / num_perms);
ar_meta_pc_z = (mean(ar_pc) - mean(meta_null_mean(:,1))) / meta_null_std(1);
nip_meta_pc_p = 1 - (sum(mean(nip_pc) > meta_null_mean(:,2)) / num_perms);
nip_meta_pc_z = (mean(nip_pc) - mean(meta_null_mean(:,2))) / meta_null_std(2);

% Print number of participants who performed above chance
fprintf('%d / %d subjects performed better than chance (areola). Z = %0.3f; p = %0.3f\n', ...
    sum(ar_pc_h), nSubjects, ar_meta_pc_z, ar_meta_pc_p)
fprintf('%d / %d subjects performed better than chance (nipple). Z = %0.3f; p = %0.3f\n', ...
    sum(nip_pc_h), nSubjects, nip_meta_pc_z, nip_meta_pc_p)

%% Load measurement from 2PD analysis
load(fullfile(DataPath(), '2PD_processed'), 'subjectData')
% Cross reference performance of each subject with measurements
meas_table = NaN(length(subjectData), 4);
sd_sl = {subjectData.Subject};
for s = 1:length(subjectData)
    %%%%%%%%% Missing subject ids?
end
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

%% Make figure
clf;
set(gcf, 'Units', 'Inches', 'Position', [30 1 6.45 2.25])
axes('Position', [0.0 0.125 0.325 0.75]); hold on
    % Plot cross
    plot([-1 1], [-1 1], 'Color', [.6 .6 .6])
    plot([-1 1], [1 -1], 'Color', [.6 .6 .6])
    % Plot circles
    x = linspace(0, 2*pi, 500);
    r = 0.4;
    plot(sin(x) .* r, cos(x) .* r, 'color', [.6 .6 .6])
    r = 1;
    plot(sin(x) .* r, cos(x) .* r, 'color', [.6 .6 .6])
    % Scatter points
    x = linspace(0, 1.5*pi, 4);
    r = 0.2; % Half radius of inner ring
    scatter(sin(x) .* r, cos(x) .* r, 30, [0.26 0.63 0.28], 'filled')
    r = 0.7; % Half way between inner and outer ring
    scatter(sin(x) .* r, cos(x) .* r, 30, [0.26 0.28 0.63], 'filled')
    % Add text
    text(1.25, 0, 'Medial', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', -90)
    text(-1.25, 0, 'Lateral', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90)
    text(0, 1.25, 'Superior', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    text(0, -1.25, 'Inferior', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')

    set(gca, 'DataAspectRatio', [1 1 1], ...
             'XLim', [-2 2], ...
             'XColor', 'none', ...
             'YColor', 'none')

axes('Position', [0.425 0.15 0.25 0.75]); hold on
    x = [0.5 0.5 1.5 1.5];
    y = [0 1 1 0];
    ww = 0.25;

    % Areola bar
    am = mean(ar_pc);
    as = std(ar_pc);
    patch(x, y .* am, [.8 .8 .8], 'EdgeColor','none')
    plot(x, y .* am, 'Color', 'k')
    wc = 1;
    plot([wc-ww, wc+ww, wc, wc, wc-ww, wc+ww], [am-as, am-as, am-as, am+as, am+as, am+as], 'Color', 'k')

    % Nipple bar
    nm = mean(nip_pc);
    ns = std(nip_pc);
    patch(x + 1.5, y .* nm, [.8 .8 .8], 'EdgeColor','none')
    plot(x + 1.5, y .* nm, 'Color', 'k')
    wc = 2.5;
    plot([wc-ww, wc+ww, wc, wc, wc-ww, wc+ww], [nm-ns, nm-ns, nm-ns, nm+ns, nm+ns, nm+ns], 'Color', 'k')

    % Lines between points
    x = repmat([1, 2.5, NaN], [nSubjects, 1])';
    y = [ar_pc, nip_pc, NaN(nSubjects, 1)]';
    plot(x(:), y(:), 'Color', [.4 .4 .4], 'LineWidth', 0.5)
    % Individual points (filled = significant)
    scatter(ones(sum(ar_pc_h),1), ar_pc(ar_pc_h), 50, [.4 .4 .4], 'MarkerFaceColor', [.4 .4 .4])
    scatter(ones(sum(~ar_pc_h),1), ar_pc(~ar_pc_h), 50, [.4 .4 .4])
    scatter(ones(sum(nip_pc_h),1) .* 2.5, nip_pc(nip_pc_h), 50, [.4 .4 .4], 'filled', 'MarkerFaceColor', [.4 .4 .4])
    scatter(ones(sum(~nip_pc_h),1) .* 2.5, nip_pc(~nip_pc_h), 50, [.4 .4 .4])

    % Plot chance
    plot([0 3.5], [0.25 0.25], 'Color', [.4 .4 .4], 'LineStyle', '--')


    set(gca, 'XLim', [0 3.5], ...
             'YLim', [0 1], ...
             'YTick', [0:.25:1], ...
             'XTick', [1, 2.5], ...
             'XTickLabel', {'Areola', 'Nipple'})
    ylabel('p(correct)')
shg