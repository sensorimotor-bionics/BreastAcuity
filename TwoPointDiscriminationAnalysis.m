%% Load 2PD data
xlsx_path = fullfile(DataPath(), 'raw_data', '2PD');
lf = {'hand', 'back', 'breast', 'NAC'};
lf_colors = [rgb(255, 160, 0); rgb(194, 24, 91); rgb(25, 118, 210); rgb(123, 31, 162)];
num_blocks = 3;

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
subjectMeta = readtable(fullfile(DataPath(), 'raw_data', 'SubjectMeta.xlsx'));
sl = [data.Subject];
ll = {data.Location};
us = unique(sl);

subjectData = struct();
for s = 1:length(us)
    subjectData(s).Subject = us{s};
    % Add responses to appropriate substructure
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
p_threshold = 0.6; % Assume participants performing worse than 60% couldn't or weren't doing the task
[jnd_table, meas_table, p_correct_table] = deal(nan(length(subjectData), length(lf)));
for s = 1:length(subjectData)
    for l = 1:length(lf)
        % Get mean response for each distance
        summary_table = ResponseTable_ConditionMean(subjectData(s).(lf{l}).data);
        if any(isnan(summary_table.distance))
            summary_table = summary_table(~isnan(summary_table.distance), :);
        end
        % Fit mean points to sigmoid (n is same across points so no need for weighting)
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

clearvars -except meas_table jnd_table subjectData lf*
save(fullfile(DataPath(), '2PD_processed'))