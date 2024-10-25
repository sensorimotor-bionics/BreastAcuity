%% Load QuadrantData from each participant
data_path = fullfile(DataPath(), 'raw_data', '3D', 'resp');
flist = dir(fullfile(data_path, '*.xlsx'));

% Get list of participants
pid = {flist.name};
for p = 1:length(pid)
    psplit = strsplit(pid{p}, '_');
    pid{p} = psplit{1};
end
pid = unique(pid);

% For each participant load data into struct
data = struct();
for p = 1:length(pid)
    data(p).SubjectID = pid{p};
    % Filter files by SubjectID
    f_filt = flist(contains({flist.name}, pid{p}));
    for f = 1:length(f_filt)
        [point_data, resp_data] = deal([]);
        ffpath = fullfile(f_filt(f).folder, f_filt(f).name);
        % Assume the sheets are X, Block 1, 2, 3
        sn = sheetnames(ffpath);
        data_cell = cell(length(sn)-1, 1); si = 1;
        for s = 1:length(sn)
            % Load the file
            fdata = readtable(ffpath, "Sheet", sn{s}, 'ReadRowNames', false);
            if contains(sn{s}, 'block', 'IgnoreCase', true)
                fdata{:, 'block'} = str2double(sn{s}(end)); % Assign block to table
                data_cell{si} = fdata; % Add tables to cell for later
                si = si + 1;
            else % Assume it's the coordinate data
                point_grid = fdata;
                center_point = mean([point_grid.x(1:2), point_grid.y(1:2), point_grid.z(1:2)],1);
                point_grid{:, 'dist'} = pdist2(center_point, ...
                                               [point_grid.x, point_grid.y, point_grid.z])';
            end
        end

        % Combine data across blocks
        resp_data = cat(1, data_cell{:});
        % Compute errors
        point_grid = compute_resp_error(point_grid, resp_data);

        % Asign based on location
        if contains(f_filt(f).name, 'back', 'IgnoreCase', true)
            data(p).back_grid = point_grid;
            data(p).back_resp = resp_data;
        elseif contains(f_filt(f).name, 'breast', 'IgnoreCase', true)
            data(p).breast_grid = point_grid;
            data(p).breast_resp = resp_data;
        else
            error('invalid location')
        end
    end
end

td_localization_data = data;
clearvars -except td_localization_data
save(fullfile(DataPath(), 'td_localization_data'))


function point_grid = compute_resp_error(point_grid, resp_grid)
    % Setup output
    [euc_err, theta_center] = deal(cell(height(point_grid), 1));
    [bias, mean_error, mean_theta, imprecision] = deal(zeros(height(point_grid), 1));
    % The 'center' is the center of the first two points
    center_xyz = mean([point_grid.x(1:2), point_grid.y(1:2), point_grid.z(1:2)], 1);
    for p = 1:height(point_grid)
        % Get matching points
        idx = resp_grid.Point == point_grid.Point(p);
        point_xyz = [point_grid.x(p), point_grid.y(p), point_grid.z(p)];
        resp_xyz = [resp_grid.x(idx), resp_grid.y(idx), resp_grid.z(idx)];
        mean_resp_xyz = mean(resp_xyz, 1);
        % Compute imprecision
        imprecision(p) = mean(pdist2(mean_resp_xyz, resp_xyz, 'euclidean'));
        % Compute error/bias for each
        euc_err{p} = pdist2(point_xyz, resp_xyz, 'euclidean');
        mean_error(p) = mean(euc_err{p});
        bias(p) = pdist2(point_xyz, mean_resp_xyz, 'euclidean');
        % Calculate angle between response vector and nipple vector
        b = point_xyz - center_xyz;
        for j = 1:size(resp_xyz, 1)
            a = point_xyz - resp_xyz(j,:);
            theta_center{p}(j) = subspace(a',b');
        end
        a = point_xyz - mean_resp_xyz;
        mean_theta(p) = subspace(a',b');
    end
    % Add values to output struct
    point_grid.err = euc_err;
    point_grid.err_mean = mean_error;
    point_grid.bias = bias;
    point_grid.imprescision = imprecision;
    point_grid.theta = theta_center;
    point_grid.theta_mean = mean_theta;
end