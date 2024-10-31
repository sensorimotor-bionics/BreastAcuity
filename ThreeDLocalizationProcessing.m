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
                point_grid.y = point_grid.y .* -1; % Fix this reference frame issues
                center_point = mean([point_grid.x(1:2), point_grid.y(1:2), point_grid.z(1:2)],1);
                point_grid{:, 'dist'} = pdist2(center_point, ...
                                               [point_grid.x, point_grid.y, point_grid.z])';
            end
        end

        % Combine data across blocks
        resp_data = cat(1, data_cell{:});
        resp_data.y = resp_data.y .* -1;  % Fix this reference frame issues
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
    [total_error, abs_theta, delta_theta] = deal(cell(height(point_grid), 1));
    [mean_error, bias_error, imprecision_error, abs_bias_theta, delta_bias_theta] = deal(zeros(height(point_grid), 1));
    % The 'center' is the center of the first two points
    center_xyz = mean([point_grid.x(1:2), point_grid.y(1:2), point_grid.z(1:2)], 1);
    for p = 1:height(point_grid)
        % Get matching points
        idx = resp_grid.Point == point_grid.Point(p);
        point_xyz = [point_grid.x(p), point_grid.y(p), point_grid.z(p)];
        resp_xyz = [resp_grid.x(idx), resp_grid.y(idx), resp_grid.z(idx)];
        mean_resp_xyz = mean(resp_xyz, 1); % Mean response for bias metric

        % Compute imprecision
        imprecision_error(p) = mean(pdist2(mean_resp_xyz, resp_xyz, 'euclidean'));
        % Compute error/bias for each
        total_error{p} = pdist2(point_xyz, resp_xyz, 'euclidean');
        mean_error(p) = mean(total_error{p});
        bias_error(p) = pdist2(point_xyz, mean_resp_xyz, 'euclidean');

        % Calculate angle between response vector and center vector
        % dimensions 1,3 are x,y so angle is only computed in 2D
        % Angle between point and center
        center_dxdy = center_xyz([1,2]) - point_xyz([1,2]);
        center_theta = atan2(center_dxdy(1), center_dxdy(2));
        for j = 1:size(resp_xyz, 1)
            % Response versus point
            dxdy = resp_xyz(j, [1,2]) - point_xyz([1,2]);
            abs_theta{p}(j) = atan2(dxdy(1), dxdy(2));
            % Angle between error and center
            delta_theta{p}(j) = center_theta - abs_theta{p}(j);
        end
        % Biased angle
        mean_dxdy = mean_resp_xyz([1,2]) - point_xyz([1,2]);
        abs_bias_theta(p) = atan2(mean_dxdy(1), mean_dxdy(2));
        % Angle between bias error and center
        delta_bias_theta(p) = center_theta - abs_bias_theta(p);
    end
    % Add values to output struct
    point_grid.tot_err = total_error;
    point_grid.err_mean = mean_error;
    point_grid.bias_err = bias_error;
    point_grid.imp_err = imprecision_error;
    point_grid.abs_theta = abs_theta;
    point_grid.abs_theta_bias = abs_bias_theta;
    point_grid.delta_theta = delta_theta;
    point_grid.delta_bias_theta = delta_bias_theta;
end