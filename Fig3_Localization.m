%%
load(fullfile(DataPath(), 'td_localization_data'))
model = readObj("U:\UserFolders\CharlesGreenspon\BB_Acuity\ModelBreast_Subj001.obj");

%% Compute error for each trial and angle of error relative to breast + means across blocks
clf; 
subplot(1,2,1); hold on
for p = 1:length(data)
    x = data(p).breast_grid.dist;
    y = data(p).breast_grid.mean_err;
    scatter(rescale(x),rescale(y))
end

subplot(1,2,2); hold on
for p = 1:length(data)
    x = data(p).back_grid.dist;
    y = data(p).back_grid.mean_err;
    scatter(rescale(x),rescale(y))
end


shg

%% Compute mean error for back and breast
mean_loc_error = zeros(length(td_localization_data), 2);
for i = 1:length(td_localization_data)
    mean_loc_error(i,1) = mean(td_localization_data(i).breast_grid.err_mean);
    mean_loc_error(i,2) = mean(td_localization_data(i).back_grid.err_mean);
end
[~, mean_err_p, ~, mean_err_t] = ttest(mean_loc_error(:,1), mean_loc_error(:,2));

%% Plot
imp_col = rgb(186, 104, 200);
bias_col = rgb(100, 181, 246);
err_col = rgb(77, 182, 172);

clf; 
set(gcf, 'Units', 'Inches', 'Position', [30 1 6.45 4.25])
axes('Position', [0.05 0.55 0.3 0.35]); hold on
    trisurf(model.f.v, model.v(:,1), model.v(:,2), model.v(:,3), 'EdgeColor', 'none', 'FaceColor', [.6 .6 .6])
    i = 1;
    point_grid = td_localization_data(i).breast_grid;
    resp_grid = td_localization_data(i).breast_resp;
    for p = 1:height(point_grid)
        % Get matching points
        idx = resp_grid.Point == point_grid.Point(p);
        point_xyz = [point_grid.x(p), point_grid.z(p), point_grid.y(p).*-1];
        resp_xyz = [resp_grid.x(idx), resp_grid.z(idx), resp_grid.y(idx).*-1];
        resp_mean = mean(resp_xyz, 1); % Get mean response
        b = point_xyz - center_xyz;
        for j = 1:size(resp_xyz, 1) % Plot individual points
            plot3([point_xyz(1), resp_xyz(j,1)], [point_xyz(2), resp_xyz(j,2)], [point_xyz(3), resp_xyz(j,3)], 'Color', imp_col)
        end
        % Plot mean
        plot3([point_xyz(1), resp_mean(1)], [point_xyz(2), resp_mean(2)], [point_xyz(3), resp_mean(3)], 'Color', bias_col)
    end
    % Plot target
    scatter3(point_grid.x, ...
             point_grid.z, ...
             point_grid.y*-1, 30, 'k', 'Marker', 'x', 'LineWidth', 2)
    set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none', 'View', [0, 90], 'DataAspectRatio', [1 1 1])
    set(gca, 'DataAspectRatio', [1 1 1], ...
             'View', [160, 5], ...
             'XColor', 'none', ...
             'YColor', 'none', ...
             'ZColor', 'none',...
             'XLim', [-150 150], ...
             'YLim', [50, 150], ...
             'ZLim', [-150, 125])
    camlight(30, 0)

% axes('Position', [0.4 0.55 0.25 0.35]); hold on
% 
% 
% axes('Position', [0.7 0.55 0.25 0.35]); hold on
%     i = 1;
%     point_grid = td_localization_data(i).breast_grid;
%     resp_grid = td_localization_data(i).breast_resp;
%     center_xyz = mean([point_grid.x(1:2), point_grid.y(1:2), point_grid.z(1:2)], 1);
%     for p = 1:height(point_grid)
%         % Get matching points
%         idx = resp_grid.Point == point_grid.Point(p);
%         point_xyz = [point_grid.x(p), point_grid.y(p), point_grid.z(p)];
%         resp_xyz = [resp_grid.x(idx), resp_grid.y(idx), resp_grid.z(idx)];
%         resp_mean = mean(resp_xyz, 1); % Get mean response
%         b = point_xyz - center_xyz;
%         for j = 1:size(resp_xyz, 1) % Plot individual points
%             plot3([point_xyz(1), resp_xyz(j,1)], [point_xyz(2), resp_xyz(j,2)], [point_xyz(3), resp_xyz(j,3)], 'Color', imp_col)
%         end
%         % Plot mean
%         plot3([point_xyz(1), resp_mean(1)], [point_xyz(2), resp_mean(2)], [point_xyz(3), resp_mean(3)], 'Color', bias_col)
%     end
%     % Plot target
%     scatter3(point_grid.x, point_grid.y, point_grid.z, 'MarkerEdgeColor', 'k', 'Marker', 'x', 'LineWidth', 2)
%     set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none', 'View', [0, 45], 'DataAspectRatio', [1 1 1])


axes('Position', [0.05 0.125 0.3 0.35]); hold on
    lims = [0 80];
    plot(lims, lims, 'Color', [.6 .6 .6], 'LineStyle', '--')
    scatter(mean_loc_error(:,1), mean_loc_error(:,2), 30, 'MarkerEdgeColor', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6])
    % Formatting
    set(gca, 'XLim', lims, ...
             'YLim', lims, ...
             'XTick', [0 40 80], ...
             'YTick', [0 40 80], ...
             'DataAspectRatio', [1 1 1], ...
             'XTickLabelRotation', 0)
    xlabel('Breast Error (mm)')
    ylabel('Back Error (mm)')
    % Stats
    [x,y] = GetAxisPosition(gca, 95, 5);
    text(x,y, sprintf('%s', pStr(mean_err_p)), 'HorizontalAlignment', 'right', 'VerticalAlignment','bottom')
shg



%% Quick plot of individual
% i = 5;
% c = turbo(height(data(i).back_grid));
% clf; 
% subplot(1,2,1); hold on
% scatter3(data(i).back_grid.x(1), data(i).back_grid.y(1), data(i).back_grid.z(1), 100, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
% for j = 1:height(data(i).back_grid)
%     scatter3(data(i).back_grid.x(j), data(i).back_grid.y(j), data(i).back_grid.z(j), 'MarkerEdgeColor', c(j,:), 'MarkerFaceColor', c(j,:))
%     idx = data(i).back_resp.Point == data(i).back_grid.Point(j);
%     scatter3(data(i).back_resp.x(idx), data(i).back_resp.y(idx), data(i).back_resp.z(idx), 'Marker','x', 'MarkerEdgeColor', c(j,:))
% end
% 
% subplot(1,2,2); hold on
% scatter3(data(i).breast_grid.x(1), data(i).breast_grid.y(1), data(i).breast_grid.z(1), 100, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
% for j = 1:height(data(i).breast_grid)
%     scatter3(data(i).breast_grid.x(j), data(i).breast_grid.y(j), data(i).breast_grid.z(j), 'MarkerEdgeColor', c(j,:), 'MarkerFaceColor', c(j,:))
%     idx = data(i).breast_resp.Point == data(i).breast_grid.Point(j);
%     scatter3(data(i).breast_resp.x(idx), data(i).breast_resp.y(idx), data(i).breast_resp.z(idx), 'Marker','x', 'MarkerEdgeColor', c(j,:))
% end
% 
% 
% i = 8;
% point_grid = data(i).breast_grid;
% center_xyz = [point_grid.x(1), point_grid.y(1), point_grid.z(1)];
% resp_grid = data(i).breast_resp;
% c = winter(256);
% cc = linspace(0, pi/2, 256);
% 
% clf; hold on
% scatter3(point_grid.x(1), point_grid.y(1), point_grid.z(1), 100, 'k', 'filled')
% scatter3(point_grid.x, point_grid.y, point_grid.z, 'k', 'filled')
% for p = 1:height(point_grid)
%     % Get matching points
%     idx = resp_grid.Point == point_grid.Point(p);
%     point_xyz = [point_grid.x(p), point_grid.y(p), point_grid.z(p)];
%     resp_xyz = [resp_grid.x(idx), resp_grid.y(idx), resp_grid.z(idx)];
%     % Compute error for each
%     euc_err{p} = pdist2(point_xyz, resp_xyz, 'euclidean');
%     % Calculate angle between response vector and nipple vector
%     b = point_xyz - center_xyz;
%     for j = 1:size(resp_xyz, 1)
%         a = point_xyz - resp_xyz(j,:);
%         theta_center{p}(j) = subspace(a',b');
%         [~, cidx] = min(abs(cc - theta_center{p}(j)));
%         scatter3(resp_xyz(j,1), resp_xyz(j,2), resp_xyz(j,3), 'xk')
%         plot3([point_xyz(1), resp_xyz(j,1)], [point_xyz(2), resp_xyz(j,2)], [point_xyz(3), resp_xyz(j,3)], 'Color', c(cidx,:))
%     end
% end