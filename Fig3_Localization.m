%% Load raw data
load(fullfile(DataPath(), 'td_localization_data'))
num_participants = length(td_localization_data);
breast_model = readObj(fullfile(DataPath(), 'ModelBreast_Subj001.obj'));
back_model = readObj(fullfile(DataPath(), 'ModelBack_Subj001.obj'));

% Colors
imp_col = rgb(186, 104, 200);
bias_col = rgb(100, 181, 246);
err_col = rgb(77, 182, 172);
center_col = rgb(233, 30, 99);

%% Breast vs back error
mean_loc_error = zeros(num_participants, 2);
for i = 1:num_participants
    mean_loc_error(i,1) = mean(td_localization_data(i).breast_grid.err_mean);
    mean_loc_error(i,2) = mean(td_localization_data(i).back_grid.err_mean);
end
[~, mean_err_p, ~, mean_err_t] = ttest(mean_loc_error(:,1), mean_loc_error(:,2));

%% Error vs Bias vs Imprecision
% First discretize and average for each participant
xe = [0, 5, 25, 40, 60, 90];
xc  = xe(1:end-1) + diff(xe)./2;
[ebi_breast, ebi_back] = deal(NaN(num_participants, length(xc), 3)); % Val for each participant (error, bias, imprecision)
for p = 1:num_participants % Participant
    for b = 1:length(xc) % Bin
        % Breast
        idx = td_localization_data(p).breast_grid.dist >= xe(b) & td_localization_data(p).breast_grid.dist < xe(b+1);
        ebi_breast(p,b,1) = mean(td_localization_data(p).breast_grid.err_mean(idx));
        ebi_breast(p,b,2) = mean(td_localization_data(p).breast_grid.bias_err(idx));
        ebi_breast(p,b,3) = mean(td_localization_data(p).breast_grid.imp_err(idx));
        % Back
        idx = td_localization_data(p).back_grid.dist >= xe(b) & td_localization_data(p).back_grid.dist < xe(b+1);
        ebi_back(p,b,1) = mean(td_localization_data(p).back_grid.err_mean(idx));
        ebi_back(p,b,2) = mean(td_localization_data(p).back_grid.bias_err(idx));
        ebi_back(p,b,3) = mean(td_localization_data(p).back_grid.imp_err(idx));
    end
end


% 3-way ANOVA (participant, distance from center)
[P,B,L] = meshgrid([1:num_participants], [1:length(xc)], [1:3]);
a_breast = anovan(ebi_breast(:), {P(:), B(:), L(:)}, 'varnames', {'Participant', 'Distance', 'Type'});
a_back = anovan(ebi_back(:), {P(:), B(:), L(:)}, 'varnames', {'Participant', 'Distance', 'Type'});


%% Theta dist
be = linspace(-pi,pi,12);
bc = be(1:end-1) + diff(be)./2;
[abs_theta_dist, delta_theta_dist] = deal(cell(num_participants, 2));
for p = 1:num_participants
    abs_theta_dist{p,1} = histcounts(cat(2, td_localization_data(p).breast_grid.abs_theta{:}), be, 'Normalization', 'Probability');
    abs_theta_dist{p,2} = histcounts(cat(2, td_localization_data(p).back_grid.abs_theta{:}), be, 'Normalization', 'Probability');
    delta_theta_dist{p,1} = histcounts(cat(2, td_localization_data(p).breast_grid.delta_theta{:}), be, 'Normalization', 'Probability');
    delta_theta_dist{p,2} = histcounts(cat(2, td_localization_data(p).back_grid.delta_theta{:}), be, 'Normalization', 'Probability');
end


%% Plot
SetFont('Arial', 9)
clf; 
i = 1; % Participant that consented to having model shown

set(gcf, 'Units', 'Inches', 'Position', [30 1 4.5 6.5])
% Breast example
axes('Position', [0.05 0.625 0.45 0.35]); hold on
    trisurf(breast_model.f.v, breast_model.v(:,1), breast_model.v(:,2), breast_model.v(:,3), 'EdgeColor', 'none', 'FaceColor', [.6 .6 .6])
    point_grid = td_localization_data(i).breast_grid;
    resp_grid = td_localization_data(i).breast_resp;
    center_xyz = mean([point_grid.x(1:2), point_grid.y(1:2), point_grid.z(1:2)], 1);
    for p = 1:height(point_grid)
        % Get matching points
        idx = resp_grid.Point == point_grid.Point(p);
        point_xyz = [point_grid.x(p), point_grid.z(p), point_grid.y(p)];
        resp_xyz = [resp_grid.x(idx), resp_grid.z(idx), resp_grid.y(idx)];
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
             point_grid.y, 30, 'k', 'Marker', 'x', 'LineWidth', 2)
    set(gca, 'DataAspectRatio', [1 1 1], ...
             'View', [160, 5], ...
             'XColor', 'none', ...
             'YColor', 'none', ...
             'ZColor', 'none',...
             'XLim', [-150 150], ...
             'YLim', [50, 150], ...
             'ZLim', [-150, 125])
    camlight(30, 0)

% Back example
axes('Position', [0.55 0.625 0.375 0.35]); hold on
    trisurf(back_model.f.v, back_model.v(:,1), back_model.v(:,2), back_model.v(:,3), 'EdgeColor', 'none', 'FaceColor', [.6 .6 .6])
    point_grid = td_localization_data(i).back_grid;
    resp_grid = td_localization_data(i).back_resp;
    center_xyz = mean([point_grid.x(1:2), point_grid.y(1:2), point_grid.z(1:2)], 1);
    for p = 1:height(point_grid)
        % Get matching points
        idx = resp_grid.Point == point_grid.Point(p);
        point_xyz = [point_grid.x(p), point_grid.z(p), point_grid.y(p)];
        resp_xyz = [resp_grid.x(idx), resp_grid.z(idx), resp_grid.y(idx)];
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
             point_grid.y, 30, 'k', 'Marker', 'x', 'LineWidth', 2)
    set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none', 'View', [0, 90], 'DataAspectRatio', [1 1 1])
    set(gca, 'DataAspectRatio', [1 1 1], ...
             'View', [-175, 1], ...
             'XColor', 'none', ...
             'YColor', 'none', ...
             'ZColor', 'none',...
             'XLim', [-150 150], ...
             'YLim', [50, 200], ...
             'ZLim', [-200, 150])
    camlight(30, 0)


% Breast error
axes('Position', [0.1 0.4125 0.35 0.2]); hold on
    AlphaLine(xc, ebi_breast(:,:,1), err_col, 'LineWidth', 1.5)
    AlphaLine(xc, ebi_breast(:,:,2), bias_col, 'LineWidth', 1.5)
    AlphaLine(xc, ebi_breast(:,:,3), imp_col, 'LineWidth', 1.5)
    set(gca, 'XLim', [0, 80], ...
             'XTick', [0:20:80], ...
             'YLim', [0 80], ...
             'YTick', [0:20:80])
    xlabel('Distance from Nipple (mm)')
    ylabel('Error (mm)')
    [x,y] = GetAxisPosition(gca, 95, 95);
    text(x,y, ColorText({'Error', 'Bias', 'Imprecision'}, [err_col; bias_col; imp_col]), ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')

% Back error
axes('Position', [0.55 0.4125 0.35 0.2]); hold on
    AlphaLine(xc, ebi_back(:,:,1), err_col, 'LineWidth', 1.5)
    AlphaLine(xc, ebi_back(:,:,2), bias_col, 'LineWidth', 1.5)
    AlphaLine(xc, ebi_back(:,:,3), imp_col, 'LineWidth', 1.5)
    set(gca, 'XLim', [0, 80], ...
             'XTick', [0:20:80], ...
             'YLim', [0 80], ...
             'YTick', [0:20:80])
    xlabel('Distance from Scapula (mm)')

% Breast polar plot 
axes('Position', [0.075 0.075 0.35 0.2], 'XColor', 'none', 'YColor', 'none'); hold on
    % Background
    rho_max = 0.2;
    rho_ticks = 0.2;
    rho_labels = 0.2;
    theta_ticks = deg2rad([0:45:315]);
    theta_labels = deg2rad([0:90:270]);

    set(gca, 'XLim', [-rho_max rho_max],...
             'YLim', [-rho_max rho_max],...
             'DataAspectRatio', [1 1 1])
    
    % Theta ticks
    [x, y] = pol2cart(theta_ticks, repmat(rho_max*1.1, [1, length(theta_ticks)]));
    for i = 1:length(theta_ticks)
       plot([0,x(i)], [0, y(i)], 'Color', [0.8 0.8 0.8], 'LineWidth', 0.1)
    end
    
    % Theta labels
    if isempty(theta_labels) == 0
        [x, y] = pol2cart(theta_labels, repmat(rho_max*1.25, [1, length(theta_labels)]));
        for t = 1:length(theta_labels)
            text(x(t),y(t), [num2str(rad2deg(theta_labels(t))), char(176)],...
                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center')
        end
    end
    
    % Rho ticks
    t = linspace(0, 2*pi, 1000);
    for r = rho_ticks
        [x, y] = pol2cart(t, repmat(r, [1, 1000]));
        plot(x, y, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.1)
    end
    
    % Rho labels
    rt = theta_ticks(2) / 2;
    ro = rho_max * 0.05;
    if isempty(rho_labels) == 0
        [x, y] = pol2cart(repmat(rt, [1,length(rho_labels)]), rho_labels+ro);
        for t = 1:length(rho_labels)
            text(x(t),y(t), num2str(rho_labels(t)),...
                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left')
        end
    end
    
    % Plot error data
    theta = [bc, bc(1)];
    rho = cat(1, abs_theta_dist{:,1});
    rho_mean = mean(rho,1); rho_mean = [rho_mean, rho_mean(1)];
    rho_std = std(rho, [], 1); rho_std = [rho_std, rho_std(1)];
    num_lines = size(theta,1);
    for l = 1:num_lines
        rho_outer = rho_mean(l,:) + rho_std(l,:);
        rho_inner = rho_mean(l,:) - rho_std(l,:);
        
        [x, y] = pol2cart(theta(l,:), rho_outer);
        [x2, y2] = pol2cart(theta(l,:), rho_inner);
        fill([x, fliplr(x2)],[y, fliplr(y2)], [.6 .6 .6],...
            'EdgeColor', [.6 .6 .6],'FaceAlpha', 0.1, 'EdgeAlpha', 0)
        [x3, y3] = pol2cart(theta(l,:), rho_mean(l,:));
        plot(x3,y3, 'Color', [.6 .6 .6], 'LineWidth', 1)
    end

    % Plot bias data
    theta = [bc, bc(1)];
    rho = cat(1, delta_theta_dist{:,1});
    rho_mean = mean(rho,1); rho_mean = [rho_mean, rho_mean(1)];
    rho_std = std(rho, [], 1); rho_std = [rho_std, rho_std(1)];
    num_lines = size(theta,1);
    for l = 1:num_lines
        rho_outer = rho_mean(l,:) + rho_std(l,:);
        rho_inner = rho_mean(l,:) - rho_std(l,:);
        
        [x, y] = pol2cart(theta(l,:), rho_outer);
        [x2, y2] = pol2cart(theta(l,:), rho_inner);
        fill([x, fliplr(x2)],[y, fliplr(y2)], center_col,...
            'EdgeColor', center_col,'FaceAlpha', 0.1, 'EdgeAlpha', 0)
        [x3, y3] = pol2cart(theta(l,:), rho_mean(l,:));
        plot(x3,y3, 'Color', center_col, 'LineWidth', 1)
    end


% Back polar plot 
axes('Position', [0.575 0.075 0.3 0.2], 'XColor', 'none', 'YColor', 'none'); hold on
    % Background
    rho_max = 0.2;
    rho_ticks = 0.2;
    rho_labels = 0.2;
    theta_ticks = deg2rad([0:45:315]);
    theta_labels = deg2rad([0:90:270]);

    set(gca, 'XLim', [-rho_max rho_max],...
             'YLim', [-rho_max rho_max],...
             'DataAspectRatio', [1 1 1])
    
    % Theta ticks
    [x, y] = pol2cart(theta_ticks, repmat(rho_max*1.1, [1, length(theta_ticks)]));
    for i = 1:length(theta_ticks)
       plot([0,x(i)], [0, y(i)], 'Color', [0.8 0.8 0.8], 'LineWidth', 0.1)
    end
    
    % Theta labels
    if isempty(theta_labels) == 0
        [x, y] = pol2cart(theta_labels, repmat(rho_max*1.25, [1, length(theta_labels)]));
        for t = 1:length(theta_labels)
            text(x(t),y(t), [num2str(rad2deg(theta_labels(t))), char(176)],...
                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center')
        end
    end
    
    % Rho ticks
    t = linspace(0, 2*pi, 1000);
    for r = rho_ticks
        [x, y] = pol2cart(t, repmat(r, [1, 1000]));
        plot(x, y, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.1)
    end
    
    % Rho labels
    rt = theta_ticks(2) / 2;
    ro = rho_max * 0.05;
    if isempty(rho_labels) == 0
        [x, y] = pol2cart(repmat(rt, [1,length(rho_labels)]), rho_labels+ro);
        for t = 1:length(rho_labels)
            text(x(t),y(t), num2str(rho_labels(t)),...
                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left')
        end
    end
    
    % Plot error data
    theta = [bc, bc(1)];
    rho = cat(1, abs_theta_dist{:,2});
    rho_mean = mean(rho,1); rho_mean = [rho_mean, rho_mean(1)];
    rho_std = std(rho, [], 1); rho_std = [rho_std, rho_std(1)];
    num_lines = size(theta,1);
    for l = 1:num_lines
        rho_outer = rho_mean(l,:) + rho_std(l,:);
        rho_inner = rho_mean(l,:) - rho_std(l,:);
        
        [x, y] = pol2cart(theta(l,:), rho_outer);
        [x2, y2] = pol2cart(theta(l,:), rho_inner);
        fill([x, fliplr(x2)],[y, fliplr(y2)], [.6 .6 .6],...
            'EdgeColor', [.6 .6 .6],'FaceAlpha', 0.1, 'EdgeAlpha', 0)
        [x3, y3] = pol2cart(theta(l,:), rho_mean(l,:));
        plot(x3,y3, 'Color', [.6 .6 .6], 'LineWidth', 1)
    end

    % Plot bias data
    theta = [bc, bc(1)];
    rho = cat(1, delta_theta_dist{:,2});
    rho_mean = mean(rho,1); rho_mean = [rho_mean, rho_mean(1)];
    rho_std = std(rho, [], 1); rho_std = [rho_std, rho_std(1)];
    num_lines = size(theta,1);
    for l = 1:num_lines
        rho_outer = rho_mean(l,:) + rho_std(l,:);
        rho_inner = rho_mean(l,:) - rho_std(l,:);
        
        [x, y] = pol2cart(theta(l,:), rho_outer);
        [x2, y2] = pol2cart(theta(l,:), rho_inner);
        fill([x, fliplr(x2)],[y, fliplr(y2)], center_col,...
            'EdgeColor', center_col,'FaceAlpha', 0.1, 'EdgeAlpha', 0)
        [x3, y3] = pol2cart(theta(l,:), rho_mean(l,:));
        plot(x3,y3, 'Color', center_col, 'LineWidth', 1)
    end

annotation("textbox", [0.35 0. 0.3 0.1], 'String', ColorText({'Absolute Angle', 'Center Angle'}, [[.6 .6 .6]; center_col]),...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none')
char_offset = 64;
annotation("textbox", [0.025 0.925 .05 .05], 'String', char(char_offset+1), ...
'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.475 0.925 .05 .05], 'String', char(char_offset+2), ...
'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.025 0.61 .05 .05], 'String', char(char_offset+3), ...
'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.475 0.61 .05 .05], 'String', char(char_offset+4), ...
'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.025 0.285 .05 .05], 'String', char(char_offset+5), ...
'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')
annotation("textbox", [0.475 0.285 .05 .05], 'String', char(char_offset+6), ...
'VerticalAlignment','top', 'HorizontalAlignment','left', 'EdgeColor', 'none', 'FontWeight','bold')

shg
print(gcf, fullfile(FigurePath, "Fig3_Localization.png"), '-dpng', '-r300')

shg
%% Supplement?
clf;
axes('Position', [0.1 0.5 0.3 0.35]); hold on
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

% Breast theta distribution
axes('Position', [0.075 0.075 0.35 0.2]); hold on
    % absolute
    abs_rho = cat(1, abs_theta_dist{:,1});
    for l = 1:size(abs_rho,1)
        plot(bc, abs_rho(l,:), 'Color', [.6 .6 .6], 'LineWidth', 1)
    end
    % delta
    delta_rho = cat(1, delta_theta_dist{:,1});
    for l = 1:size(delta_rho,1)
        plot(bc,delta_rho(l,:), 'Color', rgb(239, 154, 154), 'LineWidth', 1)
    end
    % means
    plot(bc, mean(abs_rho, 1), 'Color', 'k', 'LineWidth', 2)
    plot(bc, mean(delta_rho, 1), 'Color', center_col, 'LineWidth', 2)
    set(gca, 'XLim', [-pi pi],...
             'XTick', [-pi 0 pi], ...
             'XTickLabels', {sprintf('-%s', GetUnicodeChar('pi')), '0', sprintf('%s', GetUnicodeChar('pi'))})

% back theta distribution
axes('Position', [0.5 0.075 0.35 0.2]); hold on
    % absolute
    abs_rho = cat(1, abs_theta_dist{:,2});
    for l = 1:size(abs_rho,1)
        plot(bc, abs_rho(l,:), 'Color', [.6 .6 .6], 'LineWidth', 1)
    end
    % delta
    delta_rho = cat(1, delta_theta_dist{:,2});
    for l = 1:size(delta_rho,1)
        plot(bc,delta_rho(l,:), 'Color', rgb(239, 154, 154), 'LineWidth', 1)
    end
    % means
    plot(bc, mean(abs_rho, 1), 'Color', 'k', 'LineWidth', 2)
    plot(bc, mean(delta_rho, 1), 'Color', center_col, 'LineWidth', 2)
    set(gca, 'XLim', [-pi pi],...
             'XTick', [-pi 0 pi], ...
             'XTickLabels', {sprintf('-%s', GetUnicodeChar('pi')), '0', sprintf('%s', GetUnicodeChar('pi'))})

shg
%%
    % Plot bias data
    theta = [bc, bc(1)];
    rho = cat(1, delta_theta_dist{:,1});
    rho_mean = mean(rho,1); rho_mean = [rho_mean, rho_mean(1)];
    rho_std = std(rho, [], 1); rho_std = [rho_std, rho_std(1)];
    num_lines = size(theta,1);
    for l = 1:num_lines
        rho_outer = rho_mean(l,:) + rho_std(l,:);
        rho_inner = rho_mean(l,:) - rho_std(l,:);
        
        [x, y] = pol2cart(theta(l,:), rho_outer);
        [x2, y2] = pol2cart(theta(l,:), rho_inner);
        fill([x, fliplr(x2)],[y, fliplr(y2)], center_col,...
            'EdgeColor', center_col,'FaceAlpha', 0.1, 'EdgeAlpha', 0)
        [x3, y3] = pol2cart(theta(l,:), rho_mean(l,:));
        plot(x3,y3, 'Color', center_col, 'LineWidth', 1)
    end


% Back polar plot 
axes('Position', [0.575 0.075 0.3 0.2], 'XColor', 'none', 'YColor', 'none'); hold on
    % Background
    rho_max = 0.2;
    rho_ticks = 0.2;
    rho_labels = 0.2;
    theta_ticks = deg2rad([0:45:315]);
    theta_labels = deg2rad([0:90:270]);

    set(gca, 'XLim', [-rho_max rho_max],...
             'YLim', [-rho_max rho_max],...
             'DataAspectRatio', [1 1 1])
    
    % Theta ticks
    [x, y] = pol2cart(theta_ticks, repmat(rho_max*1.1, [1, length(theta_ticks)]));
    for i = 1:length(theta_ticks)
       plot([0,x(i)], [0, y(i)], 'Color', [0.8 0.8 0.8], 'LineWidth', 0.1)
    end
    
    % Theta labels
    if isempty(theta_labels) == 0
        [x, y] = pol2cart(theta_labels, repmat(rho_max*1.25, [1, length(theta_labels)]));
        for t = 1:length(theta_labels)
            text(x(t),y(t), [num2str(rad2deg(theta_labels(t))), char(176)],...
                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center')
        end
    end
    
    % Rho ticks
    t = linspace(0, 2*pi, 1000);
    for r = rho_ticks
        [x, y] = pol2cart(t, repmat(r, [1, 1000]));
        plot(x, y, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.1)
    end
    
    % Rho labels
    rt = theta_ticks(2) / 2;
    ro = rho_max * 0.05;
    if isempty(rho_labels) == 0
        [x, y] = pol2cart(repmat(rt, [1,length(rho_labels)]), rho_labels+ro);
        for t = 1:length(rho_labels)
            text(x(t),y(t), num2str(rho_labels(t)),...
                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left')
        end
    end
    
    % Plot error data
    theta = [bc, bc(1)];
    rho = cat(1, abs_theta_dist{:,2});
    rho_mean = mean(rho,1); rho_mean = [rho_mean, rho_mean(1)];
    rho_std = std(rho, [], 1); rho_std = [rho_std, rho_std(1)];
    num_lines = size(theta,1);
    for l = 1:num_lines
        rho_outer = rho_mean(l,:) + rho_std(l,:);
        rho_inner = rho_mean(l,:) - rho_std(l,:);
        
        [x, y] = pol2cart(theta(l,:), rho_outer);
        [x2, y2] = pol2cart(theta(l,:), rho_inner);
        fill([x, fliplr(x2)],[y, fliplr(y2)], [.6 .6 .6],...
            'EdgeColor', [.6 .6 .6],'FaceAlpha', 0.1, 'EdgeAlpha', 0)
        [x3, y3] = pol2cart(theta(l,:), rho_mean(l,:));
        plot(x3,y3, 'Color', [.6 .6 .6], 'LineWidth', 1)
    end

    % Plot bias data
    theta = [bc, bc(1)];
    rho = cat(1, delta_theta_dist{:,2});
    rho_mean = mean(rho,1); rho_mean = [rho_mean, rho_mean(1)];
    rho_std = std(rho, [], 1); rho_std = [rho_std, rho_std(1)];
    num_lines = size(theta,1);
    for l = 1:num_lines
        rho_outer = rho_mean(l,:) + rho_std(l,:);
        rho_inner = rho_mean(l,:) - rho_std(l,:);
        
        [x, y] = pol2cart(theta(l,:), rho_outer);
        [x2, y2] = pol2cart(theta(l,:), rho_inner);
        fill([x, fliplr(x2)],[y, fliplr(y2)], center_col,...
            'EdgeColor', center_col,'FaceAlpha', 0.1, 'EdgeAlpha', 0)
        [x3, y3] = pol2cart(theta(l,:), rho_mean(l,:));
        plot(x3,y3, 'Color', center_col, 'LineWidth', 1)
    end

shg

%%
clf;
i = 7;
subplot(1,2,1); hold on
    point_grid = td_localization_data(i).breast_grid;
    resp_grid = td_localization_data(i).breast_resp;
    center_xyz = mean([point_grid.x(1:2), point_grid.y(1:2), point_grid.z(1:2)], 1);
    scatter(center_xyz(1), center_xyz(2), 50, 'k')
    for p = 1:height(point_grid)
        % Get matching points
        idx = resp_grid.Point == point_grid.Point(p);
        point_xyz = [point_grid.x(p), point_grid.y(p), point_grid.z(p)];
        resp_xyz = [resp_grid.x(idx), resp_grid.y(idx), resp_grid.z(idx)];
        resp_mean = mean(resp_xyz, 1); % Get mean response
        b = point_xyz - center_xyz;
        for j = 1:size(resp_xyz, 1) % Plot individual points
            plot([point_xyz(1), resp_xyz(j,1)], [point_xyz(2), resp_xyz(j,2)], 'Color', imp_col)
        end
        % Plot mean
        plot([point_xyz(1), resp_mean(1)], [point_xyz(2), resp_mean(2)], 'Color', bias_col)
    end

subplot(1,2,2); hold on
    point_grid = td_localization_data(i).back_grid;
    resp_grid = td_localization_data(i).back_resp;
    center_xyz = mean([point_grid.x(1:2), point_grid.y(1:2), point_grid.z(1:2)], 1);
    scatter(center_xyz(1), center_xyz(2), 50, 'k')
    for p = 1:height(point_grid)
        % Get matching points
        idx = resp_grid.Point == point_grid.Point(p);
        point_xyz = [point_grid.x(p), point_grid.y(p), point_grid.z(p)];
        resp_xyz = [resp_grid.x(idx), resp_grid.y(idx), resp_grid.z(idx)];
        resp_mean = mean(resp_xyz, 1); % Get mean response
        b = point_xyz - center_xyz;
        for j = 1:size(resp_xyz, 1) % Plot individual points
            plot([point_xyz(1), resp_xyz(j,1)], [point_xyz(2), resp_xyz(j,2)], 'Color', imp_col)
        end
        % Plot mean
        plot([point_xyz(1), resp_mean(1)], [point_xyz(2), resp_mean(2)], 'Color', bias_col)
    end
