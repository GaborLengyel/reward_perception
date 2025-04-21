%% Figure ploting for Figure S6
load('simulation_data/FigureS6A.mat');
k=[];
k(1,:) = [0 1 0];
k(2,:) = [0 0 1];
k(3,:) = [1 0 0];
kk(1, :) = [0.8, 1, 0.8];
kk(2, :) = [0.8, 0.8, 1];
kk(3, :) = [1, 0.8, 0.8];
subplot(3, 1, 1)
for kkk = 1

    b = B_s{1, 2, kkk};
    c = C_s{1, 2, kkk};
    
    for i = 1: 3
        
        % mean
        plot(b(:,i), 'Color', k(i, :), 'LineWidth', 2);

        hold on
        % error band
        p = fill([1: length(b(:,i)) length(b(:,i)): -1: 1], ...
        [c(:,i, 2)', c(end:-1:1, i, 1)'], 'black');
        p.FaceColor = kk(i, :);
        p.EdgeColor = 'none';
        p.FaceAlpha= 0.5;
        
    end
    % ground truth
    plot(1: length(b(:,i)), -10 * ones(length(b(:,i)), 1), 'g--',  'LineWidth', 2)
    plot(1: length(b(:,i)), 10 * ones(length(b(:,i)), 1), 'b--',  'LineWidth', 2) % for constant decision bias
    plot(1: length(b(:,i)), 20 * ones(length(b(:,i)), 1), 'r--',  'LineWidth', 2)
    % zero
    plot(1: length(b(:,i)), 0 * ones(length(b(:,i)), 1), 'k--',  'LineWidth', 2)
end
xlim([20, 1000])
set(gca, 'xscale', 'log')
xlabel('Trial Number')
ylabel('Perceptual/Decision Bias (deg)')
set(gca, 'LineWidth', 1, 'FontSize', 10)
box off

%% Figure ploting for Figure S6B

clear
clc

gt1 = -10;

color_plot(3, :) = [0, 1/3, 0]; 
color_plot(2, :) = [0, 2/3, 0]; 
color_plot(1, :) = [0, 1, 0]; 


% data process for Figure S6BC
bb = [];
load('simulation_data/Figure5CE_informative.mat');
for i = 1: 5
    for j = 1: 12
        bb(i, j, :) = B_s{j, 2, i}(:,1);
    end
end
bb = reshape(bb, 20, 3, []);
bb = permute(bb, [3, 2, 1]);
for i = 1: 3
    dis_P_1_m(i,2,:, :) = bb(:,i,:);
end

% data process for Figure S6BC
bb = [];
load('simulation_data/FigureS6BC.mat');
for i = 1: 5
    for j = 1: 12
        bb(i, j, :) = B_s{j, 2, i}(:,1);
    end
end
bb = reshape(bb, 20, 3, []);
bb = permute(bb, [3, 2, 1]);
for i = 1: 3
    se_P_1_m(i,2,:, :) = bb(:,i,:);
end

windowSize = 40;
subplot(3, 1, 2)
% non sequential effect in figure S6C
for i = 1: 3
    hold on
    plot(33: 990, smoothdata(squeeze(mean(sqrt((dis_P_1_m(i, 2, 33: end, :) -gt1).^2), 4)), 'gaussian', windowSize), '--','Color', color_plot(i, :), 'LineWidth', 1)
end
% sequential effect in figure S6C
for i = 1: 3
    hold on
    plot(33: 990, smoothdata(squeeze(mean(sqrt((se_P_1_m(i, 2, 33: end, :) -gt1).^2), 4)), 'gaussian', windowSize), '-','Color', color_plot(i, :), 'LineWidth', 1)
end
set(gca, 'xscale', 'log')
xlim([33, 1000])
ylabel('RMSE(deg)')
xlabel('Trial Number')
% non sequential effect in figure S6C
subplot(3, 1, 3)
for i = 1: 3
    hold on
    plot(33: 990, smoothdata(squeeze(std(sqrt((dis_P_1_m(i, 2, 33: end, :) -gt1).^2),[], 4)), 'gaussian', windowSize), '--','Color', color_plot(i, :), 'LineWidth', 1)
end
% sequential effect in figure S6C
for i = 1: 3
    hold on
    plot(33: 990, smoothdata(squeeze(std(sqrt((se_P_1_m(i, 2, 33: end, :) -gt1).^2),[], 4)), 'gaussian', windowSize), '-','Color', color_plot(i, :), 'LineWidth', 1)
end
set(gca, 'xscale', 'log')
xlim([33, 1000])
ylabel('SD(deg)')
xlabel('Trial Number')
