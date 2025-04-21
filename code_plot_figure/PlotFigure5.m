%% Figure ploting for Figure 5
%% Figure ploting for Figure 5A
load('simulation_data/Figure5A.mat');
subplot(3, 2, 1)
k=[];
k(1,:) = [0 1 0];
k(2,:) = [0 0 1];
k(3,:) = [1 0 0];
kk(1, :) = [0.8, 1, 0.8];
kk(2, :) = [0.8, 0.8, 1];
kk(3, :) = [1, 0.8, 0.8];
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
set(gca, 'xscale', 'log')
xlabel('Trial')
ylabel('Perceptual/Decision Bias (deg)')
set(gca, 'LineWidth', 1, 'FontSize', 10)
xlim([20, 1000])
box off
%% Figure ploting for Figure 5B
load('simulation_data/Figure5B.mat');
subplot(3, 2, 2)
k=[];
k(1,:) = [0 1 0];
k(2,:) = [0 0 1];
k(3,:) = [1 0 0];
kk(1, :) = [0.8, 1, 0.8];
kk(2, :) = [0.8, 0.8, 1];
kk(3, :) = [1, 0.8, 0.8];
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
    plot(1: length(b(:,i)), 10 * sin((1: length(b(:,i)))/2000 * pi), 'b--',  'LineWidth', 2) % for changning decision bias
    plot(1: length(b(:,i)), 20 * ones(length(b(:,i)), 1), 'r--',  'LineWidth', 2)
    % zero
    plot(1: length(b(:,i)), 0 * ones(length(b(:,i)), 1), 'k--',  'LineWidth', 2)
end
xlim([20, 1000]);
set(gca, 'xscale', 'log')
xlabel('Trial Number')
ylabel('Perceptual/Decision Bias (deg)')
set(gca, 'LineWidth', 1, 'FontSize', 10)
box off
%% Figure ploting for Figure 5CDEF
clear
clc
gt1 = -10;

color_plot(3, :) = [0, 1/3, 0]; 
color_plot(2, :) = [0, 2/3, 0]; 
color_plot(1, :) = [0, 1, 0]; 

load('simulation_data/Figure5CE_informative.mat');
% data process for Figure 5CE
for i = 1: 25 
    for j = 1: 12
        bb(i, j, :) = B_s{j, 2, i}(:,1);
    end
end
bb = reshape(bb, 100, 3, []);
bb = permute(bb, [3, 2, 1]);
for i = 1: 3
    dis_P_1_m(i,2,:, :) = bb(:,i,:);
end
load('simulation_data/Figure5CE_noninformative.mat');
for i = 1: 100
flat_P_1_m(1,2,:, i) = B_s{1, 2, i}(:,1);
end

% Figure5 CD
% our method in figure 5C
subplot(3, 2, 3)
for i = 1: 3
    hold on
    plot(33: 990, squeeze(mean(sqrt((dis_P_1_m(i, 2, 33: end, :) -gt1).^2), 4)), '-','Color', color_plot(i, :), 'LineWidth', 1)
end
% flat in figure 5C
plot(33: 990, squeeze(mean(sqrt((flat_P_1_m(1, 2, 33: end, :) -gt1).^2), 4)), '--','Color', [0,0,0], 'LineWidth', 1) 
set(gca, 'xscale', 'log')
xlim([33, 1000])
ylabel('RMSE(deg)')
xlabel('Trial Number')
% our method in figure 5D
subplot(3, 2, 5)
for i = 1: 3
    hold on
    plot(33: 990, squeeze(std(sqrt((dis_P_1_m(i, 2, 33: end, :) -gt1).^2),[], 4)), '-','Color', color_plot(i, :), 'LineWidth', 1)
end
% flat in figure 5D
plot(33: 990, squeeze(std(sqrt((flat_P_1_m(1, 2, 33: end, :) -gt1).^2),[], 4)), '--','Color', [0,0,0], 'LineWidth', 1) 
set(gca, 'xscale', 'log')
xlim([33, 1000])
ylabel('SD(deg)')
xlabel('Trial Number')

clear
clc
gt1 = -10;

color_plot(1, :) = [0, 1/3, 0]; 
color_plot(2, :) = [0, 2/3, 0]; 
color_plot(3, :) = [0, 1, 0]; 

% data process for Figure 5DF
load('simulation_data/Figure5DF_stan.mat');
for i = 1: 100
    for j = 1: 3 
        dis_P_1_m(j,2,:, i) = B_s{j, 2, i}(:,1);
    end
end
% data process for Figure 5DF, psignifit
load('simulation_data/Figure5DF_pisignifitDefault.mat');
for i = 1: 100
    for j = 1: 3
        b = B_raw{j, 2, i};
        b_right = b{1,1};
        b_neu = b{2,1};
        psig_P_1_m(j,2,:, i) = b_right - b_neu;
    end
end

% our method in figure 5D
subplot(3, 2, 4)
for i = 1: 3
    hold on
    plot(33: 990, squeeze(mean(sqrt((dis_P_1_m(i, 2, 33: end, :) -gt1).^2), 4)), '-','Color', color_plot(i, :), 'LineWidth', 1)
end
% psignifit in figure 5D
for i = 1: 3
    hold on
    plot(33: 990, squeeze(mean(sqrt((psig_P_1_m(i, 2, 33: end, :) -gt1).^2), 4)), '--','Color', color_plot(i, :), 'LineWidth', 1)
end
set(gca, 'xscale', 'log')
xlim([33, 1000])
ylabel('RMSE(deg)')
xlabel('Trial Number')

% our method in figure 5F
subplot(3, 2, 6)
for i = 1: 3
    hold on
    plot(33: 990, squeeze(std(sqrt((dis_P_1_m(i, 2, 33: end, :) -gt1).^2),[], 4)), '-','Color', color_plot(i, :), 'LineWidth', 1)
end
% psignifit in figure 5F
for i = 1: 3
    hold on
    plot(33: 990, squeeze(std(sqrt((psig_P_1_m(i, 2, 33: end, :) -gt1).^2),[], 4)), '--','Color', color_plot(i, :), 'LineWidth', 1)
end
set(gca, 'xscale', 'log')
xlim([33, 1000])
ylabel('SD(deg)')
xlabel('Trial Number')