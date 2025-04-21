%% Figure6
%% Figure 6B upper
load('analysis_data/Figure6B1');
color_plot(1, :) = [0, 1, 0]; 
color_plot(2, :) = [0, 0, 1]; 
color_plot(3, :) = [1, 0, 0]; 

% plot raw data
subplot(2, 2, 1)
plot(data_stan(:, 1,1), data_stan(:, 2,1)./data_stan(:, 3,1), 'go', 'LineWidth', 1, 'MarkerSize', 5)
hold on
plot(data_stan(:, 1,2), data_stan(:, 2,2)./data_stan(:, 3,2), 'bsquare', 'LineWidth', 1, 'MarkerSize', 5)
plot(data_stan(:, 1,3), data_stan(:, 2,3)./data_stan(:, 3,3), 'r^', 'LineWidth', 1, 'MarkerSize', 5)


% plot fitted curve
curve_x = linspace(min(data_stan(:, 1,1)), max(data_stan(:, 1,1)), 100);
for i = 1: 3
    curve_y = lapse1 + (1 - lapse1 - lapse2) * normcdf(curve_x, P_b(i) + D_b, thre_f(i));
    plot(curve_x, curve_y, 'Color', color_plot(i, :), 'LineWidth', 1)
    hold on
end
plot(curve_x, 0.5 * ones(length(curve_x), 1), 'k--', 'LineWidth', 1)
plot(zeros(length(curve_x), 1), linspace(0, 1,length(curve_x)), 'k--', 'LineWidth', 1)

xlabel('Object Motion Direction')
ylabel('Proportion Rightward Choice')
set(gca, 'LineWidth', 1, 'FontSize', 10)
box off

%% Figure 6C upper
load('analysis_data/Figure6C1');
subplot(2, 2, 2)
k=[];
k(1,:) = [0 1 0];
k(2,:) = [0 0 1];
k(3,:) = [1 0 0];
kk(1, :) = [0.8, 1, 0.8];
kk(2, :) = [0.8, 0.8, 1];
kk(3, :) = [1, 0.8, 0.8];
for kkk = 1

    b = B_s{kkk};
    c = C_s{kkk};

    for i = 1: 3
        
        % mean
        plot(b(:,i), 'Color', k(i, :), 'LineWidth', 1);

        hold on
        % error band
        p = fill([1: length(b(:,i)) length(b(:,i)): -1: 1], ...
        [c(:,i, 2)', c(end:-1:1, i, 1)'], 'black');
        p.FaceColor = kk(i, :);
        p.EdgeColor = 'none';
        p.FaceAlpha= 0.5;
        
    end
    
    b = B_f{kkk};
    for i = 1: 3
        
        
        plot(b(:,i),'.-','Color', 0.8 * k(i, :), 'LineWidth', 0.75);
    end
   
    plot(1: length(b(:,i)), 0 * ones(length(b(:,i)), 1), 'k--',  'LineWidth', 1)
end
set(gca, 'xscale', 'log')
xlabel('Trial')
ylabel('Perceptual/Decision Bias (deg)')
set(gca, 'LineWidth', 1, 'FontSize', 10)
box off
xlim([20, 850]);
ylim([-45, 30]);
%% Figure 6B lower

load('analysis_data/Figure6B2');

color_plot(1, :) = [0, 1, 0]; 
color_plot(2, :) = [0, 0, 1]; 
color_plot(3, :) = [1, 0, 0]; 

% plot raw data
subplot(2, 2, 3)
plot(data_stan(:, 1,1), data_stan(:, 2,1)./data_stan(:, 3,1), 'go', 'LineWidth', 1, 'MarkerSize', 5)
hold on
plot(data_stan(:, 1,2), data_stan(:, 2,2)./data_stan(:, 3,2), 'bsquare', 'LineWidth', 1, 'MarkerSize', 5)
plot(data_stan(:, 1,3), data_stan(:, 2,3)./data_stan(:, 3,3), 'r^', 'LineWidth', 1, 'MarkerSize', 5)


% plot fitted curve
curve_x = linspace(min(data_stan(:, 1,1)), max(data_stan(:, 1,1)), 100);
for i = 1: 3
    curve_y = lapse1 + (1 - lapse1 - lapse2) * normcdf(curve_x, P_b(i) + D_b, thre_f(i));
    plot(curve_x, curve_y, 'Color', color_plot(i, :), 'LineWidth', 1)
    hold on
end
plot(curve_x, 0.5 * ones(length(curve_x), 1), 'k--', 'LineWidth', 1)
plot(zeros(length(curve_x), 1), linspace(0, 1,length(curve_x)), 'k--', 'LineWidth', 1)

xlabel('Object Motion Direction')
ylabel('Proportion Rightward Choice')
set(gca, 'LineWidth', 1, 'FontSize', 10)
box off

%% Figure 6C lower
load('analysis_data/Figure6C2');
subplot(2, 2, 4)
k=[];
k(1,:) = [0 1 0];
k(2,:) = [0 0 1];
k(3,:) = [1 0 0];
kk(1, :) = [0.8, 1, 0.8];
kk(2, :) = [0.8, 0.8, 1];
kk(3, :) = [1, 0.8, 0.8];
for kkk = 1

    b = B_s{kkk};
    c = C_s{kkk};

    for i = 1: 3
        
        % mean
        plot(b(:,i), 'Color', k(i, :), 'LineWidth', 1);

        hold on
        % error band
        p = fill([1: length(b(:,i)) length(b(:,i)): -1: 1], ...
        [c(:,i, 2)', c(end:-1:1, i, 1)'], 'black');
        p.FaceColor = kk(i, :);
        p.EdgeColor = 'none';
        p.FaceAlpha= 0.5;
        
    end
    b = B_f{kkk};
    for i = 1: 3
        
        
        plot(b(:,i), 'k-.', 'Color', 0.8 * k(i, :), 'LineWidth', 0.75);
    end
    plot(1: length(b(:,i)), 0 * ones(length(b(:,i)), 1), 'k--',  'LineWidth', 1)
end
set(gca, 'xscale', 'log')
xlabel('Trial Number')
ylabel('Perceptual/Decision Bias (deg)')
set(gca, 'LineWidth', 1, 'FontSize', 10)
box off
xlim([20, 850]);
ylim([-30, 30]);