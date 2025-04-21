%% Figure ploting for Figure S5
load('simulation_data/FigureS5.mat');
k = 1;
for method = 1:4
    
    % Figure for the changing of decision criteria (Figure S5 ACEG)
    subplot(4, 2, k)
    plot(1:50 * 990, 0 * ones(50 * 990,1), 'k--', 'LineWidth', 0.75)
    hold on
    bias_ci = bias_rate_learn{method}; 


    if method == 3 || method == 4
        
        p = fill([1: 50*990 50*990: -1: 1], ...
        [mean(squeeze(bias_ci(:, 1, :)), 2)' - std(squeeze(bias_ci(:, 1, :)), [], 2)', ...
        mean(squeeze(bias_ci(end: -1: 1, 1, :)), 2)' + std(squeeze(bias_ci(end: -1: 1, 1, :)), [], 2)'], 'black');
        p.FaceColor = [0.8, 0.8, 1];
        p.EdgeColor = 'none';
        p.FaceAlpha= 0.5; 
        
        p = fill([1: 50*990 50*990: -1: 1], ...
        [mean(squeeze(bias_ci(:, 1, :)), 2)' - std(squeeze(bias_ci(:, 1, :)), [], 2)', ...
        mean(squeeze(bias_ci(end: -1: 1, 1, :)), 2)' + std(squeeze(bias_ci(end: -1: 1, 1, :)), [], 2)'], 'black');
        p.FaceColor = [1, 0.8, 0.8];
        p.EdgeColor = 'none';
        p.FaceAlpha= 0.5; 
        
        p = fill([1: 50*990 50*990: -1: 1], ...
        [mean(squeeze(bias_ci(:, 1, :)), 2)' - std(squeeze(bias_ci(:, 1, :)), [], 2)', ...
        mean(squeeze(bias_ci(end: -1: 1, 1, :)), 2)' + std(squeeze(bias_ci(end: -1: 1, 1, :)), [], 2)'], 'black');
        p.FaceColor = [0.8, 1, 0.8];
        p.EdgeColor = 'none';
        p.FaceAlpha= 0.5; 

        plot(mean(squeeze(bias_ci(:, 1, :)), 2), 'g', 'LineWidth', 1)
        hold on
        plot(mean(squeeze(bias_ci(:, 1, :)), 2), 'b', 'LineWidth', 1)
        plot(mean(squeeze(bias_ci(:, 1, :)), 2), 'r', 'LineWidth', 1)
        
    elseif method == 1 || method == 2
        
        
        p = fill([1: 50*990 50*990: -1: 1], ...
        [mean(squeeze(bias_ci(:, 1, :)), 2)' - std(squeeze(bias_ci(:, 1, :)), [], 2)', ...
        mean(squeeze(bias_ci(end: -1: 1, 1, :)), 2)' + std(squeeze(bias_ci(end: -1: 1, 1, :)), [], 2)'], 'black');
        p.FaceColor = [0.8, 0.8, 1];
        p.EdgeColor = 'none';
        p.FaceAlpha= 0.5; 
        
        p = fill([1: 50*990 50*990: -1: 1], ...
        [mean(squeeze(bias_ci(:, 2, :)), 2)' - std(squeeze(bias_ci(:, 2, :)), [], 2)',...
        mean(squeeze(bias_ci(end: -1: 1, 2, :)), 2)' + std(squeeze(bias_ci(end: -1: 1, 2, :)), [], 2)'], 'black');
        p.FaceColor = [0.8, 0.8, 1];
        p.EdgeColor = 'none';
        p.FaceAlpha= 0.5;

        p = fill([1: 50*990 50*990: -1: 1], ...
        [mean(squeeze(bias_ci(:, 3, :)), 2)' - std(squeeze(bias_ci(:, 3, :)), [], 2)', ...
        mean(squeeze(bias_ci(end: -1: 1, 3, :)), 2)' + std(squeeze(bias_ci(end: -1: 1, 3, :)), [], 2)'], 'black');
        p.FaceColor = [1, 0.8, 0.8];
        p.EdgeColor = 'none';
        p.FaceAlpha= 0.5;
        plot(mean(squeeze(bias_ci(:, 1, :)), 2), 'g', 'LineWidth', 1)
        hold on
        plot(mean(squeeze(bias_ci(:, 2, :)), 2), 'b', 'LineWidth', 1)
        plot(mean(squeeze(bias_ci(:, 3, :)), 2), 'r', 'LineWidth', 1)
    end
    
    xlabel('Trail Number')
    ylabel('Decision Critera (deg)')
    ylim([-25, 25])
    xlim([1, 5 * 10^4])
    set(gca, 'LineWidth', 1, 'FontSize',10)
    box off
    k = k+ 1;
    
    % Figure for the fitted empricial bias (Figure S5 BDFH)
    bias_all(1, 1:9) = -10;
    bias_all(3, 1:9) = 20;
    bias_all(1, 10: 41) = linspace(-10, -10/3, 32);
    bias_all(3, 10: 41) = linspace(20, 20/3, 32);
    bias_all(1, 42:50) = -10/3;
    bias_all(3, 42:50) = 20/3;
    bias_all(2, 1:50) = 0;
    subplot(4, 2, k)
    bias_ci = bias_rate_psig{method}; 
    plot(1:50, bias_all(1,:), 'g--', 'LineWidth', 0.75)
    hold on
    plot(1:50, bias_all(2,:), 'b--', 'LineWidth', 0.75)
    plot(1:50, bias_all(3,:), 'r--', 'LineWidth', 0.75)
    
    p = fill([1: 50 50: -1: 1], ...
    [mean(squeeze(bias_ci(:, 1, :)), 2)' - std(squeeze(bias_ci(:, 1, :)), [], 2)', ...
    mean(squeeze(bias_ci(end: -1: 1, 1, :)), 2)' + std(squeeze(bias_ci(end: -1: 1, 1, :)), [], 2)'], 'black');
    p.FaceColor = [0.8, 1, 0.8];
    p.EdgeColor = 'none';
    p.FaceAlpha= 0.5;

    p = fill([1: 50 50: -1: 1], ...
    [mean(squeeze(bias_ci(:, 2, :)), 2)' - std(squeeze(bias_ci(:, 2, :)), [], 2)',...
    mean(squeeze(bias_ci(end: -1: 1, 2, :)), 2)' + std(squeeze(bias_ci(end: -1: 1, 2, :)), [], 2)'], 'black');
    p.FaceColor = [0.8, 0.8, 1];
    p.EdgeColor = 'none';
    p.FaceAlpha= 0.5;

    p = fill([1: 50 50: -1: 1], ...
    [mean(squeeze(bias_ci(:, 3, :)), 2)' - std(squeeze(bias_ci(:, 3, :)), [], 2)', ...
    mean(squeeze(bias_ci(end: -1: 1, 3, :)), 2)' + std(squeeze(bias_ci(end: -1: 1, 3, :)), [], 2)'], 'black');
    p.FaceColor = [1, 0.8, 0.8];
    p.EdgeColor = 'none';
    p.FaceAlpha= 0.5;

    plot(mean(squeeze(bias_ci(:, 1, :)), 2), 'g', 'LineWidth', 1)
    hold on
    plot(mean(squeeze(bias_ci(:, 2, :)), 2), 'b', 'LineWidth', 1)
    plot(mean(squeeze(bias_ci(:, 3, :)), 2), 'r', 'LineWidth', 1)
    
    xlabel('Session Number')
    ylabel('Empirical Bias (deg)')
    ylim([-25, 25])
    xlim([1, 50])
    set(gca, 'LineWidth', 1, 'FontSize', 10)
    box off
    k = k+ 1;
    
end