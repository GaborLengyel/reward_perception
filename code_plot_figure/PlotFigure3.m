%% Figure ploting for Figure 3
load('simulation_data/Figure3&FigureS3.mat');
k = 1;
figure()
for method = [1:3, 7] % method = 1:3 Figure4 A-F; method = 7 Figure 4 GH;
    
    % Figure for the changing of decision criteria (Figure 4 ACEG)
    subplot(4, 2, k)
    plot(1:50 * 990, 0 * ones(50 * 990,1), 'k--', 'LineWidth', 0.75)
    hold on
    bias_ci = bias_method_learn{method};

    p = fill([1: 50*990 50*990: -1: 1], ...
    [mean(squeeze(bias_ci(:, 1, :)), 2)' - std(squeeze(bias_ci(:, 1, :)), [], 2)', ...
    mean(squeeze(bias_ci(end: -1: 1, 1, :)), 2)' + std(squeeze(bias_ci(end: -1: 1, 1, :)), [], 2)'], 'black');
    p.FaceColor = [0.8, 1, 0.8];
    p.EdgeColor = 'none';
    p.FaceAlpha= 0.5; % there is only one 1 line with shared decision criteria (figure S5 EFGH)

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
    
    xlabel('Trail Number')
    ylabel('Decision Critera (deg)')
    ylim([-25, 25])
    xlim([1, 5 * 10^4])
    set(gca, 'LineWidth', 1, 'FontSize',10)
    box off
    k = k+ 1;
    
    % Figure for the fitted empricial bias (Figure 4 BDFH)
    subplot(4, 2, k)
    bias_ci = bias_method_psig{method};
 
    plot(1:50, -10 * ones(50,1), 'g--', 'LineWidth', 0.75)
    hold on
    plot(1:50, 0 * ones(50,1), 'b--', 'LineWidth', 0.75)
    plot(1:50, 20 * ones(50,1), 'r--', 'LineWidth', 0.75)
    
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