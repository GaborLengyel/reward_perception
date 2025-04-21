%% Figure S2
load('analysis_data/Figure7BCD&FigureS2.mat')
% FigureS2 A
subplot(1, 2, 1)
for i = 1: length(angle_all)
    plot(-angle_all(i), P1_plot_mean(i),  'go', 'LineWidth', 1, 'MarkerSize', 5, 'Color', [0, ecc_all(i)/max(ecc_all), 0])
    hold on
    plot(angle_all(i), P2_plot_mean(i),  'ro', 'LineWidth', 1, 'MarkerSize', 5, 'Color', [ecc_all(i)/max(ecc_all),0, 0])
end
xlabel('Heading Direction (deg)')
ylabel('Perceptual Bias (deg)')

% FigureS2 B
subplot(1, 2, 2)
for i = 1: length(angle_all)
    plot(-ecc_all(i), P1_plot_mean(i),  'go', 'LineWidth', 1, 'MarkerSize', 5, 'Color', [0, angle_all(i)/max(ecc_all), 0])
    hold on
    plot(ecc_all(i), P2_plot_mean(i),  'ro', 'LineWidth', 1, 'MarkerSize', 5, 'Color', [angle_all(i)/max(ecc_all), 0, 0])
end
xlabel('Eccentrcity (deg)')
ylabel('Perceptual Bias (deg)')
