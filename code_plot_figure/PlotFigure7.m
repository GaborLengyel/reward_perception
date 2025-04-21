%% Figure7 A
load('simulation_data/Figure7A');
subplot(2, 2, 1)
% data plot
errorbar([2, 4, 8, 16, 32], m_did_D(:,2,100), m_did_D(:,2,100) -  d_did_D(:,2,100), u_did_D(:,2,100) -  m_did_D(:,2,100), 'bo', 'LineWidth', 1, 'MarkerSize', 5)
hold on
errorbar([2, 4, 8, 16, 32], m_did_P_1(:,2,100), m_did_P_1(:,2,100) -  d_did_P_1(:,2,100), u_did_P_1(:,2,100) -  m_did_P_1(:,2,100), 'go', 'LineWidth', 1, 'MarkerSize', 5)
errorbar([2, 4, 8, 16, 32], m_did_P_2(:,2,100), m_did_P_2(:,2,100) -  d_did_P_2(:,2,100), u_did_P_2(:,2,100) -  m_did_P_2(:,2,100), 'ro', 'LineWidth', 1, 'MarkerSize', 5)
% starting porsition
plot(0, 0, 'bo', 'LineWidth', 1, 'MarkerSize', 5)
plot(0, -10, 'go', 'LineWidth', 1, 'MarkerSize', 5)
plot(0, 10, 'ro', 'LineWidth', 1, 'MarkerSize', 5)
% Ground Truth porsition
plot([0, 2, 4, 8, 16, 32], 10 * ones(6, 1), 'b--', 'LineWidth', 1)
plot([0, 2, 4, 8, 16, 32], -20 * ones(6, 1), 'g--', 'LineWidth', 1)
plot([0, 2, 4, 8, 16, 32], 20 * ones(6, 1), 'r--', 'LineWidth', 1)
set(gca, 'xscale', 'log')
xlim([0, 32])
set(gca, 'FontSize', 10,  'LineWidth', 1)
box off
xlabel('Prior Width (\sigma, deg)')
ylabel('Bias Estimation after 100 trials (deg')
%% Figure 7B
load('analysis_data/Figure7BCD&FigureS2.mat')
subplot(2, 2, 2)
k=[];
k(1,:) = [0 1 0];
k(2,:) = [0 0 1];
k(3,:) = [1 0 0];
kk(1, :) = [0.8, 1, 0.8];
kk(2, :) = [0.8, 0.8, 1];
kk(3, :) = [1, 0.8, 0.8];

for kkk = 1
    
    for i = 1: 3
        
        % mean
        plot(tau_m(:,i), 'Color', k(i, :), 'LineWidth', 2);

        hold on
        % error band
        p = fill([1: length(tau_m(:,i)) length(tau_m(:,i)): -1: 1], ...
        [tau_ul(:,i, 2)', tau_ul(end:-1:1, i, 1)'], 'black');
        p.FaceColor = kk(i, :);
        p.EdgeColor = 'none';
        p.FaceAlpha= 0.5;
        
    end
    
end
xlabel('Session Number')
ylabel('Prior Width (deg)')
set(gca, 'LineWidth', 1, 'FontSize', 10)
box off

%% Figure 7C

subplot(2, 2, 3)
errorbar(betae1_plot_mean * ecc_all + betan1_plot_mean * angle_all, P1_plot_mean, ...
    P1_plot_mean - P1_plot_low, P1_plot_up - P1_plot_mean, ...
    betae1_plot_mean * ecc_all + betan1_plot_mean * angle_all - ...
    betae1_plot_low * ecc_all - betan1_plot_low * angle_all,...
    betae1_plot_up * ecc_all + betan1_plot_up * angle_all - ...
    betae1_plot_mean * ecc_all - betan1_plot_mean * angle_all, 'go', 'LineWidth', 1, 'MarkerSize', 5)
hold on
errorbar(betae2_plot_mean * ecc_all + betan2_plot_mean * angle_all, P2_plot_mean, ...
    P2_plot_mean - P2_plot_low, P2_plot_up - P2_plot_mean, ...
    betae2_plot_mean * ecc_all + betan2_plot_mean * angle_all - ...
    betae2_plot_low * ecc_all - betan2_plot_low * angle_all,...
    betae2_plot_up * ecc_all + betan2_plot_up * angle_all - ...
    betae2_plot_mean * ecc_all - betan2_plot_mean * angle_all, 'ro', 'LineWidth', 1, 'MarkerSize', 5)

errorbar(zeros(length(P1_plot_mean), 1), D_plot_mean, ...
    D_plot_mean - D_plot_low, D_plot_up - D_plot_mean, zeros(length(D_plot_mean), 1),...
    zeros(length(D_plot_mean), 1), 'bo', 'LineWidth', 1, 'MarkerSize', 5)

plot(-45: 40, -45:40, 'k--', 'LineWidth', 1)
xlabel('Prior Mean (deg)')
ylabel('Posterior Mean (deg)')

%% Figure 7D
subplot(2, 2, 4)
plot(P1_plot_mean - betae1_plot_mean * ecc_all - betan1_plot_mean * angle_all,  'go', 'LineWidth', 1, 'MarkerSize', 5)
hold on
plot(P2_plot_mean - betae2_plot_mean * ecc_all - betan2_plot_mean * angle_all,  'ro', 'LineWidth', 1, 'MarkerSize', 5)

plot(0:55, zeros(56, 1), 'k--', 'LineWidth', 1)

xlabel('Session Number')
ylabel('Perceptual Bias - Long-Term Mean')