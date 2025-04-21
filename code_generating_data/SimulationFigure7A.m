%% Figure7A
clear
clc
% prior mean for perceptual bias
intial_p1 = -10; % leftward
intial_p2 = 10; % rightward

intial_d = 0;

% prior std for decision & perceptual bias (Figure7 A)
std_in_d = [2, 4, 8, 16, 32]; % decision
std_in_p = [2, 4, 8, 16, 32]; % perceptual

num_trial = 100; % number of trials
gt_b = [-20, 0, 20]; % ground truth percpetual bias (Figure7 A)

for k = 1:20
    seed = sum(round(clock));
    for j = 2 % constant or changning decision bias: 1: changning, 2: constant

        for num_all = 1: length(std_in_p)

                std_p1 = std_in_p(num_all);
                std_p2 = std_in_p(num_all);
                std_d = std_in_d(num_all);
                [data, heading_index] = simu_stan(gt_b, 16, j - 1, num_trial, seed); % generating stimulus & choice data

                iter = 5000; %number of sampling plus warm-up
                options = struct;
                options.poolMaxGap     = inf;
                options.poolMaxLength  = inf;
                options.poolxTol       = 0;
                for i = 1: num_trial
                    i;
                    if i >= 33 % update after 33 trials
                         for hra = 1:3
                            num_deading = length(find(heading_index(1: i) == hra));
                            data_stan(:,:,hra) = poolData(data{hra}(1: num_deading, :), options);           
                         end
                    % initialize for stan fitting for infromative prior
                     monk_dat = struct('ndir',11, ...
                                      'direc',squeeze(data_stan(:,1,:)),...
                                      'choi',squeeze(data_stan(:,2,:)),...
                                      'n', squeeze(data_stan(:,3,:)),...
                                      'P0_1', intial_p1(1),...
                                      'P0_2', intial_p2(1),...
                                      'D0', intial_d, ...
                                      'lapse_alpha_1', 1,...
                                      'lapse_beta_1', 10,...
                                      'lapse_alpha_2', 1,...
                                      'lapse_beta_2', 10,...
                                      'phi_alpha', 8,...
                                      'phi_beta', 0.5,...
                                      'tau_1', std_p1,...
                                      'tau_2', std_p2,...
                                      'tau_d', std_d);    
                                  
                        params = struct('file','StanSimulation_Decision.stan','data',monk_dat,'iter',iter,'chains',1);
                        fit = stan(params, 'init', struct('P_1', intial_p1(1), 'P_2', intial_p2(1), 'D',intial_d, 'phi', [16, 16, 16]));
                        waitfor(fit,'exit_value',0);
                        para_dis = fit.extract();
                        b(i, 1) = quantile(para_dis.P_1, 0.5);
                        b(i, 3) = quantile(para_dis.P_2, 0.5);
                        b(i, 2) = quantile(para_dis.D, 0.5);
                        c(i, 1, 1) = quantile(para_dis.P_1, 0.16);
                        c(i, 3, 1) = quantile(para_dis.P_2, 0.16);
                        c(i, 2, 1) = quantile(para_dis.D, 0.16);
                        c(i, 1, 2) = quantile(para_dis.P_1, 0.84);
                        c(i, 3, 2) = quantile(para_dis.P_2, 0.84);
                        c(i, 2, 2) = quantile(para_dis.D, 0.84);
                        
                        para_dis_all{i} = para_dis;
                        
                    else % before 33 trials, it will use the piror value
                        b(i, 1) = intial_p1(1);
                        b(i, 3) = intial_p2(1);
                        b(i, 2) = intial_d;
                        c(i, 1, 1) = intial_p1(1) - std_p1;
                        c(i, 3, 1) = intial_p2(1) - std_p2;
                        c(i, 2, 1) = intial_d-std_d;
                        c(i, 1, 2) = intial_p1(1) + std_p1;
                        c(i, 3, 2) = intial_p2(1) + std_p2;
                        c(i, 2, 2) = intial_d+std_d;
                    end
                end
                
                name_save = sprintf('save_onlie_stan_%d_%d_%d', num_all, j, k);
                save(name_save, 'para_dis_all');

        end
    end
end

for i = 1:5
for j = 2
    for k = 1:20
        name = sprintf('save_onlie_stan_%d_%d_%d.mat', i, j, k);
        result = load(name);
        overall = result.para_dis_all;
            for t = 33: 100
                dis_P_1(i, j, t, (k -1)* 2500 + 1 : 2500 * k) = overall{t}.P_1;
                dis_D(i, j, t, (k -1)* 2500 + 1 : 2500 * k)  = overall{t}.D;
                dis_P_2(i, j, t, (k -1)* 2500 + 1 : 2500 * k)  = overall{t}.P_2;
            end                
    end
    for t = 33: 100
        m_did_P_1(i, j, t) = quantile(dis_P_1(i, j, t, :), 0.5);
        m_did_D(i, j, t) = quantile(dis_D(i, j, t, :), 0.5);
        m_did_P_2(i, j, t) = quantile(dis_P_2(i, j, t, :), 0.5);

        u_did_P_1(i, j, t) = quantile(dis_P_1(i, j, t, :), 0.16);
        u_did_D(i, j, t) = quantile(dis_D(i, j, t, :), 0.16);
        u_did_P_2(i, j, t) = quantile(dis_P_2(i, j, t, :), 0.16);

        d_did_P_1(i, j, t) = quantile(dis_P_1(i, j, t, :), 0.84);
        d_did_D(i, j, t) = quantile(dis_D(i, j, t, :), 0.84);
        d_did_P_2(i, j, t) = quantile(dis_P_2(i, j, t, :), 0.84);
    
    end
end
end