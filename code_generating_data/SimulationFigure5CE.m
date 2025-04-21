%% Data Generating for Figure5CE
clear
clc
%% informative prior
% prior mean for perceptual bias (Figure5 CE)
intial_p1 = [-10, -10, -10, -10, -15, -15, -5, -5, -20, -20, 0, 0]; % leftward
intial_p2 = [20, 20, 20, 20, 15, 25, 15, 25, 10, 30, 10, 30]; % rightward 

intial_d = 10; % prior mean for decision bias

% prior std for decision & perceptual bias (Figure5 ABCE)
std_in_d = 10; % decision
std_in_p = 5; % perceptual

num_trial = 990; % number of trials
gt_b = [-10, 0, 20]; % ground truth percpetual bias

for k = 1:25 
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
                B_s{num_all, j, k} = b;
                C_s{num_all, j, k} = c;
                D_s{num_all, j, k} = data;
                name_save = sprintf('save_onlie_infor_%d_%d_%d', num_all, j, k); % you can change the name for saving
                save(name_save,'B_s','C_s', 'D_s');
        end
    end
end

%% noninformative prior
% prior mean for perceptual bias
intial_p1 = -10; % leftward
intial_p2 = 20; % rightward

intial_d = 0; % dummy prior mean for decision bias in flat prior

% prior std for decision & perceptual bias (Figure5 ABCE)
std_in_d = 10; % decision
std_in_p = 5; % perceptual

num_trial = 990; % number of trials
gt_b = [-10, 0, 20]; % ground truth percpetual bias

for k = 1:100 
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
                    %initialize for stan fitting for flat prior        
                    monk_dat = struct('ndir',11, ...
                                      'direc',squeeze(data_stan(:,1,:)),...
                                      'choi',squeeze(data_stan(:,2,:)),...
                                      'n', squeeze(data_stan(:,3,:)),...
                                      'DP_minus', -180,...
                                      'DP_plus', 180,...
                                      'lapse_alpha_1', 1,...
                                      'lapse_beta_1', 10,...
                                      'lapse_alpha_2', 1,...
                                      'lapse_beta_2', 10,...
                                      'phi_alpha', 8,...
                                      'phi_beta', 0.5,...
                                      'tau_1', std_p1,...
                                      'tau_2', std_p2,...
                                      'tau_d', std_d);
                               
                        params = struct('file','StanSimulation_Uniform.stan','data',monk_dat,'iter',iter,'chains',1); %flat
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
                B_s{num_all, j, k} = b;
                C_s{num_all, j, k} = c;
                D_s{num_all, j, k} = data;
                name_save = sprintf('save_onlie_flat_%d_%d_%d', num_all, j, k);
                save(name_save,'B_s','C_s', 'D_s');
        end
    end
end