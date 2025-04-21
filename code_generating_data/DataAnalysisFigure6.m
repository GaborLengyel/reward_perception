%% Monkey data anaysis for Figure6
clear
clc
data_fit = cell(3, 1);
% load data
% load('monkey_data/monkey_training_data_session_11.mat') % lower session
load('monkey_data/monkey_training_data_session_14.mat') % upper session
heading_condition = unique(data.heading_direction);
heading_condition = heading_condition(1: 3);
direction_condition = unique(data.stimulus_direction);
for i = 1: length(data.misc_params)
    if data.misc_params(i)>=0 % valide trial only
        heading = find(data.heading_direction(i) == heading_condition); % heading direction (left, neutral, right) for each trial
        direction = find(data.stimulus_direction(i) == direction_condition); % stimulus direction for each trial
        choice = data.choice(i); % monkey's choice for each trial        
        data_up = [direction_condition(direction), choice, 1];
        data_fit{heading}(end + 1, :) = data_up;
        heading_index(i) = heading;
    end
end
%% informative prior
intial = 10; 
% intial = 16 for lower panel
prior_width_p = 7;
prior_width_d = 10;
for num_all = 1: length(prior_width_p)



    num_trial = length(data_fit{1}(:,1)) + length(data_fit{1}(:,2))+  length(data_fit{1}(:,3));
    

    initi = intial;
    iter = 5000;
    options = struct;
    options.poolMaxGap     = inf;
    options.poolMaxLength  = inf;
    options.poolxTol       = 0;
    for i = 1: num_trial
        if i >= 33
             for hra = 1:3
                num_deading = length(find(heading_index(1: i) == hra));
                data_stan(:,:,hra) = poolData(data_fit{hra}(1: num_deading, :), options);           
             end
             
        % for infromative prior 
         monk_dat = struct('ndir',11, ...
                          'direc',squeeze(data_stan(:,1,:)),...
                          'choi',squeeze(data_stan(:,2,:)),...
                          'n', squeeze(data_stan(:,3,:)),...
                          'P0_1', -initi,...
                          'P0_2', initi,...
                          'lapse_alpha_1', 1,...
                          'lapse_beta_1', 10,...
                          'lapse_alpha_2', 1,...
                          'lapse_beta_2', 10,...
                          'phi_alpha', 8,...
                          'phi_beta', 0.5,...
                          'tau_1', prior_width_p,...
                          'tau_2', prior_width_p,...
                          'tau_d', prior_width_d);
                      
            params = struct('file','StanSimulation.stan','data',monk_dat,'iter',iter,'chains',1);

            fit = stan(params, 'init', struct('P_1', -initi, 'P_2', initi, 'D',0, 'phi', [16, 16, 16]));
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

        else
            b(i, 1) = -initi;
            b(i, 3) = initi;
            b(i, 2) = 0;
            c(i, 1, 1) = -initi - prior_width_p;
            c(i, 3, 1) = initi - prior_width_p;
            c(i, 2, 1) = -prior_width_d;
            c(i, 1, 2) = -initi + prior_width_p;
            c(i, 3, 2) = initi + prior_width_p;
            c(i, 2, 2) = prior_width_d;
        end
    end
    B_s{num_all} = b;
    C_s{num_all} = c;
    p_s{num_all} = para_dis;
    D_s{num_all} = data_stan;
end

lapse1 = quantile(p_s{1}.lapse_1, 0.5);
lapse2 = quantile(p_s{1}.lapse_2, 0.5);
P_b = [quantile(p_s{1}.P_1, 0.5), 0, quantile(p_s{1}.P_2, 0.5)];
D_b = quantile(p_s{1}.D, 0.5);
thre_f = quantile(p_s{1}.phi, 0.5);
%% flat prior
intial = 10; 
% intial = 16 for lower panel
prior_width_p = 7;
prior_width_d = 10;
for num_all = 1: length(prior_width_p)



    num_trial = length(data_fit{1}(:,1)) + length(data_fit{1}(:,2))+  length(data_fit{1}(:,3));
    

    initi = intial;
    iter = 5000;
    options = struct;
    options.poolMaxGap     = inf;
    options.poolMaxLength  = inf;
    options.poolxTol       = 0;
    for i = 1: num_trial
        if i >= 33
             for hra = 1:3
                num_deading = length(find(heading_index(1: i) == hra));
                data_stan(:,:,hra) = poolData(data_fit{hra}(1: num_deading, :), options);           
             end
                      
        % for flat prior 
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
                          'tau_1', prior_width_p,...
                          'tau_2', prior_width_p,...
                          'tau_d', prior_width_d);

            params = struct('file','StanSimulation_Uniform.stan','data',monk_dat,'iter',iter,'chains',1); % flat prior

            fit = stan(params, 'init', struct('P_1', -initi, 'P_2', initi, 'D',0, 'phi', [16, 16, 16]));
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

        else
            b(i, 1) = -initi;
            b(i, 3) = initi;
            b(i, 2) = 0;
            c(i, 1, 1) = -initi - prior_width_p;
            c(i, 3, 1) = initi - prior_width_p;
            c(i, 2, 1) = -prior_width_d;
            c(i, 1, 2) = -initi + prior_width_p;
            c(i, 3, 2) = initi + prior_width_p;
            c(i, 2, 2) = prior_width_d;
        end
    end
    B_f{num_all} = b;
    C_f{num_all} = c;
    p_f{num_all} = para_dis;
end

