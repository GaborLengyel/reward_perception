%% Data Generating for FigureS5
% It would be safer to run 4 of them seperatly, then form them together as 
% 'bias_rate_learn' and 'bias_rate_psig' to generatet the figure.
%% Data Generating for FigureS5AB: independent decision critera, learning rate = 0.04
clear
clc
intial_p1 = -10; % prior mean for perceptual bias
intial_p2 = 20;
alpha = 0.04; % learning rate
num_trial = 990;
num_session = 50;
std_in_d = 10; % prior std for decision bias
std_in_p = 5; % prior std for perceptual bias
bias_d0 = 0;
% changing ground truth percpetual bias
bias_all(1, 1:9) = -10;
bias_all(3, 1:9) = 20;
bias_all(1, 10: 41) = linspace(-10, -10/3, 32);
bias_all(3, 10: 41) = linspace(20, 20/3, 32);
bias_all(1, 42:50) = -10/3;
bias_all(3, 42:50) = 20/3;
bias_all(2, 1:50) = 0;
thre = 16;
for k = 1
    seed = sum(round(clock));

    for num_all = 1: num_session
        bias = bias_all(:, num_all);
        % update prior after each session
        if num_all == 1
            std_p1 = std_in_p(1);
            std_p2 = std_in_p(1);
            std_d = std_in_d(1);
            mean_p1(num_all) = intial_p1(1);
            mean_p2(num_all) = intial_p2(1);
        else
            std_p1 = tau_intial(num_all - 1, 1);
            std_p2 = tau_intial(num_all - 1, 2);
            std_d = tau_intial(num_all - 1, 3);
            mean_p1(num_all) = p_intial(num_all - 1, 1);
            mean_p2(num_all) = p_intial(num_all - 1, 2);        
        end
        % update decision critera
        % independent decision critera 
        if num_all == 1
            bias_dd = zeros(1, 3);
        else
            bias_dd = bias_learn{k, num_all - 1}(end, :);
        end

        rng(seed)

        heading_condition = [-10, 0, 10];
        direction_condition = -40: 8: 40;
        stimuli_list = [direction_condition, direction_condition, direction_condition];


        data = cell(3, 1);
        b = nan(num_trial, 3);
        c = nan(num_trial, 3, 2);
        heading_index = nan(num_trial, 1);
        data_stan = nan(length(direction_condition), 3, 3);
        bias_updata = nan(num_trial, 3); % independent
        for i = 1: num_trial

            if mod(i, 3 * length(direction_condition)) == 1
                sti_index = 1;
                stimuli = randperm(3 * length( direction_condition));
            end
            if stimuli(sti_index) > (2 * length( direction_condition))
                heading = 1;
            elseif stimuli(sti_index) > (1 * length( direction_condition))
                heading = 2;
            else 
                heading = 3;
            end

            if i == 1
                bias_d = bias_dd;
            else
                bias_d = bias_updata(i - 1, :);  % independent
            end

            perception_ob = normrnd(stimuli_list(stimuli(sti_index)) - bias(heading), thre);
            choice_prob = normcdf(perception_ob, bias_d(heading), thre);
            choice = perception_ob > bias_d(heading);
            data_up = [stimuli_list(stimuli(sti_index)), choice, 1];
            sti_index = sti_index + 1;
            heading_index(i) = heading;
            data{heading}(end + 1, :) = data_up;

            iter = 5000;
            options = struct;
            options.poolMaxGap     = inf;
            options.poolMaxLength  = inf;
            options.poolxTol       = 0;

            if i >= 33
                 for hra = 1:3
                    num_deading = length(find(heading_index(1: i) == hra));
                    data_stan(:,:,hra) = poolData(data{hra}(1: num_deading, :), options);           
                 end

             monk_dat = struct('ndir',11, ...
                          'direc',squeeze(data_stan(:,1,:)),...
                          'choi',squeeze(data_stan(:,2,:)),...
                          'n', squeeze(data_stan(:,3,:)),...
                          'P0_1', mean_p1(num_all),...
                          'P0_2', mean_p2(num_all),...
                          'D0', bias_d0, ...
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
                fit = stan(params, 'init', struct('P_1', mean_p1(num_all), 'P_2', mean_p2(num_all), 'D',bias_d0, 'phi', [16, 16, 16]));
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
                b(i, 1) = mean_p1(num_all);
                b(i, 3) = mean_p2(num_all);
                b(i, 2) = bias_d0;
                c(i, 1, 1) = mean_p1(num_all) - std_p1;
                c(i, 3, 1) = mean_p2(num_all) - std_p2;
                c(i, 2, 1) = bias_d0-std_d;
                c(i, 1, 2) = mean_p1(num_all) + std_p1;
                c(i, 3, 2) = mean_p2(num_all) + std_p2;
                c(i, 2, 2) = bias_d0+std_d;
            end
                bias_track = b;
                bias_track(:, 2) = 0;
            %%
            if i == 1
                correct_flag = sign(sign(choice - 0.5) * (stimuli_list(stimuli(sti_index - 1)) - bias_track(1, heading)));
            else
                correct_flag = sign(sign(choice - 0.5) * (stimuli_list(stimuli(sti_index - 1)) - bias_track(i - 1, heading)));
            end

            if correct_flag == 0
                correct_flag = sign(rand - 0.5);
            end
            if correct_flag == -1
                correct_flag = 0;
            end
            % update decision criteria after each trial for independent
            if i == 1
                bias_updata(i, :) = bias_d;
            else
                for heading_all = 1: 3
                    if heading_all == heading
                        bias_updata(i, heading_all) = bias_updata(i-1, heading_all) - alpha * (correct_flag -  max(choice_prob, 1-choice_prob)) * sign(choice - 0.5);
                    else
                        bias_updata(i, heading_all) = bias_updata(i-1, heading_all);
                    end
                end
            end
        end
    %%
    % multisession model fitting
    bias_estimate{k, num_all} = b;
    bias_learn{k, num_all} = bias_updata;    
    prior_estimate{k, num_all} = c; 
    exp_all{k, num_all} = data_stan;
    data_stan_session(num_all, :,:,:) = data_stan;
    if num_all >= 2
    session_dat = struct('nse',num_all, ...
                      'ndir',11, ...
                      'direc1',squeeze(data_stan_session(:,:,1,1)),...
                      'choi1',squeeze(data_stan_session(:,:,2,1)),...                
                      'n1', squeeze(data_stan_session(:,:,3,1)),...
                      'direc2',squeeze(data_stan_session(:,:,1,2)),...
                      'choi2',squeeze(data_stan_session(:,:,2,2)),...                
                      'n2', squeeze(data_stan_session(:,:,3,2)),...
                      'direc3',squeeze(data_stan_session(:,:,1,3)),...
                      'choi3',squeeze(data_stan_session(:,:,2,3)),...                
                      'n3', squeeze(data_stan_session(:,:,3,3)),...
                      'lapse_alpha_1', 1,...
                      'lapse_beta_1', 10,...
                      'lapse_alpha_2', 1,...
                      'lapse_beta_2', 10,...
                      'phi_alpha', 8,...
                      'phi_beta', 0.5);

        params = struct('file','StanSession.stan', 'data', session_dat,'iter',iter,'chains',1);
        fit = stan(params, 'init', struct('P_1', mean_p1, 'P_2', mean_p2, 'D',zeros(num_all, 1), 'phi', 16 * ones(num_all, 3),...
            'P0_1', mean_p1(num_all), 'P0_2', mean_p2(num_all), 'tau_1', std_p1, 'tau_2', std_p2, 'tau_d', std_d));

        fit.verbose = true;

        waitfor(fit,'exit_value',0);
        para_dis = fit.extract();
        p_intial(num_all, 1) = quantile(para_dis.P0_1, 0.5);
        p_intial(num_all, 2) = quantile(para_dis.P0_2, 0.5);
        tau_intial(num_all, 1) = quantile(para_dis.tau_1, 0.5);
        tau_intial(num_all, 2) = quantile(para_dis.tau_2, 0.5);
        tau_intial(num_all, 3) = quantile(para_dis.tau_d, 0.5);
    elseif num_all == 1
        p_intial(num_all, 1) = b(end, 1);
        p_intial(num_all, 2) = b(end, 3);
        tau_intial(num_all, 1) = (c(end, 1, 2) - c(end, 1, 1))/2;
        tau_intial(num_all, 2) = (c(end, 3, 2) - c(end, 3, 1))/2;
        tau_intial(num_all, 3) = (c(end, 2, 2) - c(end, 2, 1))/2;
    end
    name_save = sprintf('save_learner_stan_track_50_rep1_11_share');
    save(name_save,'bias_estimate',...
         'bias_learn','prior_estimate', 'exp_all','p_intial','tau_intial');

    end
    options.sigmoidName = 'norm';  
    options.expType     = 'YesNo'; 
    options.estimateType = 'mean';
    for session = 1: 50
        for con_hea= 1: 3
            data_fit_sesssion = exp_all{k, session};
            data_fit = data_fit_sesssion(:,:,con_hea);
            result = psignifit(data_fit, options);
            bias_psig(session, con_hea, k) = result.Fit(1);
        end
    end
    bias_learn_all = [];

    for sess_num = 1: num_session
        bias_learn_all = [bias_learn_all; bias_learn{k, sess_num}];
    end

end
    bias_rate_learn{1} = bias_learn_all;
    bias_rate_psig{1} = bias_psig;

%% Data Generating for FigureS5CD: independent decision critera, learning rate = 0.0004

intial_p1 = -10; % prior mean for perceptual bias
intial_p2 = 20;
alpha = 0.0004; % learning rate
num_trial = 990;
num_session = 50;
std_in_d = 10; % prior std for decision bias
std_in_p = 5; % prior std for perceptual bias
bias_d0 = 0;
% changing ground truth percpetual bias
bias_all(1, 1:9) = -10;
bias_all(3, 1:9) = 20;
bias_all(1, 10: 41) = linspace(-10, -10/3, 32);
bias_all(3, 10: 41) = linspace(20, 20/3, 32);
bias_all(1, 42:50) = -10/3;
bias_all(3, 42:50) = 20/3;
bias_all(2, 1:50) = 0;
thre = 16;
for k = 1
    seed = sum(round(clock));

    for num_all = 1: num_session
        bias = bias_all(:, num_all);
        % update prior after each session
        if num_all == 1
            std_p1 = std_in_p(1);
            std_p2 = std_in_p(1);
            std_d = std_in_d(1);
            mean_p1(num_all) = intial_p1(1);
            mean_p2(num_all) = intial_p2(1);
        else
            std_p1 = tau_intial(num_all - 1, 1);
            std_p2 = tau_intial(num_all - 1, 2);
            std_d = tau_intial(num_all - 1, 3);
            mean_p1(num_all) = p_intial(num_all - 1, 1);
            mean_p2(num_all) = p_intial(num_all - 1, 2);        
        end
        % update decision critera
        % independent decision critera 
        if num_all == 1
            bias_dd = zeros(1, 3);
        else
            bias_dd = bias_learn{k, num_all - 1}(end, :);
        end

        rng(seed)

        heading_condition = [-10, 0, 10];
        direction_condition = -40: 8: 40;
        stimuli_list = [direction_condition, direction_condition, direction_condition];


        data = cell(3, 1);
        b = nan(num_trial, 3);
        c = nan(num_trial, 3, 2);
        heading_index = nan(num_trial, 1);
        data_stan = nan(length(direction_condition), 3, 3);
        bias_updata = nan(num_trial, 3); % independent
        for i = 1: num_trial

            if mod(i, 3 * length(direction_condition)) == 1
                sti_index = 1;
                stimuli = randperm(3 * length( direction_condition));
            end
            if stimuli(sti_index) > (2 * length( direction_condition))
                heading = 1;
            elseif stimuli(sti_index) > (1 * length( direction_condition))
                heading = 2;
            else 
                heading = 3;
            end

            if i == 1
                bias_d = bias_dd;
            else
                bias_d = bias_updata(i - 1, :);  % independent
            end

            perception_ob = normrnd(stimuli_list(stimuli(sti_index)) - bias(heading), thre);
            choice_prob = normcdf(perception_ob, bias_d(heading), thre);
            choice = perception_ob > bias_d(heading);
            data_up = [stimuli_list(stimuli(sti_index)), choice, 1];
            sti_index = sti_index + 1;
            heading_index(i) = heading;
            data{heading}(end + 1, :) = data_up;

            iter = 5000;
            options = struct;
            options.poolMaxGap     = inf;
            options.poolMaxLength  = inf;
            options.poolxTol       = 0;

            if i >= 33
                 for hra = 1:3
                    num_deading = length(find(heading_index(1: i) == hra));
                    data_stan(:,:,hra) = poolData(data{hra}(1: num_deading, :), options);           
                 end

             monk_dat = struct('ndir',11, ...
                          'direc',squeeze(data_stan(:,1,:)),...
                          'choi',squeeze(data_stan(:,2,:)),...
                          'n', squeeze(data_stan(:,3,:)),...
                          'P0_1', mean_p1(num_all),...
                          'P0_2', mean_p2(num_all),...
                          'D0', bias_d0, ...
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
                fit = stan(params, 'init', struct('P_1', mean_p1(num_all), 'P_2', mean_p2(num_all), 'D',bias_d0, 'phi', [16, 16, 16]));
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
                b(i, 1) = mean_p1(num_all);
                b(i, 3) = mean_p2(num_all);
                b(i, 2) = bias_d0;
                c(i, 1, 1) = mean_p1(num_all) - std_p1;
                c(i, 3, 1) = mean_p2(num_all) - std_p2;
                c(i, 2, 1) = bias_d0-std_d;
                c(i, 1, 2) = mean_p1(num_all) + std_p1;
                c(i, 3, 2) = mean_p2(num_all) + std_p2;
                c(i, 2, 2) = bias_d0+std_d;
            end
                bias_track = b;
                bias_track(:, 2) = 0;
            %%
            if i == 1
                correct_flag = sign(sign(choice - 0.5) * (stimuli_list(stimuli(sti_index - 1)) - bias_track(1, heading)));
            else
                correct_flag = sign(sign(choice - 0.5) * (stimuli_list(stimuli(sti_index - 1)) - bias_track(i - 1, heading)));
            end

            if correct_flag == 0
                correct_flag = sign(rand - 0.5);
            end
            if correct_flag == -1
                correct_flag = 0;
            end
            % update decision criteria after each trial for independent
            if i == 1
                bias_updata(i, :) = bias_d;
            else
                for heading_all = 1: 3
                    if heading_all == heading
                        bias_updata(i, heading_all) = bias_updata(i-1, heading_all) - alpha * (correct_flag -  max(choice_prob, 1-choice_prob)) * sign(choice - 0.5);
                    else
                        bias_updata(i, heading_all) = bias_updata(i-1, heading_all);
                    end
                end
            end
        end
    %%
    % multisession model fitting
    bias_estimate{k, num_all} = b;
    bias_learn{k, num_all} = bias_updata;    
    prior_estimate{k, num_all} = c; 
    exp_all{k, num_all} = data_stan;
    data_stan_session(num_all, :,:,:) = data_stan;
    if num_all >= 2
    session_dat = struct('nse',num_all, ...
                      'ndir',11, ...
                      'direc1',squeeze(data_stan_session(:,:,1,1)),...
                      'choi1',squeeze(data_stan_session(:,:,2,1)),...                
                      'n1', squeeze(data_stan_session(:,:,3,1)),...
                      'direc2',squeeze(data_stan_session(:,:,1,2)),...
                      'choi2',squeeze(data_stan_session(:,:,2,2)),...                
                      'n2', squeeze(data_stan_session(:,:,3,2)),...
                      'direc3',squeeze(data_stan_session(:,:,1,3)),...
                      'choi3',squeeze(data_stan_session(:,:,2,3)),...                
                      'n3', squeeze(data_stan_session(:,:,3,3)),...
                      'lapse_alpha_1', 1,...
                      'lapse_beta_1', 10,...
                      'lapse_alpha_2', 1,...
                      'lapse_beta_2', 10,...
                      'phi_alpha', 8,...
                      'phi_beta', 0.5);

        params = struct('file','StanSession.stan', 'data', session_dat,'iter',iter,'chains',1);
        fit = stan(params, 'init', struct('P_1', mean_p1, 'P_2', mean_p2, 'D',zeros(num_all, 1), 'phi', 16 * ones(num_all, 3),...
            'P0_1', mean_p1(num_all), 'P0_2', mean_p2(num_all), 'tau_1', std_p1, 'tau_2', std_p2, 'tau_d', std_d));

        fit.verbose = true;

        waitfor(fit,'exit_value',0);
        para_dis = fit.extract();
        p_intial(num_all, 1) = quantile(para_dis.P0_1, 0.5);
        p_intial(num_all, 2) = quantile(para_dis.P0_2, 0.5);
        tau_intial(num_all, 1) = quantile(para_dis.tau_1, 0.5);
        tau_intial(num_all, 2) = quantile(para_dis.tau_2, 0.5);
        tau_intial(num_all, 3) = quantile(para_dis.tau_d, 0.5);
    elseif num_all == 1
        p_intial(num_all, 1) = b(end, 1);
        p_intial(num_all, 2) = b(end, 3);
        tau_intial(num_all, 1) = (c(end, 1, 2) - c(end, 1, 1))/2;
        tau_intial(num_all, 2) = (c(end, 3, 2) - c(end, 3, 1))/2;
        tau_intial(num_all, 3) = (c(end, 2, 2) - c(end, 2, 1))/2;
    end
    name_save = sprintf('save_learner_stan_track_50_rep1_11_share');
    save(name_save,'bias_estimate',...
         'bias_learn','prior_estimate', 'exp_all','p_intial','tau_intial');

    end
    options.sigmoidName = 'norm';  
    options.expType     = 'YesNo'; 
    options.estimateType = 'mean';
    for session = 1: 50
        for con_hea= 1: 3
            data_fit_sesssion = exp_all{k, session};
            data_fit = data_fit_sesssion(:,:,con_hea);
            result = psignifit(data_fit, options);
            bias_psig(session, con_hea, k) = result.Fit(1);
        end
    end
    bias_learn_all = [];
    for sess_num = 1: num_session
        bias_learn_all = [bias_learn_all; bias_learn{k, sess_num}];
    end

end
    bias_rate_learn{2} = bias_learn_all;
    bias_rate_psig{2} = bias_psig;
%% Data Generating for FigureS5EF: shared decision critera, learning rate = 0.04

intial_p1 = -10; % prior mean for perceptual bias
intial_p2 = 20;
alpha = 0.04; % learning rate
num_trial = 990;
num_session = 50;
std_in_d = 10; % prior std for decision bias
std_in_p = 5; % prior std for perceptual bias
bias_d0 = 0;
% changing ground truth percpetual bias
bias_all(1, 1:9) = -10;
bias_all(3, 1:9) = 20;
bias_all(1, 10: 41) = linspace(-10, -10/3, 32);
bias_all(3, 10: 41) = linspace(20, 20/3, 32);
bias_all(1, 42:50) = -10/3;
bias_all(3, 42:50) = 20/3;
bias_all(2, 1:50) = 0;
thre = 16;

for k = 1
    seed = sum(round(clock));

    for num_all = 1: num_session
        bias = bias_all(:, num_all);
        % update prior after each session
        if num_all == 1
            std_p1 = std_in_p(1);
            std_p2 = std_in_p(1);
            std_d = std_in_d(1);
            mean_p1(num_all) = intial_p1(1);
            mean_p2(num_all) = intial_p2(1);
        else
            std_p1 = tau_intial(num_all - 1, 1);
            std_p2 = tau_intial(num_all - 1, 2);
            std_d = tau_intial(num_all - 1, 3);
            mean_p1(num_all) = p_intial(num_all - 1, 1);
            mean_p2(num_all) = p_intial(num_all - 1, 2);        
        end
        % update decision critera
        % shared decision critera 
        if num_all == 1
            bias_dd = 0;
        else
            bias_dd = bias_learn{k, num_all - 1}(end);
        end

        rng(seed)

        heading_condition = [-10, 0, 10];
        direction_condition = -40: 8: 40;
        stimuli_list = [direction_condition, direction_condition, direction_condition];


        data = cell(3, 1);
        b = nan(num_trial, 3);
        c = nan(num_trial, 3, 2);
        heading_index = nan(num_trial, 1);
        data_stan = nan(length(direction_condition), 3, 3);
        bias_updata = nan(num_trial, 1); % shared
        for i = 1: num_trial

            if mod(i, 3 * length(direction_condition)) == 1
                sti_index = 1;
                stimuli = randperm(3 * length( direction_condition));
            end
            if stimuli(sti_index) > (2 * length( direction_condition))
                heading = 1;
            elseif stimuli(sti_index) > (1 * length( direction_condition))
                heading = 2;
            else 
                heading = 3;
            end

            if i == 1
                bias_d = bias_dd;
            else
                bias_d = bias_updata(i - 1); % shared
            end

            perception_ob = normrnd(stimuli_list(stimuli(sti_index)) - bias(heading), thre);
            choice_prob = normcdf(perception_ob, bias_d, thre);
            choice = perception_ob > bias_d;
            data_up = [stimuli_list(stimuli(sti_index)), choice, 1];
            sti_index = sti_index + 1;
            heading_index(i) = heading;
            data{heading}(end + 1, :) = data_up;

            iter = 5000;
            options = struct;
            options.poolMaxGap     = inf;
            options.poolMaxLength  = inf;
            options.poolxTol       = 0;

            if i >= 33
                 for hra = 1:3
                    num_deading = length(find(heading_index(1: i) == hra));
                    data_stan(:,:,hra) = poolData(data{hra}(1: num_deading, :), options);           
                 end

             monk_dat = struct('ndir',11, ...
                          'direc',squeeze(data_stan(:,1,:)),...
                          'choi',squeeze(data_stan(:,2,:)),...
                          'n', squeeze(data_stan(:,3,:)),...
                          'P0_1', mean_p1(num_all),...
                          'P0_2', mean_p2(num_all),...
                          'D0', bias_d0, ...
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
                fit = stan(params, 'init', struct('P_1', mean_p1(num_all), 'P_2', mean_p2(num_all), 'D',bias_d0, 'phi', [16, 16, 16]));
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
                b(i, 1) = mean_p1(num_all);
                b(i, 3) = mean_p2(num_all);
                b(i, 2) = bias_d0;
                c(i, 1, 1) = mean_p1(num_all) - std_p1;
                c(i, 3, 1) = mean_p2(num_all) - std_p2;
                c(i, 2, 1) = bias_d0-std_d;
                c(i, 1, 2) = mean_p1(num_all) + std_p1;
                c(i, 3, 2) = mean_p2(num_all) + std_p2;
                c(i, 2, 2) = bias_d0+std_d;
            end
                bias_track = b;
                bias_track(:, 2) = 0;
            %%
            if i == 1
                correct_flag = sign(sign(choice - 0.5) * (stimuli_list(stimuli(sti_index - 1)) - bias_track(1, heading)));
            else
                correct_flag = sign(sign(choice - 0.5) * (stimuli_list(stimuli(sti_index - 1)) - bias_track(i - 1, heading)));
            end

            if correct_flag == 0
                correct_flag = sign(rand - 0.5);
            end
            if correct_flag == -1
                correct_flag = 0;
            end
            % update decision criteria after each trial for shared
            if i == 1
                bias_updata(i) = bias_d;
            else
                bias_updata(i) = bias_updata(i-1) - alpha * (correct_flag -  max(choice_prob, 1-choice_prob)) * sign(choice - 0.5);
            end
            
        end
    %%
    % multisession model fitting
    bias_estimate{k, num_all} = b;
    bias_learn{k, num_all} = bias_updata;    
    prior_estimate{k, num_all} = c; 
    exp_all{k, num_all} = data_stan;
    data_stan_session(num_all, :,:,:) = data_stan;
    if num_all >= 2
    session_dat = struct('nse',num_all, ...
                      'ndir',11, ...
                      'direc1',squeeze(data_stan_session(:,:,1,1)),...
                      'choi1',squeeze(data_stan_session(:,:,2,1)),...                
                      'n1', squeeze(data_stan_session(:,:,3,1)),...
                      'direc2',squeeze(data_stan_session(:,:,1,2)),...
                      'choi2',squeeze(data_stan_session(:,:,2,2)),...                
                      'n2', squeeze(data_stan_session(:,:,3,2)),...
                      'direc3',squeeze(data_stan_session(:,:,1,3)),...
                      'choi3',squeeze(data_stan_session(:,:,2,3)),...                
                      'n3', squeeze(data_stan_session(:,:,3,3)),...
                      'lapse_alpha_1', 1,...
                      'lapse_beta_1', 10,...
                      'lapse_alpha_2', 1,...
                      'lapse_beta_2', 10,...
                      'phi_alpha', 8,...
                      'phi_beta', 0.5);

        params = struct('file','StanSession.stan', 'data', session_dat,'iter',iter,'chains',1);
        fit = stan(params, 'init', struct('P_1', mean_p1, 'P_2', mean_p2, 'D',zeros(num_all, 1), 'phi', 16 * ones(num_all, 3),...
            'P0_1', mean_p1(num_all), 'P0_2', mean_p2(num_all), 'tau_1', std_p1, 'tau_2', std_p2, 'tau_d', std_d));

        fit.verbose = true;

        waitfor(fit,'exit_value',0);
        para_dis = fit.extract();
        p_intial(num_all, 1) = quantile(para_dis.P0_1, 0.5);
        p_intial(num_all, 2) = quantile(para_dis.P0_2, 0.5);
        tau_intial(num_all, 1) = quantile(para_dis.tau_1, 0.5);
        tau_intial(num_all, 2) = quantile(para_dis.tau_2, 0.5);
        tau_intial(num_all, 3) = quantile(para_dis.tau_d, 0.5);
    elseif num_all == 1
        p_intial(num_all, 1) = b(end, 1);
        p_intial(num_all, 2) = b(end, 3);
        tau_intial(num_all, 1) = (c(end, 1, 2) - c(end, 1, 1))/2;
        tau_intial(num_all, 2) = (c(end, 3, 2) - c(end, 3, 1))/2;
        tau_intial(num_all, 3) = (c(end, 2, 2) - c(end, 2, 1))/2;
    end
    name_save = sprintf('save_learner_stan_track_50_rep1_11_share');
    save(name_save,'bias_estimate',...
         'bias_learn','prior_estimate', 'exp_all','p_intial','tau_intial');

    end
    options.sigmoidName = 'norm';  
    options.expType     = 'YesNo'; 
    options.estimateType = 'mean';
    for session = 1: 50
        for con_hea= 1: 3
            data_fit_sesssion = exp_all{k, session};
            data_fit = data_fit_sesssion(:,:,con_hea);
            result = psignifit(data_fit, options);
            bias_psig(session, con_hea, k) = result.Fit(1);
        end
    end
    
    bias_learn_all = [];
    for sess_num = 1: num_session
        bias_learn_all = [bias_learn_all; bias_learn{k, sess_num}];
    end

end
    bias_rate_learn{3} = bias_learn_all;
    bias_rate_psig{3} = bias_psig;
%% Data Generating for FigureS5GH: shared decision critera, learning rate = 4

intial_p1 = -10; % prior mean for perceptual bias
intial_p2 = 20;
alpha = 4; % learning rate
num_trial = 990;
num_session = 50;
std_in_d = 10; % prior std for decision bias
std_in_p = 5; % prior std for perceptual bias
bias_d0 = 0;
% changing ground truth percpetual bias
bias_all(1, 1:9) = -10;
bias_all(3, 1:9) = 20;
bias_all(1, 10: 41) = linspace(-10, -10/3, 32);
bias_all(3, 10: 41) = linspace(20, 20/3, 32);
bias_all(1, 42:50) = -10/3;
bias_all(3, 42:50) = 20/3;
bias_all(2, 1:50) = 0;
thre = 16;

for k = 1
    seed = sum(round(clock));

    for num_all = 1: num_session
        bias = bias_all(:, num_all);
        % update prior after each session
        if num_all == 1
            std_p1 = std_in_p(1);
            std_p2 = std_in_p(1);
            std_d = std_in_d(1);
            mean_p1(num_all) = intial_p1(1);
            mean_p2(num_all) = intial_p2(1);
        else
            std_p1 = tau_intial(num_all - 1, 1);
            std_p2 = tau_intial(num_all - 1, 2);
            std_d = tau_intial(num_all - 1, 3);
            mean_p1(num_all) = p_intial(num_all - 1, 1);
            mean_p2(num_all) = p_intial(num_all - 1, 2);        
        end
        % update decision critera
        % shared decision critera 
        if num_all == 1
            bias_dd = 0;
        else
            bias_dd = bias_learn{k, num_all - 1}(end);
        end

        rng(seed)

        heading_condition = [-10, 0, 10];
        direction_condition = -40: 8: 40;
        stimuli_list = [direction_condition, direction_condition, direction_condition];


        data = cell(3, 1);
        b = nan(num_trial, 3);
        c = nan(num_trial, 3, 2);
        heading_index = nan(num_trial, 1);
        data_stan = nan(length(direction_condition), 3, 3);
        bias_updata = nan(num_trial, 1); % shared
        for i = 1: num_trial

            if mod(i, 3 * length(direction_condition)) == 1
                sti_index = 1;
                stimuli = randperm(3 * length( direction_condition));
            end
            if stimuli(sti_index) > (2 * length( direction_condition))
                heading = 1;
            elseif stimuli(sti_index) > (1 * length( direction_condition))
                heading = 2;
            else 
                heading = 3;
            end

            if i == 1
                bias_d = bias_dd;
            else
                bias_d = bias_updata(i - 1); % shared
            end

            perception_ob = normrnd(stimuli_list(stimuli(sti_index)) - bias(heading), thre);
            choice_prob = normcdf(perception_ob, bias_d, thre);
            choice = perception_ob > bias_d;
            data_up = [stimuli_list(stimuli(sti_index)), choice, 1];
            sti_index = sti_index + 1;
            heading_index(i) = heading;
            data{heading}(end + 1, :) = data_up;

            iter = 5000;
            options = struct;
            options.poolMaxGap     = inf;
            options.poolMaxLength  = inf;
            options.poolxTol       = 0;

            if i >= 33
                 for hra = 1:3
                    num_deading = length(find(heading_index(1: i) == hra));
                    data_stan(:,:,hra) = poolData(data{hra}(1: num_deading, :), options);           
                 end

             monk_dat = struct('ndir',11, ...
                          'direc',squeeze(data_stan(:,1,:)),...
                          'choi',squeeze(data_stan(:,2,:)),...
                          'n', squeeze(data_stan(:,3,:)),...
                          'P0_1', mean_p1(num_all),...
                          'P0_2', mean_p2(num_all),...
                          'D0', bias_d0, ...
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
                fit = stan(params, 'init', struct('P_1', mean_p1(num_all), 'P_2', mean_p2(num_all), 'D',bias_d0, 'phi', [16, 16, 16]));
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
                b(i, 1) = mean_p1(num_all);
                b(i, 3) = mean_p2(num_all);
                b(i, 2) = bias_d0;
                c(i, 1, 1) = mean_p1(num_all) - std_p1;
                c(i, 3, 1) = mean_p2(num_all) - std_p2;
                c(i, 2, 1) = bias_d0-std_d;
                c(i, 1, 2) = mean_p1(num_all) + std_p1;
                c(i, 3, 2) = mean_p2(num_all) + std_p2;
                c(i, 2, 2) = bias_d0+std_d;
            end
                bias_track = b;
                bias_track(:, 2) = 0;
            %%
            if i == 1
                correct_flag = sign(sign(choice - 0.5) * (stimuli_list(stimuli(sti_index - 1)) - bias_track(1, heading)));
            else
                correct_flag = sign(sign(choice - 0.5) * (stimuli_list(stimuli(sti_index - 1)) - bias_track(i - 1, heading)));
            end

            if correct_flag == 0
                correct_flag = sign(rand - 0.5);
            end
            if correct_flag == -1
                correct_flag = 0;
            end
            % update decision criteria after each trial for shared
            if i == 1
                bias_updata(i) = bias_d;
            else
                bias_updata(i) = bias_updata(i-1) - alpha * (correct_flag -  max(choice_prob, 1-choice_prob)) * sign(choice - 0.5);
            end
            
        end
    %%
    % multisession model fitting
    bias_estimate{k, num_all} = b;
    bias_learn{k, num_all} = bias_updata;    
    prior_estimate{k, num_all} = c; 
    exp_all{k, num_all} = data_stan;
    data_stan_session(num_all, :,:,:) = data_stan;
    if num_all >= 2
    session_dat = struct('nse',num_all, ...
                      'ndir',11, ...
                      'direc1',squeeze(data_stan_session(:,:,1,1)),...
                      'choi1',squeeze(data_stan_session(:,:,2,1)),...                
                      'n1', squeeze(data_stan_session(:,:,3,1)),...
                      'direc2',squeeze(data_stan_session(:,:,1,2)),...
                      'choi2',squeeze(data_stan_session(:,:,2,2)),...                
                      'n2', squeeze(data_stan_session(:,:,3,2)),...
                      'direc3',squeeze(data_stan_session(:,:,1,3)),...
                      'choi3',squeeze(data_stan_session(:,:,2,3)),...                
                      'n3', squeeze(data_stan_session(:,:,3,3)),...
                      'lapse_alpha_1', 1,...
                      'lapse_beta_1', 10,...
                      'lapse_alpha_2', 1,...
                      'lapse_beta_2', 10,...
                      'phi_alpha', 8,...
                      'phi_beta', 0.5);

        params = struct('file','StanSession.stan', 'data', session_dat,'iter',iter,'chains',1);
        fit = stan(params, 'init', struct('P_1', mean_p1, 'P_2', mean_p2, 'D',zeros(num_all, 1), 'phi', 16 * ones(num_all, 3),...
            'P0_1', mean_p1(num_all), 'P0_2', mean_p2(num_all), 'tau_1', std_p1, 'tau_2', std_p2, 'tau_d', std_d));

        fit.verbose = true;

        waitfor(fit,'exit_value',0);
        para_dis = fit.extract();
        p_intial(num_all, 1) = quantile(para_dis.P0_1, 0.5);
        p_intial(num_all, 2) = quantile(para_dis.P0_2, 0.5);
        tau_intial(num_all, 1) = quantile(para_dis.tau_1, 0.5);
        tau_intial(num_all, 2) = quantile(para_dis.tau_2, 0.5);
        tau_intial(num_all, 3) = quantile(para_dis.tau_d, 0.5);
    elseif num_all == 1
        p_intial(num_all, 1) = b(end, 1);
        p_intial(num_all, 2) = b(end, 3);
        tau_intial(num_all, 1) = (c(end, 1, 2) - c(end, 1, 1))/2;
        tau_intial(num_all, 2) = (c(end, 3, 2) - c(end, 3, 1))/2;
        tau_intial(num_all, 3) = (c(end, 2, 2) - c(end, 2, 1))/2;
    end
    name_save = sprintf('save_learner_stan_track_50_rep1_11_share');
    save(name_save,'bias_estimate',...
         'bias_learn','prior_estimate', 'exp_all','p_intial','tau_intial');

    end
    options.sigmoidName = 'norm';  
    options.expType     = 'YesNo'; 
    options.estimateType = 'mean';
    for session = 1: 50
        for con_hea= 1: 3
            data_fit_sesssion = exp_all{k, session};
            data_fit = data_fit_sesssion(:,:,con_hea);
            result = psignifit(data_fit, options);
            bias_psig(session, con_hea, k) = result.Fit(1);
        end
    end
    
    bias_learn_all = [];
    for sess_num = 1: num_session
        bias_learn_all = [bias_learn_all; bias_learn{k, sess_num}];
    end

end
    bias_rate_learn{4} = bias_learn_all;
    bias_rate_psig{4} = bias_psig;