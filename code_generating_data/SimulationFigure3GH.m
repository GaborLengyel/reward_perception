clear
clc
reward_strategy = 'track';%'all' 'track' 'veridical' 'random' 'none'
intial_p1 = -10;
intial_p2 = 20;
alpha = 0.04;
num_trial = 990;
num_session = 50;
std_in_d = 10;
std_in_p = 5;
bias_d0 = 0;
bias = [-10, 0, 20];
thre = 16;


for k = 1:10
    seed = sum(round(clock));
    for num_all = 1: num_session

        if num_all == 1
            std_p1 = std_in_p(1);
            std_p2 = std_in_p(1);
            std_d = std_in_d(1);
        else
            std_p1 = abs(prior_estimate{k, num_all - 1}(end, 1, 2) - prior_estimate{k, num_all - 1}(end, 1, 1))/2;
            std_p2 = abs(prior_estimate{k, num_all - 1}(end, 3, 2) - prior_estimate{k, num_all - 1}(end, 3, 1))/2;
            std_d = abs(prior_estimate{k, num_all - 1}(end, 2, 2) - prior_estimate{k, num_all - 1}(end, 2, 1))/2;
        end

        if num_all == 1
            bias_dd = zeros(1, 3);
        else
            bias_dd = bias_learn{k, num_all - 1}(end, :);
        end

        rng(seed)

        heading_condition = [-10, 0, 10];
        direction_condition = -40: 8: 40;
        stimuli_list = [direction_condition, direction_condition, direction_condition];

        options = struct; 
        options.poolMaxGap     = inf;            
        options.poolMaxLength  = inf;            
        options.poolxTol       = 0;
        data = cell(3, 1);
        b = nan(num_trial, 3);
        c = nan(num_trial, 3, 2);
        heading_index = nan(num_trial, 1);
        data_stan = nan(length(direction_condition), 3, 3);
        bias_updata = nan(num_trial, 3);
        choice = nan(num_trial, 1);
        correct_flag = nan(num_trial, 1);
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
                bias_d = bias_updata(i - 1, :);
            end
                perception_ob = normrnd(stimuli_list(stimuli(sti_index)) - bias(heading), thre);
                choice_prob = normcdf(perception_ob, bias_d(heading), thre);
                choice(i) = perception_ob > bias_d(heading);
                data_up = [stimuli_list(stimuli(sti_index)), choice(i), 1];
                sti_index = sti_index + 1;
                heading_index(i) = heading;
                data{heading}(end + 1, :) = data_up;


            iter = 5000;
            options = struct;
            options.poolMaxGap     = inf;
            options.poolMaxLength  = inf;
            options.poolxTol       = 0;

            i
            if i >= 33
                 for hra = 1:3
                    num_deading = length(find(heading_index(1: i) == hra));
                    data_stan(:,:,hra) = poolData(data{hra}(1: num_deading, :), options);           
                 end

                 monk_dat = struct('ndir',11, ...
                              'direc',squeeze(data_stan(:,1,:)),...
                              'choi',squeeze(data_stan(:,2,:)),...
                              'n', squeeze(data_stan(:,3,:)),...
                              'P0_1', intial_p1(1),...
                              'P0_2', intial_p2(1),...
                              'D0', 0, ...
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
                fit = stan(params, 'init', struct('P_1', intial_p1(1), 'P_2', intial_p2(1), 'D',0, 'phi', [16, 16, 16]));
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
                b(i, 1) = intial_p1(1);
                b(i, 3) = intial_p2(1);
                b(i, 2) = bias_d0;
                c(i, 1, 1) = intial_p1(1) - std_p1;
                c(i, 3, 1) = intial_p2(1) - std_p2;
                c(i, 2, 1) = bias_d0-std_d;
                c(i, 1, 2) = intial_p1(1) + std_p1;
                c(i, 3, 2) = intial_p2(1) + std_p2;
                c(i, 2, 2) = bias_d0+std_d;
            end
            switch reward_strategy
                case 'track'
                    bias_track = b;
                    bias_track(:, 2) = 0;
                case 'veridical'
                    bias_track = zeros(size(b));
                case 'random'
                    bias_track = zeros(size(b));
                    if heading == 1
                        if stimuli_list(stimuli(sti_index - 1)) >= intial_p1 && stimuli_list(stimuli(sti_index - 1))<0
                            bias_track(i - 1, heading) = stimuli_list(stimuli(sti_index - 1));

                        end
                    elseif heading == 3
                        if stimuli_list(stimuli(sti_index - 1)) <= intial_p2 && stimuli_list(stimuli(sti_index - 1))>0 
                            bias_track(i - 1, heading) = stimuli_list(stimuli(sti_index - 1));
                        end
                    end
                case 'all'
                    bias_track = zeros(size(b));
                    if heading == 1
                        if stimuli_list(stimuli(sti_index - 1)) >= intial_p1 && stimuli_list(stimuli(sti_index - 1))<0
                            bias_track(i - 1, heading) = stimuli_list(stimuli(sti_index - 1)) - sign(choice(i) - 0.5);

                        end
                    elseif heading == 3
                        if stimuli_list(stimuli(sti_index - 1)) <= intial_p2 && stimuli_list(stimuli(sti_index - 1))>0 
                            bias_track(i - 1, heading) = stimuli_list(stimuli(sti_index - 1)) - sign(choice(i) - 0.5);
                        end
                    end
                case 'none'
                    bias_track = zeros(size(b));
                    if heading == 1
                        if stimuli_list(stimuli(sti_index - 1)) >= intial_p1 && stimuli_list(stimuli(sti_index - 1))<0
                            bias_track(i - 1, heading) = stimuli_list(stimuli(sti_index - 1)) + sign(choice(i) - 0.5);

                        end
                    elseif heading == 3
                        if stimuli_list(stimuli(sti_index - 1)) <= intial_p2 && stimuli_list(stimuli(sti_index - 1))>0 
                            bias_track(i - 1, heading) = stimuli_list(stimuli(sti_index - 1)) + sign(choice(i) - 0.5);
                        end
                    end
            end
            %%
            if i == 1
                correct_flag(i) = sign(sign(choice(i) - 0.5) * (stimuli_list(stimuli(sti_index - 1)) - bias_track(1, heading)));
            else
                correct_flag(i) = sign(sign(choice(i) - 0.5) * (stimuli_list(stimuli(sti_index - 1)) - bias_track(i - 1, heading)));
            end

            if correct_flag(i) == 0
                correct_flag(i) = sign(rand - 0.5);
            end
            if correct_flag(i) == -1
                correct_flag(i) = 0;
            end
            if i == 1
                bias_updata(i, :) = bias_d;
            else
                for heading_all = 1: 3
                    if heading_all == heading
                        bias_updata(i, heading_all) = bias_updata(i-1, heading_all) - alpha * (correct_flag(i) -  max(choice_prob, 1-choice_prob)) * sign(choice(i) - 0.5);
                    else
                        bias_updata(i, heading_all) = bias_updata(i-1, heading_all);
                    end
                end
            end
        end
        bias_estimate{k, num_all} = b;
        bias_learn{k, num_all} = bias_updata;    
        prior_estimate{k, num_all} = c; 
        exp_all{k, num_all} = data_stan;

        name_save = sprintf('save_learner_stan_track');
        save(name_save,'bias_estimate',...
             'bias_learn','prior_estimate', 'exp_all');

    end
    bias_learn_all = [];
    for session = 1: num_session
        bias_learn_all = [bias_learn_all; bias_learn{k, session}];
        
        for con_hea= 1: 3
            data_fit_sesssion = exp_all{k, session};
            data_fit = data_fit_sesssion(:,:,con_hea);
            result = psignifit(data_fit, options);
            bias_psig(session, con_hea, k) = result.Fit(1);
        end
    end
    bias_savve(:, :, k) = bias_learn_all; 
end
bias_method_learn{7} = bias_savve;
bias_method_psig{7} = bias_psig;
