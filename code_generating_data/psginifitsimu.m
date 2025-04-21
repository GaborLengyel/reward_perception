
%%
function [bias_raw, conf_raw] = psginifitsimu(thre, bias, mean_fit, pr, num_trial, cons) 
rng(sum(round(clock)));

bias_raw = cell(3, length(pr(1, :)) + 1);
result_list = cell(3, length(pr(1, :)) + 1);

bias_real = cell(3, length(pr(1, :)) + 1);
conf_real = cell(3, length(pr(1, :)) + 1);
num_upd = cell(3, length(pr(1, :)) + 1);
clevel = 0.68; % credic interval
for pp = 1: length(pr(1,:))



    direction_condition = -40: 8: 40;

    stimuli_list = [direction_condition, direction_condition, direction_condition];

    trial_index = 1;
    fit_flag = [0 0 0];
    data = cell(3, 1);
    post = cell(3, 1);
    if pp <= length(pr(1,:)) 
        for kk = 1: 3
            prior_input{kk, pp}.x = min(direction_condition): max(direction_condition);
            prior_input{kk, pp}.y = normpdf(prior_input{kk, pp}.x, mean_fit(kk), pr(kk, pp));
        end
    end
    pr_real = [sqrt(pr(1, pp)^2 - pr(2, pp)^2), pr(2, pp), sqrt(pr(3, pp)^2 - pr(2, pp)^2)];
    % constant decision bias or not
    if cons == 0
        bias_d0 = 10 * sin((1: num_trial)/2000 * pi);
    elseif cons == 1
        bias_d0 = 10 * ones(1, num_trial);
    end


    for i = 1: num_trial
        if mod(i, 3 * length(direction_condition)) == 1
            sti_index = 1;
            stimuli = randperm(3 * length( direction_condition));
        end
        if stimuli(sti_index) > (2 * length( direction_condition))
            bias_trial = bias(1) + bias_d0(i);
            heading = 1;
        elseif stimuli(sti_index) > (1 * length( direction_condition))
            bias_trial = bias_d0(i);
            heading = 2;
        else 
            bias_trial = bias(3) + bias_d0(i);
            heading = 3;
        end
            choice_prob = normcdf(stimuli_list(stimuli(sti_index)), bias_trial, thre);
            choice = binornd(1, choice_prob);
            data_up = [stimuli_list(stimuli(sti_index)), choice, 1];
            sti_index = sti_index + 1;
            data{heading}(end + 1, :) = data_up;
            if length(data{heading}) >= length(direction_condition)
                if pp == length(pr(1, :)) + 1
                    [bias_trail, porst_trial, result] = calcu_bias(data{heading});
                else
                    [bias_trail, porst_trial, result] = calcu_bias(data{heading}, prior_input{heading, pp});
                end
                for heading_index = 1: 3
                    if heading_index == heading
                        bias_raw{heading_index, pp}(trial_index) = bias_trail;
                        conf_raw{heading_index, pp, trial_index} = porst_trial;
                    else
                        bias_raw{heading_index, pp}(trial_index) = bias_raw{heading_index, pp}(trial_index - 1);
                        conf_raw{heading_index, pp, trial_index} = conf_raw{heading_index, pp, trial_index - 1};
                    end
                end
                result_list{heading, pp} = result;
                fit_flag(heading) = 1;
                post{heading} = prior_up(porst_trial.x, porst_trial.y, min(direction_condition): max(direction_condition));
                if sum(fit_flag) >= 3
                    heading_conv_range = min(direction_condition): 0.5: max(direction_condition);
                    p_distri = conv(post{1}, fliplr(post{3})); % distribution of perceptual bias1
                    bias_p = sum(p_distri .* heading_conv_range)/sum(p_distri); % mean of perceptual bias
                    bias_real{1, pp}(end + 1) = bias_p;
                    conf_real{1, pp}(end + 1, :) = ...
                        CIdist(heading_conv_range, p_distri/sum(p_distri), clevel); % CI for perceptual bias
                    bias_real{3, pp}(end + 1) = -bias_p;
                    conf_real{3, pp}(end + 1, :) = [-conf_real{1, pp}(end, 2), -conf_real{1, pp}(end, 1)];
                    d_distri = conv(post{1}, post{3}); % distribution of decision bias
                    bias_d = sum(d_distri .* heading_conv_range)/sum(d_distri); % mean of decision bias1
                    bias_real{2, pp}(end + 1) = bias_d;
                    conf_real{2, pp}(end + 1, :) = ...
                        CIdist(heading_conv_range, d_distri/sum(d_distri), clevel); % CI for decision bias
                    num_upd{1, pp}(end + 1) = trial_index; 
                    num_upd{2, pp}(end + 1) = trial_index; 
                    num_upd{3, pp}(end + 1) = trial_index;

                else
                    bias_real{1, pp}(trial_index) = mean_fit(1);
                    bias_real{2, pp}(trial_index) = mean_fit(2);
                    bias_real{3, pp}(trial_index) = mean_fit(3);
                    conf_real{1, pp}(trial_index, 1: 2) = [-pr_real(1) + mean_fit(1), pr_real(1) + mean_fit(1)];
                    conf_real{2, pp}(trial_index, 1: 2) = [-pr_real(2) + mean_fit(2), pr_real(2) + mean_fit(2)];
                    conf_real{3, pp}(trial_index, 1: 2) = [-pr_real(3) + mean_fit(3), pr_real(3) + mean_fit(3)];
                    num_upd{1, pp}(trial_index) = trial_index; 
                    num_upd{2, pp}(trial_index) = trial_index; 
                    num_upd{3, pp}(trial_index) = trial_index;
                end
            else
                    bias_real{1, pp}(trial_index) = mean_fit(1);
                    bias_real{2, pp}(trial_index) = mean_fit(2);
                    bias_real{3, pp}(trial_index) = mean_fit(3);
                    for fl = 1: 3
                        if fit_flag(fl) == 1
                            bias_raw{fl, pp}(trial_index) = bias_raw{fl, pp}(trial_index-1);
                        elseif fit_flag(fl) == 0
                            bias_raw{fl, pp}(trial_index) = mean_fit(fl);
                        end
                    end                
                    conf_raw{1, pp, trial_index} = [];
                    conf_raw{2, pp, trial_index} = [];
                    conf_raw{3, pp, trial_index} = [];
                    conf_real{1, pp}(trial_index, 1: 2) = [-pr_real(1) + mean_fit(1), pr_real(1) + mean_fit(1)];
                    conf_real{2, pp}(trial_index, 1: 2) = [-pr_real(2) + mean_fit(2), pr_real(2) + mean_fit(2)];
                    conf_real{3, pp}(trial_index, 1: 2) = [-pr_real(3) + mean_fit(3), pr_real(3) + mean_fit(3)];
                    num_upd{1, pp}(trial_index) = trial_index; 
                    num_upd{2, pp}(trial_index) = trial_index; 
                    num_upd{3, pp}(trial_index) = trial_index;
            end

            trial_index = trial_index + 1;
    end

end


%     save(ss, 'porst')



