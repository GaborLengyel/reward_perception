%% Data Generating for Figure3 & FigureS3
clear
clc
alpha = 0.04; % learning rate
num_trial = 50 * 990; % session times trials per session 
thre = 16; % psychometrical threshold of the agent
bias_p = [-10, 0, 20]; % ground truth perceptual bias
bias_dd = 0; % initial decision criteria
direction_condition = -40: 8: 40; % stimulus direction
stimuli_list = [direction_condition, direction_condition, direction_condition]; % stimulus direction for 3 conditions
% loop through different reward methods.
% 1: reward vericially (Fig.4 AB); 
% 2: random reward (Fig.4 CD); 
% 3: reward subjectivly (Fig.4 EF);
% 4: always reward (Fig.S3 AB); 
% 5: never reward (Fig.S3 CD).
% 6: balanced stimulus (Fig.S3 EF).
for method = 1: 6     
    bias_psig = nan(50, 3, 10); % session perceptual bias fitted by psignifit
    bias_savve = nan(num_trial, 3, 10); % decision criteria after each trial
    
    for k = 1: 10
        % initialize random seed
        seed = sum(round(clock));
        rng(seed)
        data = cell(3, 1);
        % initialize random seed
        heading = nan(num_trial, 1);
        choice = nan(num_trial, 1);
        stimulus = nan(num_trial, 1);
        bias_updata = nan(num_trial, 3);
              
        for i = 1: num_trial
            % initialize block randomlization for stimulus direction
            if mod(i, 3 * length(direction_condition)) == 1
                sti_index = 1;
                stimuli = randperm(3 * length( direction_condition));
            end
            if stimuli(sti_index) > (2 * length( direction_condition))
                heading(i) = 1;
            elseif stimuli(sti_index) > (1 * length( direction_condition))
                heading(i) = 2;
            else 
                heading(i) = 3;
            end
            if method == 6
                stimulus(i) = stimuli_list(stimuli(sti_index)) + bias_p(heading(i)); % symetric stimulus for figure S3 EF
            else
                stimulus(i) = stimuli_list(stimuli(sti_index));
            end
            sti_index = sti_index + 1;
            % updating dicision cretera
            if i == 1
                bias_d = bias_dd * ones(1, 3);
            else
                bias_d = bias_updata(i - 1, :);
            end
            perception_ob = normrnd(stimulus(i) - bias_p(heading(i)), thre); % noised perception
            choice(i) = perception_ob > bias_d(heading(i)); % choice for this trial 1: right, 0: left

            % reward strategy
            if method == 1
                correct_flag = sign(sign(choice(i) - 0.5) * (stimulus(i)));
            elseif method == 2
                correct_flag = sign(sign(choice(i) - 0.5) * (stimulus(i)));
                if heading(i) == 1 && stimulus(i) >= bias_p(1) && stimulus(i) <=0
                    correct_flag = 0;
                elseif heading(i) == 3 && stimulus(i) <= bias_p(3) && stimulus(i) >=0 
                    correct_flag = 0;
                end
            elseif method == 3
                correct_flag = sign(sign(choice(i) - 0.5) * (stimulus(i) -  bias_p(heading(i))));
            elseif method == 4
                correct_flag = sign(sign(choice(i) - 0.5) * (stimulus(i)));
                if heading(i) == 1 && stimulus(i) >= bias_p(1) && stimulus(i) <=0
                    correct_flag = 1;
                elseif heading(i) == 3 && stimulus(i) <= bias_p(3) && stimulus(i) >=0 
                    correct_flag = 1;
                end    
            elseif method == 5
                correct_flag = sign(sign(choice(i) - 0.5) * (stimulus(i)));
                if heading(i) == 1 && stimulus(i) >= bias_p(1) && stimulus(i) <=0
                    correct_flag = -1;
                elseif heading(i) == 3 && stimulus(i) <= bias_p(3) && stimulus(i) >=0 
                    correct_flag = -1;
                end
            elseif method == 6
                correct_flag = sign(sign(choice(i) - 0.5) * (stimulus(i) -  bias_p(heading(i))));
            end
            
            % reassign the index for reward: correct_flag: 1: reward, -1: no reward
            if correct_flag == 0
                correct_flag = sign(rand - 0.5);
            end
            if correct_flag == -1
                correct_flag = 0;
            end
            choice_prob = normcdf(perception_ob, bias_d(heading(i)), thre); % choice confidence
            
            % update decision criteria
            if i == 1
                bias_updata(i, :) = bias_d;
            else
                for hea = 1: 3
                    if hea == heading(i)
                        delta_bias = alpha * (correct_flag -  max(choice_prob, 1-choice_prob)) * sign(choice(i) - 0.5);
                        bias_updata(i, hea) = bias_updata(i-1, hea) - delta_bias;
                    else
                        bias_updata(i, hea) = bias_updata(i-1, hea);
                    end
                end
            end

        data_up = [stimulus(i), choice(i), 1];
        data{heading(i)}(end+1, :) = data_up; % data saved for psignifit fitting

        end
        
        bias_savve(:, :, k) = bias_updata;

        % fiting psychometric cureve using psignifit
        options.sigmoidName = 'norm';  
        options.expType     = 'YesNo'; 
        options.estimateType = 'mean';
        for session = 1: 2%50
            for con_hea= 1: 3
                data_fit = data{con_hea}(330 * (session-1) + 1: 330 * session, :);
                result = psignifit(data_fit, options);
                bias_psig(session, con_hea, k) = result.Fit(1);
            end
        end
    end
    bias_method_learn{method} = bias_savve;
    bias_method_psig{method} = bias_psig;
end


