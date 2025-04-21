
function [data,heading_index] = simu_stan(bias, thre, cons, num_trial, seed)
%%
rng(seed)

direction_condition = -40: 8: 40; % stimulus direction
stimuli_list = [direction_condition, direction_condition, direction_condition];

if cons == 0 % changing decision bias
    bias_d0 = 10 * sin((1: num_trial)/2000 * pi);
elseif cons == 1 % constant decision bias
    bias_d0 = 10 * ones(1, num_trial);
end

options = struct; 
options.poolMaxGap     = inf;            
options.poolMaxLength  = inf;            
options.poolxTol       = 0;
data = cell(3, 1);

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
        heading_index(i) = heading;
        data{heading}(end + 1, :) = data_up;
        
end
end
