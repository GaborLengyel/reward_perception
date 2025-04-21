%% Monkey data anaysis for Figure7, S2

clear
clc
dt = dir('monkey_data/*.mat');

for k = 1: length(dt) - 3
    [angle_all, ecc_all, fit] = fit_stan_offline(dt(1: k+3));
    paras = fit.extract();
    tau_m(k, 1) = quantile(paras.tau1, 0.5);
    tau_m(k, 2) = quantile(paras.tau_d, 0.5);
    tau_m(k, 3) = quantile(paras.tau2, 0.5);
    tau_ul(k, 1, 1) = quantile(paras.tau1, 0.16);
    tau_ul(k, 2, 1) = quantile(paras.tau_d, 0.16);
    tau_ul(k, 3, 1) = quantile(paras.tau2, 0.16);
    tau_ul(k, 1, 2) = quantile(paras.tau1, 0.84);
    tau_ul(k, 2, 2) = quantile(paras.tau_d, 0.84);
    tau_ul(k, 3, 2) = quantile(paras.tau2, 0.84);
end
% perceptual & decision bias
P1_posterior = paras.P1;
P2_posterior = paras.P2;
D_posterior = paras.D;
% regression weight
betae1_posterior = paras.E_1;
betae2_posterior = paras.E_2;
betan1_posterior = paras.N_1;
betan2_posterior = paras.N_2;

P1_plot_mean = quantile(P1_posterior, 0.5);
P1_plot_up = quantile(P1_posterior, 0.84);
P1_plot_low = quantile(P1_posterior, 0.16);

P2_plot_mean = quantile(P2_posterior, 0.5);
P2_plot_up = quantile(P2_posterior, 0.84);
P2_plot_low = quantile(P2_posterior, 0.16);

D_plot_mean = quantile(D_posterior, 0.5);
D_plot_up = quantile(D_posterior, 0.84);
D_plot_low = quantile(D_posterior, 0.16);

betae1_plot_mean = quantile(betae1_posterior, 0.5);
betae1_plot_up = quantile(betae1_posterior, 0.84);
betae1_plot_low = quantile(betae1_posterior, 0.16);

betae2_plot_mean = quantile(betae2_posterior, 0.5);
betae2_plot_up = quantile(betae2_posterior, 0.84);
betae2_plot_low = quantile(betae2_posterior, 0.16);

betan1_plot_mean = quantile(betan1_posterior, 0.5);
betan1_plot_up = quantile(betan1_posterior, 0.84);
betan1_plot_low = quantile(betan1_posterior, 0.16);

betan2_plot_mean = quantile(betan2_posterior, 0.5);
betan2_plot_up = quantile(betan2_posterior, 0.84);
betan2_plot_low = quantile(betan2_posterior, 0.16);

 


function [angle_all, ecc_all, fit] = fit_stan_offline(dt)

for k = 1: length(dt)
    name = fullfile(dt((k)).folder, dt((k)).name);
    load(name);
    heading_condition = unique(data.heading_direction);
    heading_condition = heading_condition(heading_condition ~= -888);
    heading_condition = heading_condition(~isnan(heading_condition)); 
    angle_all(k) = abs(heading_condition(1)); % find foe angle
    direction_condition = unique(data.stimulus_direction); % find stimulus direction
    x_loc = data.x_location;
    y_loc = data.y_location;
    ecc_all(k) = sqrt(x_loc.^2 + y_loc.^2); % find eccentricity
    amp_all(k) = 0.1;
    heading = [];
    direction = [];
    choice =[];
    for i = 1: length(data.misc_params)
        if data.misc_params(i)>=0 % valide trial only
            heading(i) = find(data.heading_direction(i) == heading_condition); % heading direction (left, neutral, right) for each trial
            direction(i) = find(data.stimulus_direction(i) == direction_condition); % stimulus direction for each trial
            choice(i) = data.choice(i); % monkey's choice for each trial
        end
    end
    kk = ones(1,3);
    stimu = [];
    choii = [];
    for i = 1: length(heading)
        for j = 1: 3
            if heading(i) == j
               stimu{j}(kk(j)) = direction_condition(direction(i));
               choii{j}(kk(j)) = choice(i);
               kk(j) = kk(j) + 1;
            end
        end
    end
%%
    options = struct; 
    options.poolMaxGap     = inf;            
    options.poolMaxLength  = inf;            
    options.poolxTol       = 0;
    % organize data
    for hra = 1:3
        data = [];
        data(:,1) = stimu{hra};
        data(:,2) = choii{hra};
        data(:,3) = ones(1,length(choii{hra}));
        data_stan(:,:,hra, k) = poolData(data, options);
    end
end
%%
    % hireracical model fitting
    monk_dat = struct('s', length(amp_all),...
                      'n1', squeeze(data_stan(:,3,1,:)),...
                      'direc1',squeeze(data_stan(:,1,1, :)),...
                      'choo1',squeeze(data_stan(:,2,1, :)),...
                      'n2', squeeze(data_stan(:,3,2,:)),...
                      'direc2',squeeze(data_stan(:,1,2, :)),...
                      'choo2',squeeze(data_stan(:,2,2, :)),...
                      'n3', squeeze(data_stan(:,3,3,:)),...
                      'direc3',squeeze(data_stan(:,1,3, :)),...
                      'choo3',squeeze(data_stan(:,2,3, :)),...
                      'angle', angle_all,...
                      'amp', amp_all,...
                      'ecc', ecc_all,...
                      'ok', ones(11,3));
    params = struct('file','monkey_fit_all.stan','data',monk_dat,'iter',50000,'chains',1);
    fit = stan(params, 'init', struct('P0', ones(length(dt), 1) * 10,  'phi', 15 * ones(3, length(dt)), 'tau1', 5, 'tau2', 5, 'tau_d', 5,...
        'P1', ones(length(dt),1) * -10, 'P2',ones(length(dt),1) * 10,'D', zeros(length(dt),1), 'lapse1', 0.1 * ones(length(dt),1), 'lapse2', 0.1 * ones(length(dt),1),...
        'theta1', 0.5 * ones(11, length(dt)), 'theta2', 0.5 * ones(11, length(dt)), 'theta3', 0.5 * ones(11, length(dt)), 'N_1', -1.5, 'N_2', 1.5, 'E_1', -0.5, 'E_2', 0.5));
    fit.verbose = true;
    addlistener(fit,'exit',@exitHandler);
    fit.block();

end