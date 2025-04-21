function [bias_trail, porst, result] = calcu_bias(data, prior)

% initial options in psignifit
    options             = struct;
    options.moveBorders    = 1;
%     options.maxBorderValue = exp(-20);
    if nargin == 2
        options.priors{1} = @(x) prior_up(prior.x, prior.y, x);
    end
    options.sigmoidName = 'norm'; 
    options.expType     = 'YesNo'; 
     % pass prior 
    options.estimateType   = 'mean'; 
    result = psignifit(data, options);
    % if fit using psignifit woudl be too slow compared to our method, 
    % you can try different settings in psignifit 
    % if you are curious about the relative performance between our method to psignifit
    % also you can try psignifitFast algorithm from their library:
    % result = psignifitFast(data, options); %fit
    bias_trail = result.Fit(1);
% get posterior y
    dims = [1,2];
    P = result.Posterior .* result.weight; 
    margdims = setdiff(1:5, dims);
    for j = length(margdims):-1:1
    P = sum(P, margdims(j));
    end
    P = squeeze(P);
    bias_dis = sum(P,2);
% get postierior bias
    bias_cor = result.X1D{1};
    porst.x = bias_cor;
    porst.y = bias_dis;
% % get postierior thre
%     thre_cor = result.X1D{2};
%     porst.x{2} = thre_cor;
%     porst.y{2} = thre_dis;
end