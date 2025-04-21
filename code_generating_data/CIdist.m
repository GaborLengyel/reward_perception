%% CIdist is also adapted from getConfRegions.m from psignifit library https://psignifit.sourceforge.net/WELCOME.html https://github.com/wichmann-lab/psignifit
function conf_Intervals = CIdist(x, Mass, para)
if length(x)==1 % catch when a dimension was fixed
    start=x;
    stop=x;
else               % probability Mass for each point
    cumMass = cumsum(Mass);                       % cumulated p-Mass

    % in the intervall are all points between the alpha and the
    % 1-alpha percentile.-> where alpha<cumMass<(1-alpha)
    confRegionM = cumMass > ((1-para)/2) & cumMass < (1-(1-para)/2);

    if any(confRegionM)
        % Now we have the confidence regions
        % put the borders between the nearest contained and the first
        % not contained point
        alpha       = (1-para)/2;
        startIndex  = find(confRegionM,1,'first');
        if startIndex > 1
%             start = (x(startIndex-1)+x(startIndex))/2 +(alpha-cumMass(startIndex-1))/margin(startIndex);
            start = (x(startIndex-1)+x(startIndex))/2;
        else
%             start = x(startIndex)+ alpha/margin(startIndex);
            start = x(startIndex);
        end

        stopIndex=find(confRegionM,1,'last');
        if stopIndex < length(x)
%             stop  = (x(stopIndex)+x(stopIndex+1))/2 + (1-alpha-cumMass(stopIndex))/margin(stopIndex+1);
            stop  = (x(stopIndex)+x(stopIndex+1))/2;
        else
%             stop = x(stopIndex)-alpha/margin(stopIndex); % this actually cannot happen because cumMass(end) ==1 > (1-(1-result.options.confP)/2) -> End is never in confRegionM
            stop = x(stopIndex);
        end
    else
        idx = find(cumMass >((1-para)/2), 1,'first');
        MMid = cumMass(idx)-Mass(idx)/2; % middle of the bar which contains both start and stop
        start = x(idx);
        stop = x(idx);
    end
    
end

conf_Intervals =[start,stop];
end