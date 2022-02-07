function [data] = RemovePOutliers(data, names)
%REMOVEPOUTLIERS Removes outliers in the future simulations
%   In the future data, some september values become very high. While this
%   might be realistic, the BL-model cannot handle this and thus the
%   outliers are removed.

%% Remove outliers

if strcmp(names{2}, '_1') == 1 
    uniqyrs = unique(data(:,1));
    meanP = zeros(length(uniqyrs),1);
    
    for y= 1:length(uniqyrs)
        datasep = data(data(:,1) == uniqyrs(y) & data(:,2) == 9,6);
        meanP(y) = mean(datasep);
    end
    
    data(data(:,1) == uniqyrs(meanP == max(meanP)) & data(:,2) == 9,6) = mean(meanP);
    
end

end

