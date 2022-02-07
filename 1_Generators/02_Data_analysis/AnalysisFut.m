%% AnalysisFut
%   This script runs an analysis of the future data (MPI-rcp45 precipitation)
%   to understand the year-to-year variability of every month.
%
%   Last updated by J. Van de Velde on 03/06/'21: documentation

%% Data loading

origtmp = matload('MPI-rcp45corr.mat');
origdata=origtmp{1};

%% Initialisation

meanmat = nan(30,12);
varmat = nan(30,12);
autocovmat = nan(30,12);
zdpmat = nan(30,12);


%% Calculation

uniqyrs = unique(origdata(:,1));

for y = 1:length(uniqyrs)
    for k = 1:12
        origm = origdata(origdata(:,1) == uniqyrs(y) & origdata(:,2) == k,6);
        meanmat(y,k) = mean(origm);
        varmat(y,k) = var(origm);
        c1=cov([origm(1:length(origm)-1) origm(2:length(origm))]);
        autocovmat(y,k) = c1(1,2);
        zdpmat(y,k)=length(find(origm==0))./length(origm);
    end
end

%% Plots

labelstring = {'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'};

f1= figure(1);
boxplot(meanmat, 'labels', labelstring)
title('Future mean P [mm]')

f2= figure(2);
boxplot(varmat, 'labels', labelstring)
title('Future P variance [mm^2]')

f3= figure(3);
boxplot(autocovmat, 'labels', labelstring)
title('Future P autocovariance [mm^2]')

f4= figure(4);
boxplot(zdpmat, 'labels', labelstring)
title('Future ZDP (%)')
