function prepareVinesdata(data,vars,name)
%% prepareVinesdata This function prepares the data for the vine copula fit
%   Two vines, TPEpE and TpPT, are used. To make sure the fitting process
%   will be correct, the data has to be unique (noise is added), within-month
%   trends have to be removed and autocorrelation has to be checked.
%
%   Input:
%       data: dataset used for vines fitting
%       vars: order of the variables in the data
%       name: base filename used
%
%   Originally implemented by M. T. Pham and H. Vernieuwe (?)
%   Last updated by J. Van de Velde on 08/04/'21: documentation

%% Adding noise
dataNoise = addNoise(data,vars,0.001);
time = dataNoise(:,1:3);

%% Computing statistics
[dataDetrend, ~, pTrend, ~, hLbq] = statsData(dataNoise,vars,name); %Detrended data, result of anova test and result of autocorrelation are used

%% Checking existence of within-month trends. If so, detrended data is used.
detrend = nan(1,length(vars));
dataNew = nan(size(dataNoise));
dataNew(:,1:3) = time;
for n = 1:length(vars)
    if strcmp(vars{n},'T') == 1
        if sum((pTrend(n,:)<0.001)) > 0 % More than zero low p-values: there is a trend
            T = dataDetrend(:,n+3);
            dataNew(:,n+3) = T;
            disp('Existence of within-month trends for Temperature so detrended data will be used.');
            detrend(n) = 1; 
        else
            T = dataNoise(:,n+3); %Use of original, noise-added, data
            dataNew(:,n+3) = T;
            detrend(n) = 0;
        end
    elseif strcmp(vars{n}, 'E') == 1
        if sum((pTrend(n,:)<0.001)) > 0 % More than zero low p-values: there is a trend
            E = dataDetrend(:,n+3);
            dataNew(:,n+3) = E;
            disp('Existence of within-month trends for Evapotranspiration so detrended data will be used.');
            detrend(n) = 1;
        else
            E = dataNoise(:,n+3); %Use of original, noise-added, data 
            dataNew(:,n+3) = E;
            detrend(n) = 0;
        end
    elseif strcmp(vars{n}, 'P') == 1
        if sum((pTrend(n,:)<0.001)) > 0 % More than zero low p-values: there is a trend
            P = dataDetrend(:,n+3);
            dataNew(:,n+3) = P;
            disp('Existence of within-month trends for Precipitation so detrended data will be used.');
            detrend(n) = 1;
        else
            P = dataNoise(:,n+3); %Use of original, noise-added, data 
            dataNew(:,n+3) = P;
            detrend(n) = 0;
        end
    end
end

outputfile = [name '_detrended'];
save(['E:\Users\jpvdveld\Onderzoek\Data\2_detrended+vines\', outputfile],'dataNew');

%% Checking existence of autocorrelation
som = zeros(1,length(vars));
nyrs = length(unique(time(:,1)));
for n = 1:length(vars)
    if strcmp(vars{n},'T') == 1
        for i = 1:12
            som(n) = som(n) + sum(sum(hLbq{i,n}));
        end
        if som(n) > (nyrs*5*12)/8 %Why this?
            disp('Autocorrelation exists for Temperature.');
        end
    elseif strcmp(vars{n}, 'E') == 1
        for i = 1:12
            som(n) = som(n) + sum(sum(hLbq{i,n}));
        end
        if som(n) > (nyrs*5*12)/8
            disp('Autocorrelation exists for Evatranspiration.');
        end
    elseif strcmp(vars{n}, 'P') == 1
        for i = 1:12
            som(n) = som(n) + sum(sum(hLbq{i,n}));
        end
        if som(n) > (nyrs*5*12)/8
            disp('Autocorrelation exists for Precipitation.');
        end
    end
end

%% Prepare data structures for TpPT and TPEEp vines
prepareTPEpE([time T P E], name);
prepareTpPT([time T P ], name);

end

