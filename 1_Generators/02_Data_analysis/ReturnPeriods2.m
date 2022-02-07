function [] = ReturnPeriods2(datafile,origfile, aggdata, rep, aggplot,type,period, model, names)
%Returnperiods2 Calculates and plots return periods
%   This function calculates and plots the return periods of the
%   original and B-L-model-simulated data, using the 'showme' function of the BLBOX code set.
%
%   Input:
%       datafile: simulated data
%       origfile: original dataset, either observations or the set used for
%       calibration
%       aggdata: the aggregation level of the simulated data (string)
%       rep: the number of repetitions used in the simulated data
%       aggplot: vector with the aggregate levels to be used in the plots
%       type: 'seasonal' or 'yearly' plots
%       period: if the plot is seasonal, a vector with the seasons (1-4)
%       that have to be plotted
%       model: string with the model used
%
%   Last updated by J. Van de Velde on 05/005/'21: small updates to
%   documentation + bug fix for cell origfile

%% Set-up

addpath(genpath('E:\Users\jpvdveld\Onderzoek\Data'))

%% Loading & preprocessing

% Loading

data = matload(datafile);
orig = matload(origfile);

if iscell(orig)
    orig = orig{1};
end

% Removal of outliers in orig data

orig = RemovePOutliers(orig, names);

% Simulated data

if strcmp(aggdata,'24') == 1
    dataUcclestyle = zeros(length(data)*24*rep, 13);
else
    dataUcclestyle = zeros(length(data)*rep, 13);
end

for k = 1:rep
    switch aggdata %Transforms according to given aggregate level + aggregates so showme can handle all data
        case '24'
            p = data(:,[1:3 k+3]);
            time = p(:,1:3);
            dataUcclestyle(length(p)*24*(k-1)+1:length(p)*24*(k),:) = rcm2uccle(p,time,aggdata);        
        case '1'
            p = data(:,[1:4 k+4]);
            time = p(:,1:4);
            dataUcclestyle(length(p)*(k-1)+1:length(p)*(k),:) = rcm2uccle(p,time,aggdata);
    end
end

% Original data

switch aggdata
    case '1'
        obs=rcm2uccle(orig, orig(:,1:4), aggdata);
    case '24'
        obs=rcm2uccle(orig, orig(:,1:3), aggdata);
end



%% Calculation

legendlabel = {model};

h=showexpot(dataUcclestyle,aggplot,obs,type,period,legendlabel);

end

