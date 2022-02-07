function [] = MonthlyMeans(datafile,origfile, aggdata, rep)
%MonthlyMeans Calculates and plots the monthly mean rainfall
%   This function calculates and plots the monthly mean rainfall of both
%   original and B-L-model-simulated data, using the 'showme' function of the BLBOX code set.
%
%   Input:
%       datafile: simulated data
%       origfile: original dataset, either observations or the set used for
%       calibration
%       aggdata: the aggregation level of the simulated data
%       rep: the number of repetitions used in the simulated data
%
%   Last updated by J. Van de Velde on 18/11/'20: documentation added

%% Set-up

addpath(genpath('E:\Users\jpvdveld\Onderzoek\Data'))

%% Loading & preprocessing

% Loading

data = matload(datafile);
orig = matload(origfile);

if length(orig)> length(data)
    startdate = data(1,1:(size(data,2)-rep));
    enddate = data(end,1:(size(data,2)-rep));
    orig_start = find(sum(orig(:,1:(size(data,2)-rep)) == startdate, 2) == 4);
    orig_end = find(sum(orig(:,1:(size(data,2)-rep)) == enddate, 2) == 4);
    orig = orig(orig_start:orig_end,:);
end

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
%         case '1/6' % Not completed yet
%             dataUcclestyletmp = zeros(r,13);
%             dataUcclestyletmp(:,1) = ones(r,1)*6;
%             dataUcclestyletmp(:,2) = ones(r,1)*6100;
%             dataUcclestyletmp(:,7) = ones(r,1);
%             dataUcclestyletmp(:,[3:6,8:end]) = data{k};
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

dataUcclestyle(:,8:13) = dataUcclestyle(:,8:13)/rep;
f1 = figure(1);
subplot(1,2,1)
title('Simulations')
showme(dataUcclestyle, 'monthlymean');
ylim([0 250])
subplot(1,2,2)
title('Observations')
showme(obs, 'monthlymean');
ylim([0 250])

end

