function [] = returnperiodsPlot(datafile,origfile, agg, repetitions, fignumber, names)
%% returnperiodsPlot Calculates and plots return periods
%   This function sets up the calculation and plotting of the return
%   periods of the precipitation data.
%
%   Inputs:
%       datafile: file with the data (string)
%       origfile: file used for calibration (string)
%       agg: time aggregation level of the data (numerical)
%       repetitions: number of repetitions used for simulation
%
%   Last update by J. Van de Velde on 05/05/'21: documentation added

%% Set-up

addpath(genpath('E:\Users\jpvdveld\Onderzoek\Data'))

%% Loading & preprocessing

% Loading

data = matload(datafile);
orig = matload(origfile);
yrsdata = unique(data(:,1));
if iscell(orig) == 1
    yrsorig = unique(orig{1}(:,1));
    orig = orig{1};
else
    yrsorig = unique(orig(:,1));
end

% Removal P outliers

orig = RemovePOutliers(orig, names);

%% Calculation & Plot

figure_handle = figure(fignumber);

for j = 1:repetitions
    [T_data, Ex_data] = returnperiods(data(:,j+size(data,2)-repetitions), yrsdata, agg);
    p_data = semilogx(T_data, Ex_data, '*', 'Color', [209, 208, 208]/255);
    hold on
end

if iscell(orig) == 1
    [T_orig, Ex_orig] = returnperiods(orig{1}(:,end), yrsorig, agg);
else
    [T_orig, Ex_orig] = returnperiods(orig(:,end), yrsorig, agg);
end


p_orig = semilogx(T_orig, Ex_orig, 'r*');


legend([p_data, p_orig], {'Simulated values', 'Original values'}, 'Location', 'northwest')
xlabel('Return period [year]')
ylabel('Precipitation [mm]')
xlim([0 1500])
set(gca, 'XScale', 'log');
title([num2str(agg) ' h return period values'])

set(gca,'fontsize',14)
fullfig(figure_handle)

end

