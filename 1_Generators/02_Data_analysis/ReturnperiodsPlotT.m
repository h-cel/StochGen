function [] = ReturnperiodsPlotT(datafile,origfile, agg, repetitions, fignumber)
%% returnperiodsPlot Calculates and plots return periods
%   This function sets up the calculation and plotting of the return
%   periods of temperature data.
%
%   Inputs:
%       datafile: file with the data (string)
%       origfile: file used for calibration (string)
%       agg: time aggregation level of the data (numerical)
%       repetitions: number of repetitions used for simulation

%% Set-up

addpath(genpath('E:\Users\jpvdveld\Onderzoek\Data'))

%% Loading & preprocessing

% Loading

data = matload(datafile);
orig = matload(origfile);
yrsdata = unique(data{1}(:,1));
if iscell(orig) == 1
    yrsorig = unique(orig{1}(:,1));
    orig = orig{1};
else
    yrsorig = unique(orig(:,1));
end

%% Calculation & Plot

figure_handle = figure(fignumber);

for i=1:repetitions
    for j = 1:repetitions
        [T_data, Ex_data] = returnperiods(data{i}(:,j+size(data{1},2)-repetitions), yrsdata, agg);
        p_data = semilogx(T_data, Ex_data, '*', 'Color', [209, 208, 208]/255);
        hold on
    end
end

if iscell(orig) == 1
    [T_orig, Ex_orig] = returnperiods(orig{1}(:,5), yrsorig, agg);
else
    [T_orig, Ex_orig] = returnperiods(orig(:,5), yrsorig, agg);
end


p_orig = semilogx(T_orig, Ex_orig, 'r*');


legend([p_data, p_orig], {'Simulated values', 'Original values'}, 'Location', 'northwest')
xlabel('Return period [year]')
ylabel('Temperature [°C]')
xlim([0 500])
set(gca, 'XScale', 'log');
title([num2str(agg) ' h return period values'])

set(gca,'fontsize',14)
fullfig(figure_handle)

end

