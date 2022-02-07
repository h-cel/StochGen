function [] = returnperiodsQ(names, casenumber, agg, repetitions, fignumber)
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

filenameP = [names{1}, '_', names{3}, '_', casenumber{1}, '.mat'];
filenameE = [names{1}, '_Esim100_', casenumber{1}, '.mat'];
origfile = [names{1}, '.mat'];
PDM_par = matload('paramPDM.mat');
inputs.A = 385;

dataP = matload(filenameP);
dataE = matload(filenameE);
orig = matload(origfile);
yrsdata = unique(dataP(:,1));
ndaysFut = length(dataP);
if iscell(orig) == 1
    yrsorig = unique(orig{1}(:,1));
    orig = orig{1};
else
    yrsorig = unique(orig(:,1));
end

ndaysOrig = length(orig);

% Removal P outliers

orig = RemovePOutliers(orig, names);

%% Calculation & Plot

% Simulations

figure_handle = figure(fignumber);
for i=1:repetitions
    for j = 1:repetitions
        P_rep = dataP(:,j+size(dataP,2)-repetitions);
        E_rep = dataE{i}(:,j+size(dataE{1},2)-repetitions);
        inputs.P = reshape(kron(P_rep/24, ones(1,24))', ndaysFut*24,1); %Hourly precipitation input data
        inputs.E = reshape(kron(E_rep/24, ones(1,24))', ndaysFut*24,1); %Hourly evaporation input data
        Q_rep = PDMPieter(inputs, PDM_par);
        [T_data, Ex_data] = returnperiods(Q_rep, yrsdata, agg);
        Q_data = semilogx(T_data, Ex_data, '*', 'Color', [209, 208, 208]/255);
        hold on
    end
end

% Original

P_orig = orig(:,end);
E_orig = orig(:,4);

inputs.P = reshape(kron(P_orig/24, ones(1,24))', ndaysOrig*24,1); %Hourly precipitation input data
inputs.E = reshape(kron(E_orig/24, ones(1,24))', ndaysOrig*24,1); %Hourly evaporation input data
Q = PDMPieter(inputs, PDM_par);
[T_orig, Ex_orig] = returnperiods(Q, yrsorig, agg);
Q_orig = semilogx(T_orig, Ex_orig, 'r*');


legend([Q_data, Q_orig], {'Simulated values', 'Original values'}, 'Location', 'northwest')
xlabel('Return period [year]')
ylabel('Discharge [m^3/s]')
xlim([0 1500])
set(gca, 'XScale', 'log');
title([num2str(agg) ' h return period values'])

set(gca,'fontsize',14)
fullfig(figure_handle)

end

