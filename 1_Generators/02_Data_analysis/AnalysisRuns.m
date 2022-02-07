%% Analysis runs
%   This script runs all analyses of the stochastic Generator series
%   Settings are based on the standard setup (24 hours)
%
%   Last update by J. Van de Velde on 05/05/'21

%% Settings

names = {'MPI-rcp45corr', '', 'Psim100daily', 'time100'}; %Second argument: '_1' if BA future setup
casenumber = {'c3'};
filenameP = [names{1}, '_', names{3}, '_', casenumber{1}, '.mat'];
filenameT = [names{1}, '_Tsim100_', casenumber{1}, '.mat'];
filenameE = [names{1}, '_Esim100_', casenumber{1}, '.mat'];
origfile = [names{1}, '.mat'];
timename = [names{1}, '_', names{4}, '_', casenumber{1}, '.mat'];
casename = '100';
aggdata = '24';
rep = 20;

%% Analyses

% Moments

MomentsT(names{1}, casenumber, rep, timename, casename, 18);
MomentsE(names{1}, casenumber, rep, timename, casename, 19);
CalMomentsP(filenameP, rep, [24 48 72], origfile, 1, aggdata, 20); %Can be extended to include subdaily moments, origfile should then be specified

% Kendall's tau

kendallsTauPlot(names{1}, timename, casenumber, casename,rep, 21)

% PDFs

score = pdfPlot(names{1}, timename, casenumber, casename, 22);

% Return periods

returnperiodsPlot(filenameP, origfile, str2double(aggdata), rep, 23, names)
ReturnperiodsPlotT(filenameT, origfile, str2double(aggdata), rep, 23)
ReturnperiodsPlotE(filenameE, origfile, timename, str2double(aggdata), rep, 23)
returnperiodsPlotMonthly(filenameP, origfile, str2double(aggdata), rep, 23, names)
ReturnPeriods2(filenameP, origfile, aggdata, rep, 24, 'yearly', 1, 'RBL2', names)
ReturnPeriods2(filenameP, origfile, aggdata, rep, 24, 'seasonal', [1 2 3 4], 'RBL2', names)

% Q

returnperiodsQ(names, casenumber, 24, rep, 1)
MomentsQ(names, casenumber, rep, timename, casename, 18);


