%% Configuration file to launch simulations of P,T,E or Q via C-Vine Copula's
% Adapted by Jorn Van de Velde, 2018
%
% Last update by J. Van de Velde on 04/05/'21: save Psimtmp
%% Set-up

close all; clear all;
clc;

%% 1. Add C-Vine Copula code and its subfolders to Matlab's path

addpath(genpath('D:\Users\jpvdveld\Documents\PhD\Code\'), genpath('E:\Users\jpvdveld\Onderzoek\Data')); %Both Code and Data paths need to be added with their subfolders

%% 2. Specify observation dataset, input variable and period
% Format of the detrended observation dataset:
% year month day var1 (...) varN
% Remark: Fill datagaps with NaN
% This dataset is used as an input for all stochastic steps.
data = matload('MPI-rcp45corr_1_detrended.mat');

% Vars specifies the input variable:
% "P" for Precipitation
% "E" for Evapotranspiration
% "T" for Temperature
vars = {'E', 'T', 'P'};

%% 4. Specify PDM parameters

PDM_par = matload('paramPDM.mat');

%% 5. How many times repeating

repeat = 20; %Because of the stochastic nature of this step, repeats are needed to reduce uncertainty

%% 6. In case of scenario 3, how many years to determine?

yrs = 100;

%% 6. Define output name of simulation vectors
% Define in following order: 
% Base name - Repetition number - Evaporation - Temperature - Discharge - Precipitation
% Leave empty ('') when not applicable 

names = {'MPI-rcp45corr', '_1', 'Esim100test', 'Tsim100test', 'Qsim100test', 'Psim100aggtest','Psim100dailytest', 'time100', 'calstr2test'};

%% 7. Start the procedure

t10 = launchCases(data, vars, PDM_par, repeat, yrs, names, 0, '24');
save(['E:\Users\jpvdveld\Onderzoek\Data\3_cases\', 'MPI-rcp45.mat'],'t10') %t10? %Saves the time needed to run this.
