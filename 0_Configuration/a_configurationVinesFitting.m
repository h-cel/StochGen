%% Configuration file to create workable data formats in orde to fit vines in R 

close all; clear all;
clc;

%% 1. Add C-Vine Copula code and its subfolders to Matlab's path

addpath(genpath('D:\Users\jpvdveld\Documents\PhD\Code\'), genpath('E:\Users\jpvdveld\Onderzoek\Data')); %Both Code and Data paths need to be added with their subfolders.

%% 2. Specify observation dataset, datatype and input variable

% Format of observation dataset:
% year month day var1 (...) varN
% Remark: Fill datagaps with NaN
data = matload('MPI-rcp45corr_xfs.mat');

% Datatype makes a distinction between:
% the historical period, 'hist'
% the future period, 'fut'
datatype = 'hist'; 

% Vars specifies the input variable:
% "P" for Precipitation
% "E" for Evapotranspiration
% "T" for Temperature
vars = {'E', 'T', 'P'};

%% 2. Define basal output filename

name = 'MPI-rcp45orig';

%% 3. Start the procedure

tic
if strcmp(datatype,'hist')
    prepareVinesdata(data,vars,name)
else
    for i = 1:1
        idata = data{i};
        iname = [name sprintf('_%d',i)];
        prepareVinesdata(idata,vars,iname)
    end
end
toc


