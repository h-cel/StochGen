%% MomentsP
%   This script loads the necessary data for an analysis of the moments of
%   stochastically generated P.
%
%   Last update by J. Van de Velde on 22/10/'20

%% Setup and data
 
clear all; close all; clc;
addpath(genpath('E:\Users\jpvdveld\Onderzoek\Data')) %Both Code and Data paths need to be added with their subfolders.
time = matload('Ucclecorr_xho_time_c1.mat');

%% Comparison
% Nog eens op te ruimen -> vooral CalMomentsP zelf veel gebruikt de laatste
% tijd

% Processed data
orig_tmp = matload('Ucclecorr_xho.mat'); %Original time series used for BL parameterisation
if isa(orig_tmp, 'cell') == 1
    orig = rcm2uccle(orig_tmp{1}(:,6),time, '24');
else
    orig = rcm2uccle(orig_tmp(:,6), time, '24');
end
% Original hourly data
%orig_tmp = matload('P120_hourly.mat'); %Original time series used for BL parameterisation
% orig = rcm2uccle(orig_tmp,time, '1');

StochGen = 'Ucclecorr_xho_Psim100_RBL2_c3.mat'; %Stochastically generated time series
stats1 = CalMomentsP(StochGen,20,[24 48 72],orig,1, '24');
