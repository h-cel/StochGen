function X = addNoise(dataset,datatype,noise)
%% Add noise to precipitation, temperature, evapotranspiration
% This file adds noise. This is to make sure that no two moments are
% exactly the same, which is important for the copula fitting process.
%
% Noise calculation depends on the variable:
% Positive and negative noise for P&E BUT only positive noise to zero values
% Positive and negative noise for T
%
% INPUTS
% Dataset has to be in the following format:
% year month day var1 (...) varn
%
% Datatype specifies the input variable:
% "P" for Precipitation
% "E" for Evapotranspiration
% "T" for Temperature
%
% Noise specifies the order of magnitude of the noise
% e.g. noise = 0.001 will give noise values between 0 and 0.001
%
% Note that this file has been adapted: originally, it would also modify
% dry fraction, but this has been removed by J. Van de Velde. The original
% file is still available.

%% Initialization

X = nan(size(dataset));
X(:,1:3) = dataset(:,1:3);

%% Loop

for n = 1:length(datatype)
    data = dataset(:,3+n); %Selection of the column with the data to add noise to
    
    %Noise calculation based on variable
    
    if strcmp(datatype{n},'P') == 1 %Precipitation
        for i = 1:size(data,1)
            if data(i) <= noise
                data(i) = rand(1)*noise; %Adjust the data
            else
                data(i) = data(i) + rand(1)*2*noise -noise; %Substract off or add to the data
            end
        end
    elseif strcmp(datatype{n},'T') == 1 || strcmp(datatype{n},'E') == 1 %Temperature and Evaporation
        for i = 1:size(data,1)
            data(i)= data(i) + rand(1)*2*noise -noise; %Substract off or add to the data
        end
    else 
        error('Error. Datatype has to be P, E or T.')
    end
    X(:,3+n) = data; %Put the data back in de dataset
end
