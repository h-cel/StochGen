function t = launchCases(data, vars, PDM_par, repeat, yrs, names, calstrPresent, agg)
%% launchCases Creates the set up for all cases used
%
%   Launchcases launches the different cases used in the Stochastic
%   Generator research of Jorn Van de Velde. It is based on former code by
%   Kato De Roos (2018).
%
%   The function consists of several sections, each of which loads in or
%   sets up data. In the last section, the different cases are launched.
%   These cases are:
%   - Case 0: Input P & E is used for discharge simulation
%   - Case 1: Input P and T is used to stochastically generate E, input P
%   and generated E are used for discharge simulation
%   - case 2: Input P is used to stochastically generate T, input P and
%   generated T are used to stochastically generate E, input P and
%   generated E are used for discharge simulation
%   - case 3: Input P is used to stochastically generate P, generated P is
%   used to generate T, generated P and T are used to generate E, generated
%   E and P are used for discharge simulation.
%
%   Inputs:
%       data: detrended input time series
%       vars: order of the variables
%       PDM_par: PDM parameters (for discharge simulation). Directly used in the cases
%       repeat: base number of repetitions. Directly used in the cases
%       yrs: years generated in case 3
%       names: struct, names to save or load files
%       calstrPresent: Boolean, indicates whether or not a calibration
%       string for the Bartlett-Lewis model in case 3 has already been
%       calibrated
%       agg: string, indicates the aggregation level of the Bartlett-Lewis
%       model. Directly used in the cases.
%
%   Outputs:
%       t: time to run all the cases
%
%   Last update by J. Van de Velde on 12/05/'21: Rewrote sampling to
%   account for RVineCopula

%% Dataset
% Many variables are loaded in to use them in the cases instead of this main script 
startdate = data(1,1:3);
enddate = data(end,1:3);
enddate(1) = startdate(1) + yrs -1;
time = data(:,1:3); %Used in the cases
for i = 4:length(vars)+3
    if strcmp(vars{i-3},'E')
        E = data(:,i); %Used in the cases
    elseif strcmp(vars{i-3},'T')
        T = data(:,i); %Used in the cases
    else
        Pnoise = data(:,i); %Used in the cases
        tmp = matload(sprintf('%s.mat',names{1})); %No detrended data needed for P
        if strcmp(names{2},'_1')
            P = tmp{1}(:,i); %Used in the cases
        else
            P = tmp(:,i);
        end
    end
end

%% Calstr

if calstrPresent == 1
    tmp = load([names{1} '_' names{9}]);
    calstr = tmp.CalStr; %Nog aanpassen? Anders kan niet alles ingeladen worden
else
    calstr = [];
end

%% Ranks
RE = cell(1,12); %Rank Evaporation
for m = 1:12
    filename = sprintf('%s%s_TPEpE_%d.mat',names{1},names{2},m);%Loads monthly data used to fit vine
    dataMonthE{m} = matload(filename);
    dataPos = [dataMonthE{m}(:,1) dataMonthE{m}(:,2) dataMonthE{m}(:,3) dataMonthE{m}(:,4)]; %Makes matrix out of the data of dataMonthE
  
    RE{m} = nan(size(dataPos,1), size(dataPos,2));
    % Makes ECDF's of the data
    [CDF1, X1] = ecdf(dataPos(:,1));
    [CDF2, X2] = ecdf(dataPos(:,2));
    [CDF3, X3] = ecdf(dataPos(:,3));
    [CDF4, X4] = ecdf(dataPos(:,4));
    for i = 1:length(dataPos(:,1))
        %Setting ranks
        RE{m}(i,1) = max(CDF1(X1 <= dataPos(i,1)));
        RE{m}(i,2) = max(CDF2(X2 <= dataPos(i,2)));
        RE{m}(i,3) = max(CDF3(X3 <= dataPos(i,3)));
        RE{m}(i,4) = max(CDF4(X4 <= dataPos(i,4)));
    end
end

RT = cell(1,12); %Rank Temperature
for m = 1:12
    filename = sprintf('%s%s_TpPT_%d.mat',names{1},names{2},m); %Loads monthly data used to fit vine
    dataMonthT{m} = matload(filename);
    dataPos = [dataMonthT{m}(:,1) dataMonthT{m}(:,2) dataMonthT{m}(:,3)]; %Makes matrix out of the data of dataMonthT
    margT{m} = dataMonthT{m}(:,3);
    
    RT{m} = nan(size(dataPos,1), size(dataPos,2));
    %Makes ECDF's out of the data
    [CDF1, X1] = ecdf(dataPos(:,1));
    [CDF2, X2] = ecdf(dataPos(:,2));
    [CDF3, X3] = ecdf(dataPos(:,3));
    for i = 1:length(dataPos(:,1))
        %Setting ranks
        RT{m}(i,1) = max(CDF1(X1 <= dataPos(i,1)));
        RT{m}(i,2) = max(CDF2(X2 <= dataPos(i,2)));
        RT{m}(i,3) = max(CDF3(X3 <= dataPos(i,3)));
    end
end

%% Vines
filename = sprintf('%s%s_TPEpE_vines.mat',names{1},names{2}); %Loads vines
tmp = load(filename);
tmp2 = reshape(tmp.fams, [12,4,4]);
thetasE = reshape(tmp.pars, [12,4,4]); %Used in the cases
copMatrixE = reshape(tmp.copmatrix, [12,4,4]); %Used in the cases

familyTPEpE = cell(size(tmp2));
for i = 1:size(tmp2,1) %Months
    for j = 2:size(tmp2,2) %Variables used
        for k = 1:(size(tmp2,3)-1)
            if tmp2(i,j,k) == 1
                familyTPEpE{i,k,j} = 'Gaussian';
            elseif tmp2(i,j,k) == 3
                familyTPEpE{i,k,j} = 'Clayton';
            elseif tmp2(i,j,k) == 4
                familyTPEpE{i,k,j} = 'Gumbel';
            elseif tmp2(i,j,k) == 5
                familyTPEpE{i,k,j} = 'Frank';
            elseif tmp2(i,j,k) == 13
                familyTPEpE{i,k,j} = 'Clayton 180°';
            elseif tmp2(i,j,k) == 14
                familyTPEpE{i,k,j} = 'Gumbel 180°';
            elseif tmp2(i,j,k) == 23
                familyTPEpE{i,k,j} = 'Clayton 90°';
            elseif tmp2(i,j,k) == 24
                familyTPEpE{i,k,j} = 'Gumbel 90°';
            elseif tmp2(i,j,k) == 33
                familyTPEpE{i,k,j} = 'Clayton 270°';
            elseif tmp2(i,j,k) == 34
                familyTPEpE{i,k,j} = 'Gumbel 270°';
            end
        end
    end
end

filename = sprintf('%s%s_TpPT_vines.mat',names{1},names{2});
tmp = load(filename);
tmp2 = reshape(tmp.fams, [12,3,3]);
thetasT = reshape(tmp.pars, [12,3,3]); %Used in the cases
copMatrixT = reshape(tmp.copmatrix, [12,3,3]); %Used in the cases

familyTpPT = cell(size(tmp2));
for i = 1:size(tmp2,1) %Months
    for j = 2:size(tmp2,2) %Variables used
        for k=1:size(tmp2,3)
            if tmp2(i,j,k) == 1
                familyTpPT{i,k,j} = 'Gaussian';
            elseif tmp2(i,j,k) == 3
                familyTpPT{i,k,j} = 'Clayton';
            elseif tmp2(i,j,k) == 4
                familyTpPT{i,k,j} = 'Gumbel';
            elseif tmp2(i,j,k) == 5
                familyTpPT{i,k,j} = 'Frank';
            elseif tmp2(i,j,k) == 13
                familyTpPT{i,k,j} = 'Clayton 180°';
            elseif tmp2(i,j,k) == 14
                familyTpPT{i,k,j} = 'Gumbel 180°';
            elseif tmp2(i,j,k) == 23
                familyTpPT{i,k,j} = 'Clayton 90°';
            elseif tmp2(i,j,k) == 24
                familyTpPT{i,k,j} = 'Gumbel 90°';
            elseif tmp2(i,j,k) == 33
                familyTpPT{i,k,j} = 'Clayton 270°';
            elseif tmp2(i,j,k) == 34
                familyTpPT{i,k,j} = 'Gumbel 270°';
            end
        end
    end
end

filename = sprintf('statsData_%s%s_E.mat',names{1},names{2});
tmp = matload(filename); 
dataStand = tmp.dataStand;

%% Cases

% tic
% run('case0.m') %No generation
% t(1) = toc;
% 
% tic
% run('case1.m') %Evaporation generation 
% t(2) = toc;
% 
% tic
% run('case2.m') %Temperature and Evaporation generation
% t(3) = toc;

% tic
run('case3.m') %Precipitation, Temperature and Evaporation Generation, Evaporation is simulated in one time series
% t(4) = toc;

end

