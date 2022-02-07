function T =  generateT(P,RT,thetasT,family,time,Tini,dataMonthT, copMatrix)
%% This function generates T based on Tp and P.
% This is based on the vine copula TpPT, which is in this function used
% with a combination of input data and randomly generated variables to
% generate T.
%
%   Inputs:
%       P: Precipation input series
%       R: Ranks of the data used to fit the vine
%       thetasT: parameters of the vine copula
%       family: copula families of the vine copula
%       time: time input series
%       Tini: initial Temperature series
%       datamonthT: monthly data used to fit vine
%   Output:
%       T: temperature time series
%
%   Originally implemented by M. T. Pham and H. Vernieuwe
%   Last update by J. Van de Velde on 23/04/'21: Reverted sampling to
%   CDVineCopula matrix definition

%% Initialization of time variables
% Time variables
ndays = length(time(:,1));
DOM = [31 28 31 30 31 30 31 31 30 31 30 31;
    31 29 31 30 31 30 31 31 30 31 30 31];
CDOM = cumsum(DOM,2); %cumsum of each row
DOY = 2;
month = 1;
year = min(time(:,1));

%% Create UP
% CDF of P, so this variable can be uniformly sampled.
UP = nan(ndays,1);
for m = 1:12 % Monthly loop
    %Making cdf's
    [CDFP, XP] = ecdf(dataMonthT{m}(:,2));
    pos = find(time(:,2) == m);
    for i = 1:length(pos)
        if P(pos(i)) < min(XP)
            UP(pos(i)) = 0;
        else
            UP(pos(i)) = max(CDFP(XP <= P(pos(i) )));
        end
    end
end

%% Create UT, UTp & T
%Initialization
T = nan(ndays,1);
T(1) = Tini;
UT = nan(ndays,1);
UTp = nan(ndays,1);

%Loop over all days
for i = 1:ndays-1
    
    if(leapy(year)==1)
        k = 2;
    else
        k = 1;
    end
    
    if DOY > CDOM(k,month) %Checks if the month is over
        month = month + 1;
    end
    
    %Building of UTp
    [CDFTp, XTp] = ecdf(dataMonthT{month}(:,1));
    if T(i) < min(XTp)
        UTp(i) = 0;
    else
        UTp(i) = max(CDFTp(XTp <= T(i)));
    end
    
    %Building of UT
    UT(i) = sampleCvine([UTp(i) UP(i+1)],rand(1),family(month,:,:),thetasT(month,:,:), copMatrix(month,:,:));
    if isnan(UT(i)) %What if no data can be generated?
        if UTp(i) == 1 || UTp(i) == 0
            UTp(i) = ksdensity(dataMonthT{m}(:,1), T(i), 'function','cdf');
        end
        if UP(i+1) == 1 || UP(i+1) == 0
            UP(i+1) = ksdensity(dataMonthT{m}(:,2), P(i+1), 'function','cdf');
        end
        UT(i) = sampleCvine([UTp(i) UP(i+1)],rand(1),family(month,:,:),thetasT(month,:,:), copMatrix(month,:,:)); % New sample
    end
    T(i+1) = interp1(RT{month}(:,3), dataMonthT{month}(:,3), UT(i), 'linear', 'extrap'); %Calculate T based on UT

    DOY = DOY+1;
    if DOY > CDOM(k,end) %Checks if the year is over
        DOY = 1;
        year = year+1;
        month = 1;
    end
    
end

end



