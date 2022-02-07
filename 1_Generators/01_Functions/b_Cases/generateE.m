function E =  generateE(T,P,R,thetasE,family,time,Eini,dataMonthE, copMatrix)
%% generateE This function generates E based on T and P.
%   This is based on the vine copula TPEpE, which is in this function used
%   with a combination of input data and randomly generated variables to
%   generate E.
%
%   Inputs:
%       T: temperature input series
%       P: Precipation input series
%       R: Ranks of data used to fit the vine
%       thetasE: parameters of the vine copula
%       family: copula families of the vine copula
%       time: time input series
%       Eini: initial Evaporation series
%       datamonthE: monthly data used to fit vine
%   Output:
%       E: Evaporation time series
%
%   Originally implemented by M. T. Pham and Hilde Vernieuwe
%   Last updated by J. Van de Velde on 12/07/'21: Updated sampling and
%   calculation of UEp

%% Initialization of time variables
ndays = length(time(:,1));
DOM = [31 28 31 30 31 30 31 31 30 31 30 31;
    31 29 31 30 31 30 31 31 30 31 30 31];
CDOM = cumsum(DOM,2); %Cumsum of each row
DOY = 2;
month = 1;
year = min(time(:,1));

%% Create UP & UT
% CDF's of P and T, so these can be uniformly sampled
UP = nan(ndays,1);
UT = nan(ndays,1);
for m = 1:12 %Monthly loop
    %Making cdf's
    [CDFT, XT] = ecdf(dataMonthE{m}(:,1));
    [CDFP, XP] = ecdf(dataMonthE{m}(:,2));
    pos = find(time(:,2) == m);
    for i = 1:length(pos)
        if T(pos(i)) < min(XT)
            UT(pos(i)) = 0;
        else
            UT(pos(i)) = max(CDFT(XT <= T(pos(i)) ));
        end
        if P(pos(i)) < min(XP)
            UP(pos(i)) = 0;
        else
            UP(pos(i)) = max(CDFP(XP <= P(pos(i)) ));
        end
    end
end

%% Create UE, UEp & E
%Initialization
E = nan(ndays,1);
E(1) = Eini;
UE = nan(ndays,1);
UEp = nan(ndays,1);
%Loop over all days
for i = 1:ndays-1
    
    if(leapy(year)==1)
        k = 2;
    else
        k = 1;
    end
    
    if DOY > CDOM(k,month) %checks if the month is over
        month = month + 1;
    end
    
    %Building of UEp
    [CDFEp, XEp] = ecdf(dataMonthE{month}(:,3));
    if E(i) < min(XEp)
        UEp(i) = 0;
    else
        UEp(i) = max(CDFEp(XEp <= E(i)));
    end
   
    %Building of UE
    UE(i) = sampleCvine([UT(i+1) UP(i+1) UEp(i)],rand(1),family(month,:,:),thetasE(month,:,:), copMatrix(month,:,:));  
    if isnan(UE(i)) % What if no data can be generated?
        if UT(i+1) == 1 || UT(i+1) == 0
            UT(i+1) = ksdensity(dataMonthE{m}(:,1), T(i+1), 'function','cdf');
        end
        if UP(i+1) == 1 || UP(i+1) == 0
            UP(i+1) = ksdensity(dataMonthE{m}(:,2), P(i+1), 'function','cdf', 'support', 'positive');
        end
        if UEp(i) == 1 || UEp(i) == 0
            UEp(i) = ksdensity(dataMonthE{m}(:,3), E(i), 'function','cdf');
        end
        UE(i) = sampleCvine([UT(i+1) UP(i+1) UEp(i)],rand(1),family(month,:,:),thetasE(month,:,:), copMatrix(month,:,:)); %New sample
    end
    E(i+1) = interp1(R{month}(:,4), dataMonthE{month}(:,4), UE(i), 'linear', 'extrap'); %Calculate E based on UT
           
    DOY = DOY+1;
    if DOY > CDOM(k,end) %Checks if the year is over
        DOY = 1;
        year = year+1;
        month = 1;
    end
       
end

end



