function [Tsort, EXsort] = returnperiods(data,yrs, agg)
%% returnperiods Calculates return periods
%
%   This function calculates empirical return periods (yearly block maxima) based on the data.
%
%   Inputs:
%       data: dataset
%       yrs: list of years to use in the dataset
%       agg: the aggregation level of the dataset
%
%   Outputs:
%       Tsort: sorted (ascending) return periods of the data
%       Exsort: sorted (ascending) extreme values corresponding to the
%       return periods in Tsort
%
%   Last update by J. Van de Velde on 22/10/'21

%% Calculation

numberextrema = floor(length(data)/365/(24/agg));
extrema = nan(numberextrema,1);
cnt = 0;

for i = 1:numberextrema
    if(leapy(yrs(i))==1)
       low = 365*(i-1)*(24/agg) + 1 + cnt; %Lower and upper boundary of a year
       up = 365*i*(24/agg) + 1 + cnt;
       cnt = cnt + 1; % To make sure that leap years don't distort the data
    else
       low = 365*(i-1)*(24/agg) + 1 + cnt;
       up = 365*i*(24/agg) + cnt; 
    end
    tmp = data(low:up);
    extrema(i) = max(tmp); %Selection of extremes per year
end

% Return periods
T = (length(extrema)+1)./(1:length(extrema));       % Denominator = rank

% Sort
[EXsort, ~] = sort(extrema);
Tsort = sort(T);


end

