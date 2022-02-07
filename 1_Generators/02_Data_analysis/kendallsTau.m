function tau = kendallsTau(time,data1,data2)
%% kendallsTau Calculation of Kendall's tau
%   This function calculates Kendall's tau between two variables
%
%   Inputs:
%       time: matrix with dates for every timestep
%       data1: timeseries of variable 1 (vector)
%       data2: timeseries of variable 2 (vector)
%
%   Output:
%       tau: Kendall's tau value
%
%   Last updated by J. Van de Velde on 03/05/'21: documentation

%% Calculation

tau = nan(1,12);
for m = 1:12
    pos = find(time(:,2)==m);
    R = nan(length(pos),2);
    [CDF1, X1] = ecdf(data1(pos));
    [CDF2, X2] = ecdf(data2(pos));
    for i = 1:length(pos)
        R(i,1) = max(CDF1(X1 <= data1(pos(i))));
        R(i,2) = max(CDF2(X2 <= data2(pos(i))));
    end
    tmp = corr(R,'type','Kendall');
    tau(m) = tmp(1,2);
end
    
end

