function Y=percentiles(X,n)

% X: array of 10-minute rainfall measurements forming a storm
% n: the number of equal intervals in which to calculate the total amount
%    of rainfal. E.g. n=4: storm duration is divided in 4 equal parts. 

normt=[0 ((1:1:length(X))./length(X))];
Xcum=[0 cumsum(X)];
xi=0:1/n:1;
yi=interp1(normt,Xcum,xi,'linear');
Y=yi(2:length(yi))-yi(1:length(yi)-1);

end







