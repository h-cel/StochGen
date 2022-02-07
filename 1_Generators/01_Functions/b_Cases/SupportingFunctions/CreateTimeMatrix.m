function [time] = CreateTimeMatrix(startdate,enddate)
%CREATETIMEMATRIX Creates a time matrix 

%   Creates a time matrix based on a given startdate and
%   enddate

datestart=datenum(startdate);
datestop=datenum(enddate);
dates=datestart:1:datestop;
datesVec=datevec(dates);
time = datesVec(:,1:3);

end

