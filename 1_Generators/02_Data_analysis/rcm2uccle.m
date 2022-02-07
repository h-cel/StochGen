function dataUcclestyle = rcm2uccle(data,ymd, agg)
%% rcm2uccle Transformation of data to the style of the Uccle time series.
%
%   This function transforms precipitation data into the Uccle style format. This is
%   necessary for some data analysis tools in BLBOX, in which precipitation
%   is always simulated in the Uccle style format.
%
%   Inputs:
%       data: precipitation data to be transformed; supposed to be a matrix
%       with 3 or 4 columns of time data, and one column of precipitation
%       data
%       ymd: time matrix, necessary in case the data does not contain time
%       rows
%       agg: original aggregation of the data. Supposed to be either hourly
%       ('1') or daily ('24')
%   
%   Outputs:
%       dataUcclestyle: data in the Uccle style format
%
%   Last updated by J. Van de Velde on 21/10/'20

%% Calculation

switch agg
    case '1'
        r = size(data(:,end),1);
        dataUcclestyle = zeros(r,13);
        dataUcclestyle(:,1) = ones(r,1)*6;
        dataUcclestyle(:,2) = ones(r,1)*6100;
        dataUcclestyle(:,7) = ones(r,1);
        dataUcclestyle(:,[3:6,8]) = data;
        dataUcclestyle((dataUcclestyle(:,8) < 0.1),8) = 0;
    case '24'
        r = size(data(:,end),1)*24;
        dataUcclestyle = zeros(r,13);
        iDay = 0:24:r;
        dataUcclestyle(:,1) = ones(r,1)*6;
        dataUcclestyle(:,2) = ones(r,1)*6100;
        dataUcclestyle(:,7) = ones(r,1);
        for j = 1:size(data,1)
            dataUcclestyle(iDay(j)+1:iDay(j+1),3:5) = repmat(ymd(j,:),24,1);
            dataUcclestyle(iDay(j)+1:iDay(j+1),6) = 0:23;
            dataUcclestyle(iDay(j)+1,8) = data(j,end);
        end
        dataUcclestyle((dataUcclestyle(:,8) < 0.1),8) = 0;

end

end

