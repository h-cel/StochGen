function moments=StatCalcm3(inputfile,agg,zdpflag,outputfile,options)
%% StatCalm3 This function calculates rainfall-statistics for given aggregation levels
%
%   This files calculates rainfall statistics, although it can also be used
%   for other variables.
%
%   Inputs:
%       inputfile: matfile - Uccle-style (0.1 mm)
%       agg: array of aggregation levels (minimum is 1/6 hour)
%       zdpflag: Boolean to indicate how zdp is calculated (0 = days are
%       dry when 0 mm, 1 = days are dry when < 0.1 mm)
%       outputfile: matfile to save results to
%       options='MO'(default) : pool data according to the months before
%           calculating properties, calculated properties are displayed for each
%           month
%           '-' : no subdivision of data set, properties are calculated for the whole
%           dataset. 
%
%   Outputs:
%       moments=3D matrix; for each aggregation level a 12x6 matrix -
%           rows = months; columns = mean-var-covar-autocorr-zdp-skewness
%       
%   Last updated by J. Van de Velde, 23/09/'20

%% Load input; Rainfall time series
if ischar(inputfile)==1;
    Obs=matload(inputfile);
else
    Obs=inputfile;
end

aggl=1/length(Obs(1,8:end));

%% Calculations

if nargin<4 || isempty(options) || strcmpi(options,'mo')
    % Initialize output
    moments=zeros(5,length(agg),12);

    % Seperate data according to month and calculate properties at specified
    % levels of aggregation
    for month=1:12

        % Extract data in certain
        pos=find(Obs(:,4)==month);
        % Number of data points
        n=length(Obs(pos,:));

        % For each level of aggregation, calculate ...
        for i=1:length(agg);

            % Reorganize matrix for better handling
            if n/agg(i)==round(n/agg(i))
                data=reshape(Obs(pos,8:end)',1/aggl*agg(i),n/agg(i))';
            else 
                data=reshape(Obs(pos(1:end-agg(i)*(n/agg(i)-floor(n/agg(i)))),8:end)',1/aggl*agg(i),floor(n/agg(i)))';
            end
            data=sum(data,2);        

            % Mean
            m=nanmean(data);
            % Variance
            s2=var(data, 'omitnan');
            % Covariance (lag-1)
            c1=cov([data(1:length(data)-1) data(2:length(data))], 'omitrows');
            c1=c1(1,2);
            % Covariance (lag-2)
%             c2=cov([data(1:length(data)-2) data(3:length(data))]);
%             c2=c2(1,2);
            % Covariance (lag-2)
    %         c3=cov([data(1:length(data)-3) data(4:length(data))]);
    %         c3=c3(1,2);
            % Correlation (lag-1-3)
    %         corr1=c1/s2;
    %         corr2=c2/s2;
    %         corr3=c3/s2;    
            % Zero depth probability
            % origineel Willem-Jan (effectief nulwaarden in reeks ukkel)
            if (zdpflag ==0)
            zdp=length(find(data==0))./length(data);
            % aangepast Hilde (voor climateruns)
            else
                zdp = length(find(data<0.1))./length(data);
            end
            % Skewness
    %         skew=skewness(data);
            % Third central moment
            m3=moment(data,3);      


            moments(:,i,month)=[m s2 c1 zdp m3];
        end

    end

else
    % Initialize output
    moments=zeros(5,length(agg));

    % Number of data points
    n=length(Obs);

    % For each level of aggregation, calculate ...
    for i=1:length(agg);

        % Reorganize matrix for better handling
        if n/agg(i)==round(n/agg(i))
            data=reshape(Obs(:,8:end)',1/aggl*agg(i),n/agg(i))';
        else 
            shift=agg(i)*(n/agg(i)-floor(n/agg(i)));
            data=reshape(Obs(1:end-uint64(shift),8:end)',1/aggl*agg(i),floor(n/agg(i)))';
        end
        data=sum(data,2);        

        % Mean
        m=nanmean(data);
        % Variance
        s2=var(data, 'omitnan');
        % Covariance (lag-1)
        c1=cov([data(1:length(data)-1) data(2:length(data))], 'omitrows');
        c1=c1(1,2);
        % Covariance (lag-2)
%         c2=cov([data(1:length(data)-2) data(3:length(data))]);
%         c2=c2(1,2);
        % Covariance (lag-2)
%         c3=cov([data(1:length(data)-3) data(4:length(data))]);
%         c3=c3(1,2);
        % Correlation (lag-1-3)
%         corr1=c1/s2;
%         corr2=c2/s2;
%         corr3=c3/s2;    
        % Zero depth probability
        zdp=length(find(data==0))./length(data);
        % Skewness
%         skew=skewness(data);
        % Third central moment
        m3=moment(data,3);      


        moments(:,i)=[m s2 c1 zdp m3];
    end
end
   
    
%% Save output to file
if  nargin>3
    if ~isempty(outputfile) 
    save(outputfile,'moments');
    end
end

