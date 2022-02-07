function  [dataP, CalStr] = simBLBOX(model, numpar, data, vars, startdate, enddate, repeat, agg, calstr)
%% simBLBOX Simulates and postprocesses precipitation according to the Bartlett-Lewis process
%   simBLBOX calibrates and simulates a given Bartlett-Lewis model
%   
%   Inputs:
%       model: model used (string)
%       numpar: number of parameters used in the model (double)
%       data: timeseries used for calibration
%       vars: variables in the time series, used to define the column with
%       precipitation
%       startdate: start date of the simulation (vector)
%       enddate: end date of the simulation (vector)
%       repeat: number of repetitions used
%       agg: desired aggregation level of the output
%       calstr: (optional) information of an earlier calibration
%   Outputs:
%       dataP: simulated time series
%       CalStr: information on the calibration. If calstr as an input was
%       present, it will be used as this output.
%
%   Originally implemented by W.-J. Vanhaute
%   Last update by J. Van de Velde on 30/03/'21: documentation

%% Transforming data in Uccle style

r = size(data,1)*24;
dataUcclestyle = zeros(r,13);
ymd = data(:,1:3);
iDay = 0:24:r;
for k = 1:3
    if strcmp(vars{k}, 'P') == 1
        iP = k;
    end
end

dataUcclestyle(:,1) = ones(r,1)*6;
dataUcclestyle(:,2) = ones(r,1)*6100;
dataUcclestyle(:,7) = ones(r,1);
for j = 1:size(data,1)
    dataUcclestyle(iDay(j)+1:iDay(j+1),3:5) = repmat(ymd(j,:),24,1);
    dataUcclestyle(iDay(j)+1:iDay(j+1),6) = 0:23;
    dataUcclestyle(iDay(j)+1,8) = data(j,iP+3);
end
dataUcclestyle((dataUcclestyle(:,8) < 0.1),8) = 0; 

%% Calibration of BLBOX

if isempty(calstr)
    zdpflag=1;
    CalStr = momFit(model,dataUcclestyle,'OBJ_fn_4','SCE',1:12,20,getRules,getBoundaries(model),[],[],zdpflag);
else
    CalStr = calstr;
end

%% Defining parameters (x) corresponding with best objective function (z)

z = CalStr.z(:,:);
[A,~] = find(z==min(z));
param = nan(12,numpar);

for i = 1:12
    param(i,:) =  reshape(CalStr.x(A(i),:,i),numpar,1)';
end

%% Simulating precipitation

datestart=datenum(startdate);
datestop=datenum(enddate);
dates=datestart:1:datestop;
datesVec=datevec(dates);
time = datesVec(:,1:3);

switch agg
    % Daily level
    case '24'
        dataP = nan(length(time(:,1)),repeat+3);
        for j = 1:repeat 
            synthrain = BartlettLewis(param,startdate(1),enddate(1),model);                                                                                           
            row = size(synthrain,1);                                                                                                                                  
            step = 0:24:row;                                                                                                                                          
            for i = 1:(row/24)                                                                                                                                        
                dataP(i,j+3) = sum(sum(synthrain((step(i)+1):step(i+1),8:13)));
                if (dataP(i,j+3) < 0.001 && dataP(i,j+3) ~=0) % Low-value issues
                    dataP(i,j) = 0;
                end
            end                                                                                                                                                              
        end
        dataP(:,1:3) = time;
        % Hourly level
    case '1'
        dataP = nan(length(time(:,1))*24,repeat+4);
        for j = 1:repeat
            synthrain = BartlettLewis(param,startdate(1),enddate(1),model);
            row = size(synthrain,1);
            for i = 1:row
                dataP(i,j+4) = sum(synthrain(i,8:13));  %aggregates per hour
                dataP(i,1:3) = time(ceil(i/24),:);
                dataP(i,4) = i-floor(i/24)*24-1;
                %Clean up
                if dataP(i,4) == -1 %Unrealistic hours
                    dataP(i,4) = 23;
                end
                if (dataP(i,j+4) < 0.001 && dataP(i,j+4) ~=0) % Low-value issues
                    dataP(i,j) = 0;
                end
            end
        end
    % 10 min
    case '1/6'
        dataP = cell(1,20);
        for j = 1:repeat
            dataP{j} = nan(length(time(:,1))*24,10);
            synthrain = BartlettLewis(param,startdate(1),enddate(1),model);
            row = size(synthrain,1);
            for i = 1:row
                dataP{j}(i,:) = synthrain(i,[3:6,8:13]);
                for k = 1:6 % Low-value issues
                    if (dataP{j}(i,k+4) < 0.001 && dataP{j}(i,k+4) ~=0)
                        dataP{j}(i,k+4) = 0;
                    end
                end
            end
        end
end

%% Solving issue of low values
 
for i = 1:size(dataP,1)
    for j = 4:size(dataP,2)
        if (dataP(i,j) < 0.001 && dataP(i,j) ~=0)
            dataP(i,j) = 0;
        end
    end
end


end

