function [dataDetrend, dataStand, pTrend, pDetrend, hLbq] = statsData(data,datatype,modelname)
%% This function calculates the statistics of the dataset.

%Save location
save_loc = 'E:\Users\jpvdveld\Onderzoek\Data\2_detrended+vines\';

%% Standardisation of dataset

%Initialization
nvar = size(data,2)-3; %Without time
dataDetrend = nan(size(data,1),nvar+3);
dataDetrend(:,1:3) = data(:,1:3);
dataStand = nan(size(data,1),nvar*2+3); %nvar*2, because for each variable both the mean and standard deviation will be saved.
dataStand(:,1:3) = data(:,1:3);
dataDayTrend = cell(12,31); % Cells for each month and day
dataDayDetrend = cell(12,31);

%Loop to detrend
for i=1:size(data,1)
    pos = find(data(:,2)==data(i,2) & data(:,3)==data(i,3)); %Selection of month and day, for each year
    dataDetrend(i,4:end) = (data(i,4:end) - ... % Detrended for this day, over all years
        nanmean(data(pos,4:end)))./ nanstd(data(pos,4:end));
    dataStand(i,4:end) = [nanmean(data(pos,4:end)) nanstd(data(pos,4:end))]; %Saving of both mean and st.dev. for each variable
    dataDayTrend{data(i,2), data(i,3)} = [dataDayTrend{data(i,2), data(i,3)};...
        data(i,4:end)]; %Adding each specific day under the previous ones.
    dataDayDetrend{data(i,2), data(i,3)} = [dataDayDetrend{data(i,2), data(i,3)};...
        dataDetrend(i,4:end)]; %Adding each specific detrended day under the previous ones.
end
out.dataStand = dataStand;

%% Statistical tests on normality, homoscedasticity and autocorrelation

%Initialization
months = {'January', 'February', 'March', 'April', 'May', 'June', 'July',...
    'August', 'September', 'October', 'November', 'December'};
nyrs = length(unique(data(:,1)));
dataMonthTrend = cell(12,nvar);
dataMonthDetrend = cell(12,nvar);
hLillieTrend = cell(1,nvar);
hLillieDetrend = cell(1,nvar);
pTrend = nan(nvar,12); %p-value
pDetrend = nan(nvar,12);
pVarTrend = nan(nvar,12);
pVarDetrend = nan(nvar,12);
hLbq = cell(12,nvar);
pLbq = cell(12,nvar);

for j = 1:nvar %Loop over the variables
    hLillieTrend{j} = nan(12,31); %Fill the cells with data for each month and day
    hLillieDetrend{j} = nan(12,31);
    
    for i = 1:12 %Loop over the months
        
        % Data per month
        pos = find(dataDetrend(:,2)==i);
        if i == 2 %February!
            pos = find(dataDetrend(:,2)==i & dataDetrend(:,3)~=29);
        end
        
        % Plot (1), for each month all days will be plotted
        figure();
        plot(dataDetrend(pos,3),dataDetrend(pos,j+3),'b.');
        titlefig = sprintf('%s: x = days, y = detrended data for %s ',datatype{j},months{i});
        title(titlefig)
        
        % Normally distributed
        for k = 1:31
            if (i==2 && k~=29) %February
                if ~isempty(dataDayTrend{i,k})
                    dataMonthTrend{i,j} = [dataMonthTrend{i,j} dataDayTrend{i,k}(:,j)]; % Collects all daily data for each month, per variable
                    dataMonthDetrend{i,j} = [dataMonthDetrend{i,j} dataDayDetrend{i,k}(:,j)]; % Collects al daily detrended data for each month
                    hLillieTrend{j}(i,k) = lillietest(dataMonthTrend{i,j}(:,k),'alpha',0.001); % Performs Lilliefors goodness-of-fit test of composite normality on the data 
                    hLillieDetrend{j}(i,k) = lillietest(dataMonthDetrend{i,j}(:,k),'alpha',0.001); % Performs Lilliefors goodness-of-fit test of composite normality on the detrended data
                end
            else
                if (i~=2 && ~isempty(dataDayTrend{i,k}))
                    dataMonthTrend{i,j} = [dataMonthTrend{i,j} dataDayTrend{i,k}(:,j)]; % Collects all daily data for each month
                    dataMonthDetrend{i,j} = [dataMonthDetrend{i,j} dataDayDetrend{i,k}(:,j)]; % Collects al daily detrended data for each month
                    hLillieTrend{j}(i,k) = lillietest(dataMonthTrend{i,j}(:,k),'alpha',0.001);  % Performs Lilliefors goodness-of-fit test of composite normality on the data, returns 1 if the null hypothesis will be rejected 
                    hLillieDetrend{j}(i,k) = lillietest(dataMonthDetrend{i,j}(:,k),'alpha',0.001); % Performs Lilliefors goodness-of-fit test of composite normality on the detrended data, returns 1 if the null hypothesis will be rejected
                end   
            end
        end
            
        % Data with trend: homo- or heteroscedastic?
        y1 = 1:1:size(dataMonthTrend{i,j},2); %Uses the number of days
        y2 = repmat(y1,nyrs,1); %Retiled matrix: jaren * dagen 
        pVarTrend(j,i) = vartestn(dataMonthTrend{i,j},'TestType','BrownForsythe','display','off'); %Compares the different days for each month
        if sum(hLillieTrend{j}(i,:))~=0
            % Not-normally distributed
            if (pVarTrend(j,i) >= 0.001) %Probability under null hypothesis is large enough
                % Homoscedastic
                pTrend(j,i) = anova1(dataMonthTrend{i,j},y1,'off'); % anovatest, grouped per column
            else
                % Heteroscedastic
                pTrend(j,i) = kruskalwallis(dataMonthTrend{i,j},[],'off'); %Non-parametric test
                %tmp = [dataMonthTrend{i,j}(:) y2(:)]; 
                %pTrend(j,i) = welchAnova(tmp,0.001);
            end       
        else
            % Normally distributed
            if (pVarTrend(j,i) >= 0.001)
                % Homoscedastic
                pTrend(j,i) = anova1(dataMonthTrend{i,j},y1,'off');
            else
                % Heteroscedastic
                tmp = [dataMonthTrend{i,j}(:) y2(:)]; %colum 1: data, column 2: sample code
                pTrend(j,i) = welchAnova(tmp,0.001); % Test to compare the means of normally distributed populations 
            end
        end
        
        % Detrended data: homo- or heteroscedastic?
        y1 = 1:1:size(dataMonthDetrend{i,j},2); %Uses the number of days
        y2 = repmat(y1,nyrs,1); %Retiled matrix: jaren * dagen 
        pVarDetrend(j,i) = vartestn(dataMonthDetrend{i,j},'TestType','BrownForsythe','display','off'); %Compares the different days for each month    
        if sum(hLillieDetrend{j}(i,:))~=0
            % Not-normal distributed
            if (pVarDetrend(j,i) >= 0.001) %Probability under null hypothesis is large enough
                % Homoscedastic
                pDetrend(j,i) = anova1(dataMonthDetrend{i,j},y1,'off'); % anovatest, grouped per column
            else
                % Heteroscedastic
                pDetrend(j,i) = kruskalwallis(dataMonthDetrend{i,j},[],'off'); %Non-parametric test
                %tmp = [dataMonthTrend{i,j}(:) y2(:)];
                %pDetrend(j,i) = welchAnova(tmp,0.001);
            end
        else
            % Normal distributed
            if (pVarDetrend(j,i) >= 0.001)
                % Homoscedastic
                pDetrend(j,i) = anova1(dataMonthDetrend{i,j},y1,'off');
            else
                % Heteroscedastic
                tmp = [dataMonthTrend{i,j}(:) y2(:)]; %colum 1: data, column 2: sample code
                pDetrend(j,i) = welchAnova(tmp,0.001); % Test to compare the means of normally distributed populations 
            end
        end
        
        % Autocorrelation?
        tmp = dataDetrend(pos,j+3);
        maxday = max(dataDetrend(pos,3));
        for k = 1:nyrs
            [hLbq{i,j}(k,:), pLbq{i,j}(k,:)] = lbqTest(tmp((k-1)*maxday+1:(k*maxday)), 'lags', [1:5], 'alpha', 0.001); %Ljung-box Q-test for residual autocorrelation
        end
        
    end
    
    if strcmp(datatype{j},'P') == 1
        outputfile = sprintf('statsData_%s_%s.mat',modelname, 'P');
    elseif strcmp(datatype{j}, 'E') == 1
        outputfile = sprintf('statsData_%s_%s.mat',modelname, 'E');
    elseif strcmp(datatype{j},'T') == 1
        outputfile = sprintf('statsData_%s_%s.mat',modelname, 'T');
    end
    
    %Out per variable
    out.dataDetrend = [dataDetrend(:,1:3) dataDetrend(:,3+j)];
    out.hLillieTrend = hLillieTrend{j};
    out.hLillieDetrend = hLillieDetrend{j};
    out.pTrend = pTrend(j,:);
    out.pDetrend = pDetrend(j,:);
    out.pVarTrend = pVarTrend(j,:);
    out.pVarDetrend = pVarDetrend(j,:);
    out.hLbq = hLbq(:,j);
    out.pLbq = pLbq(:,j);
    save(strcat(save_loc, outputfile),'out');
    
    close all;
end

end

