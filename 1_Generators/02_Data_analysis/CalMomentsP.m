function statistics=CalMomentsP(datafile,rep,aggreq,origfile,zdpflag, aggdata, fignumber, locP)
%% CalMomentsP Function to calculate the moments of P and to make a plot of these
%   moments, based on the aggregation level.
%
%   Inputs:
%       datafile: stochastically generated time series
%       rep: number of repetitions used
%       aggreq: vector with aggregation levels to be used
%       origfile: original time series used for calibration of BL-model
%       zdpflag: Boolean to indicate how zdp has to be calculated (1: dry
%       days are all days <0.1 mm, 0: dry days are all days with 0 mm)
%       aggdata: aggregation of the data
%       locP: if the orig file contains multiple variables, this indicates
%       the location
%       fignumber: number of the figure, to allow for multiple figures to
%       be plotted in the same script.
%
%   Outputs:
%       statistics: 4-D matrix with the statistcs
%
%   Last updated by J. Van de Velde at 10/06/'21: fignumber added
%   Other contributors to earlier versions of this file:
%       M. T. Pham
%       K. De Roos
%       H. Vernieuwe?

%% Set-up

addpath(genpath('E:\Users\jpvdveld\Onderzoek\Data'))

%% Loading & preprocessing

data = matload(datafile);
orig = matload(origfile);

if iscell(orig)
    orig=orig{1};
end

uniqyrs = unique(orig(:,1));
meanP = zeros(length(uniqyrs),1);

for y= 1:length(uniqyrs)
    datasep = orig(orig(:,1) == uniqyrs(y) & orig(:,2) == 9,6);
    meanP(y) = mean(datasep);
end

orig(orig(:,1) == uniqyrs(meanP == max(meanP)) & orig(:,2) == 9,6) = mean(meanP); 

%% Statistics stochastically generated time series

% 
nragg=length(aggreq);
nrstat=5;
statistics=nan(rep,12,nragg,nrstat);

for k = 1:rep
    switch aggdata %Transforms according to given aggregate level
        case '24'
            p = data(:,[1:3 k+3]);
            time = p(:,1:3);
            dataUcclestyle = rcm2uccle(p,time,aggdata);
        case '1'
            p = data(:,[1:4 k+4]);
            time = p(:,1:4);
            dataUcclestyle = rcm2uccle(p,time,aggdata);
    end

    tmp=StatCalcm3(dataUcclestyle,aggreq,zdpflag);
    for m=1:12
        statistics(k,m,:,:)=tmp(:,:,m)';
    end
end

%% Statistics original time series

if nargin == 8
    switch aggdata
        case '1'
            obs=rcm2uccle(orig(:,[1:4,locP]), orig(:,1:4), aggdata);
        case '24'
            obs=rcm2uccle(orig(:,[1:3,locP]), orig(:,1:3), aggdata);
    end
else
    switch aggdata
        case '1'
            obs=rcm2uccle(orig, orig(:,1:4), aggdata);
        case '24'
            obs=rcm2uccle(orig, orig(:,1:3), aggdata);
    end
end

tmp=StatCalcm3(obs,aggreq,zdpflag);
statobs=zeros(1,12,nragg,nrstat);
for m=1:12
    statobs(1,m,:,:)=tmp(:,:,m)';
end

%% Plots

font=14;
t=0;
figure_handle = figure(fignumber);
set(gcf,'color','white')
ylabels={'Mean [mm]','Variance [mm^2]','Autocovariance [mm^2]','ZDP', 'Third central moment'};
for i=1:nrstat
    set(gca,'FontSize',font);
    for a=1:nragg
        
       % Plot
        t=t+1;
        subplot(5,nragg,t)
        boxplot(statistics(:,:,a,i),'color',[0.5 0.5 0.5],'symbol','.k');
        hold on
        plot(1:12,statobs(:,:,a,i),'*','linewidth',1.25,'color',[51/256 153/256 1]);
        
        % Labels
        
        if rem((t-1), nragg) == 0
            ylabel(ylabels{i},'FontSize',font);
        end
        
        if t>nragg*4
            xlabel('Month','FontSize',font);
        end
        
        set(gca,'FontSize',font);
        
        % Title
        
        if i == 1
            if aggreq(a)==1/6
                title('10 min','FontSize',font,'FontWeight','Normal');
            else
                title([num2str(aggreq(a)) ' h'],'FontWeight','Normal');
            end
        end
        
        % Limits
        
        maxdata = max(max(statistics(:,:,a,i)));
        maxobs =  max(max(statobs(:,:,a,i)));
        
        mindata = min(min(statistics(:,:,a,i)));
        minobs =  min(min(statobs(:,:,a,i)));
        
        maxplot = max(maxdata,maxobs);
        minplot = min(mindata, minobs);
        
        ylim([floor(minplot),ceil(maxplot)]);
            
        xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
        
    end
end
fullfig(figure_handle)

end

