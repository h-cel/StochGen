function statistics=CalMomentsT(data,rep,agg,orig,zdpflag, fignumber)
%% CalMomentsT Function to calculate the moments of T and to make a plot 
%
%   Inputs:
%       data: stochastically generated time series
%       rep: number of repetitions used
%       agg: vector with aggregation levels to be used
%       orig: original time series used for calibration of the vine copula
%       model
%       zdpflag: Boolean to indicate how zdp has to be calculated (1: dry
%       days are all days <0.1 mm, 0: dry days are all days with 0 mm)
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

%% Statistics stochastically generated time series

nragg=length(agg);
nrstat=5;
statistics=nan(rep,12,nragg,nrstat);

for k = 1:rep
    p = data(:,[1:3 k+3]);
    r = size(p,1)*24;
    dataUcclestyle = zeros(r,13);
    ymd = p(:,1:3);
    iDay = 0:24:r;
    dataUcclestyle(:,1) = ones(r,1)*6;
    dataUcclestyle(:,2) = ones(r,1)*6100;
    dataUcclestyle(:,7) = ones(r,1);
    for j = 1:size(p,1)
        dataUcclestyle(iDay(j)+1:iDay(j+1),3:5) = repmat(ymd(j,:),24,1);
        dataUcclestyle(iDay(j)+1:iDay(j+1),6) = 0:23;
        dataUcclestyle(iDay(j)+1,8) = p(j,end);
    end
    dataUcclestyle((dataUcclestyle(:,8) < 0.1),8) = 0;
    
    tmp=StatCalcm3(dataUcclestyle,agg,zdpflag);
    for m=1:12
        statistics(k,m,:,:)=tmp(:,:,m)';
    end
end

%% Statistics original time series

tmp=StatCalcm3(orig,agg,zdpflag);
statorig=zeros(1,12,nragg,nrstat);
for m=1:12
    statorig(1,m,:,:)=tmp(:,:,m)';
end

%% Plots

font=14;
t=0;
figure_handle = figure(fignumber);
set(gcf,'color','white')
ylabels={'Mean [°C]','Variance [°C²]','Autocovariance [°C²]'};
for i = 1:3
    set(gca,'FontSize',font);
    for a = 1:nragg
        t = t+1;
        subplot(2,3,t)
        boxplot(statistics(:,:,a,i),'color',[0.5 0.5 0.5],'symbol','.k');hold on
        plot(1:12,statorig(:,:,a,i),'*','linewidth',1.25,'color',[51/256 153/256 1]);
        hold on
        ylabel(ylabels{i},'FontSize',18);
        set(gca,'FontSize',font);
        
        if i == 1
            ylim([0 28])
        elseif i == 2
            ylim([5 40])
            % Legend support plots
            subplot(2,3,t+3)
            plot(-100,-100,'*','linewidth',1.25,'color',[51/256 153/256 1]);
            hold on
            plot(-10:-1,-10:-1,'-','color',[0.5 0.5 0.5]);
            ylim([5 100])
            axis off
            h = legend({'Original time series', 'Stochastically generated time series'},'FontSize',15);
            legend('Location','northoutside','Orientation','horizontal')
            legend('boxoff')
        else
            ylim([0 20])
        end
        xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
        
        hold on
        xlabel('Maand','FontSize',font);
       
    end
end
set(gca,'fontsize',14)
fullfig(figure_handle)

end

