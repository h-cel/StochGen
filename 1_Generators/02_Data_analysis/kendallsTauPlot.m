function ktvalues = kendallsTauPlot(file, timefile, casenumber, casename, rep, fignumber)
%% KendallstauPlot Calculates Kendall's tau for the different cases and combinations of variables
%
%   This function calculates Kendall's tau values for the different
%   combinations of variables in StochGen, i.e. PT, ET, PE. These values
%   are, when relevant, plotted.
%
%   Inputs:
%       file: the filename to be loaded
%       timefile: the timefile corresponding to the file (which contains
%       the time values of the input time series).
%       casenumber: cell array indicating for which cases the Kendall's tau
%       values have to be compared
%       casename: string, allows to specify the way case 3 was run
%       rep: number of repetitions used
%       fignumber: number of the figure, to allow for multiple figures to
%       be plotted in the same script.

%   Output:
%       ktvalues: values of Kendall's tau for each combination that
%       contains at least 1 generated variable.
%
%   Last updated by J. Van de Velde on 10/06/'21: fignumber added

%% Setup 
% General

addpath(genpath('D:\Users\jpvdveld\Documents\PhD\Code\StochasticModelling'), genpath('E:\Users\jpvdveld\Onderzoek\Data')) %Both Code and Data paths need to be added with their subfolders.

% Data
time = matload(timefile);

%% Reference

orig = matload([file, '.mat']);
if isa(orig, 'cell') == 1
    e_orig = orig{1}(:,4);
    t_orig = orig{1}(:,5);
    p_orig = orig{1}(:,6);
    time_orig = orig{1}(:,1:3);
else
    e_orig = orig(:,4);
    t_orig = orig(:,5);
    p_orig = orig(:,6);
    time_orig = orig(:,1:3);
end


 
pt_orig = kendallsTau(time_orig,p_orig,t_orig);
pe_orig = kendallsTau(time_orig,p_orig,e_orig);
et_orig = kendallsTau(time_orig,e_orig,t_orig);

%% File preparation

pec1 = NaN(rep,12);
etc1 = NaN(rep,12);
ptc2 = NaN(rep,12);
pec2 = NaN(rep^2, 12);
etc2 = NaN(rep^2, 12);
ptc3 = NaN(rep^2,12);
pec3 = NaN(rep^3, 12);
etc3 = NaN(rep^3, 12);
    

%% Case 1

if any(strcmp(casenumber, 'c1'))
    
    % Setup
    
    StochGenE = matload([file,'_Esim_c1.mat']);
    
    % Calculation
    
    for i = 1:rep
        pec1(i,:) = kendallsTau(time,p_orig,StochGenE(:, i+3));
        etc1(i,:) = kendallsTau(time,StochGenE(:, i+3),t_orig);
    end
    
end

%% Case 2

if any(strcmp(casenumber, 'c2'))
    
    % Setup
    
    StochGenE = matload([file,'_Esim_c2.mat']);
    StochGenT = matload([file,'_Tsim_c2.mat']);
    
    % Calculation
    
    t=0;
    r1 = randi([1 rep], [1 rep]);
    r2 = randi([1 rep], [1 rep]);
    for i = 1:rep
        ptc2(i,:) = kendallsTau(time,p_orig, StochGenT(:,i+3));
        
        for j = 1:rep
            t = t+1;
            pec2(t,:) = kendallsTau(time,p_orig,StochGenE{r1(i)}(:,r2(j)+3));
            etc2(t,:) = kendallsTau(time,StochGenE{r1(i)}(:,r2(j)+3),StochGenT(:,r1(i)+3));
        end
    end
    
end

%% Case3

if any(strcmp(casenumber, 'c3'))
    
    % Setup
    
    StochGenE = matload([file,'_Esim' casename '_c3.mat']);
    StochGenT = matload([file,'_Tsim' casename '_c3.mat']);
    StochGenP = matload([file,'_Psim100daily_c3.mat']);
    
    % Calculation
    
    r1 = randi([1 rep], [1 rep]);
    r2 = randi([1 rep], [1 rep]);
    r3 = randi([1 rep], [1 rep]);
    
    t=0;
    y=0;
    for k = 1:rep
        for i = 1:rep
            y=y+1;
            ptc3(y,:) = kendallsTau(time,StochGenP(:,r1(k)+3), StochGenT{r1(k)}(:,r2(i)+3));
            
            for j = 1:rep
                t = t+1;
                pec3(t,:) = kendallsTau(time,StochGenP(:,r1(k)+3),StochGenE{r1(k),r2(i)}(:,r3(j)+3));
                etc3(t,:) = kendallsTau(time,StochGenE{r1(k),r2(i)}(:,r3(j)+3),StochGenT{r1(k)}(:,r2(i)+3));
                
            end
        end
    end
    
end

%% Save values

ktvalues = {NaN, ptc2, ptc3; pec1, pec2, pec3; etc1, etc2, etc3};

%% Plots

figure_handle = figure(fignumber);

if any(strcmp(casenumber, 'c1'))
    
    % C1 - PT: does not exist in simulations
    subplot(3,3,1)
    axis off
    title('Case 1', 'FontWeight','Normal','FontSize',15)
    ylabel('PT','Fontsize', 15)
    
    % C1 - PE
    subplot(3,3,4)
    hold on
    plot(1:12,pe_orig,'*','linewidth', 1.5,'Color',[51/256 153/256 1])
    boxplot(pec1,'Color',[0.7 0.7 0.7],'symbol','+')
    plot(0:13,zeros(1,14),'Color',[0.5 0.5 0.5])
    xticks([1:12])
    xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
    ylabel({'Kendall''s tau'; 'PE'},'FontSize',15)
    ylim([-0.6 0.6])
    set(gca,'fontsize',14)
    
    % C1 - ET
    subplot(3,3,7)
    hold on
    plot(1:12,et_orig,'*','linewidth', 1.5,'Color',[51/256 153/256 1])
    boxplot(etc1,'Color',[0.7 0.7 0.7],'symbol','+')
    plot(0:13,zeros(1,14),'Color',[0.5 0.5 0.5])
    xticks([1:12])
    xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
    ylabel('ET','Fontsize', 15)
    ylim([-0.6 0.6])
    set(gca,'fontsize',14)
    
end

if any(strcmp(casenumber, 'c2'))
    
    % C2 - PT
    subplot(3,3,2)
    hold on
    plot(1:12,pt_orig,'*','linewidth', 1.5,'Color',[51/256 153/256 1])
    boxplot(ptc2,'Color',[0.7 0.7 0.7],'symbol','+')
    plot(0:13,zeros(1,14),'Color',[0.5 0.5 0.5])
    xticks([1:12])
    xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
    title('Case 2', 'FontWeight','Normal','FontSize',15)
    ylabel('PT','Fontsize', 15)
    ylim([-0.6 0.6])
    set(gca,'fontsize',14)
    
    % C2 - PE: 20 values randomly selected from 400
    subplot(3,3,5)
    hold on
    plot(1:12,pe_orig,'*','linewidth', 1.5,'Color',[51/256 153/256 1])
    boxplot(pec2,'Color',[0.7 0.7 0.7],'symbol','+')
    plot(0:13,zeros(1,14),'Color',[0.5 0.5 0.5])
    xticks([1:12])
    xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
    ylim([-0.6 0.6])
    set(gca,'fontsize',14)
    
    % C2 - ET: 20 values randomly selected from 400
    subplot(3,3,8)
    hold on
    plot(1:12,et_orig,'*','linewidth', 1.5,'Color',[51/256 153/256 1])
    boxplot(etc2,'Color',[0.7 0.7 0.7],'symbol','+')
    plot(0:13,zeros(1,14),'Color',[0.5 0.5 0.5])
    xticks([1:12])
    xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
    %xlabel('Maand')
    ylim([-0.6 0.6])
    set(gca,'fontsize',14)
    
end

if any(strcmp(casenumber, 'c3'))
    
    % C3 - PT: 20 values randomly selected from 400
    subplot(3,3,3)
    hold on
    plot(1:12,pt_orig,'*','linewidth', 1.5,'Color',[51/256 153/256 1])
    boxplot(ptc3,'Color',[0.7 0.7 0.7],'symbol','+')
    plot(0:13,zeros(1,14),'Color',[0.5 0.5 0.5])
    xticks([1:12])
    xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
    title('Case 3', 'FontWeight','Normal','FontSize',15)
    ylim([-1 1])
    set(gca,'fontsize',14)
    
    % C3 - PE: 20 values randomly selected from 8000
    subplot(3,3,6)
    hold on
    plot(1:12,pe_orig,'*','linewidth', 1.5,'Color',[51/256 153/256 1])
    boxplot(pec3,'Color',[0.7 0.7 0.7],'symbol','+')
    plot(0:13,zeros(1,14),'Color',[0.5 0.5 0.5])
    xticks([1:12])
    xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
    ylim([-1 1])
    set(gca,'fontsize',14)
    
    % C3 - ET: 20 values randomly selected from 8000
    subplot(3,3,9)
    hold on
    plot(1:12,et_orig,'*','linewidth', 1.5,'Color',[51/256 153/256 1])
    boxplot(etc3,'Color',[0.7 0.7 0.7],'symbol','+')
    plot(0:13,zeros(1,14),'Color',[0.5 0.5 0.5])
    xticks([1:12])
    xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
    ylim([-1 1])
    set(gca,'fontsize',14)
    
end

fullfig(figure_handle)

end