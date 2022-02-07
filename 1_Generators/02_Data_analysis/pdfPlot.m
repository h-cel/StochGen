function score = pdfPlot(file, timefile, casenumber, casename, fignumber)
%% pdfPlot Plots the pdfs for all cases
%
%   This function draws pdfs (based on the historgram- for the different
%   variables in StochGen: E (case 1-2-3), T (case 2-3) and P (case 3). In
%   addition, the Perkins' Skill Score is calculated for a quantitative
%   analysis.
%
%   Inputs:
%       file: the filename to be loaded
%       timefile: the timefile corresponding to the file (which contains
%       the time values of the input time series).
%       casenumber: cell array with cases that have to be compared
%       casename: for case 3, the specific name (e.g. time series length)
%       that has been run
%       fignumber: number of the figure, to allow for multiple figures to
%       be plotted in the same script.
%
%   Output:
%       score: Perkins' Skill Score
%
%   Last updated by J. Van de Velde on 06/10/'21: fignumber added

%% Setup 
% General

addpath(genpath('D:\Users\jpvdveld\Documents\PhD\Code\StochasticModelling'), genpath('E:\Users\jpvdveld\Onderzoek\Data')) %Both Code and Data paths need to be added with their subfolders.

% Data
time = matload(timefile);

% Score
score = zeros(3,1);

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


%% Case 1

if any(strcmp(casenumber, 'c1'))
    
    % Setup
    
    StochGenEc1 = matload([file,'_Esim_c1.mat']);
    
end

%% Case 2

if any(strcmp(casenumber, 'c2'))
    
    % Setup
    
    StochGenEc2 = matload([file,'_Esim_c2.mat']);
    StochGenTc2 = matload([file,'_Tsim_c2.mat']);
    
    StochGenEc2All = NaN(length(StochGenEc2{1}), 400);
    
    
    % Ordering
    
    t=0;
    for i = 1:20
        for j = 1:20
            t = t+1;
            StochGenEc2All(:,t) = StochGenEc2{i}(:,j+3);
        end
    end
    
end

%% Case3

if any(strcmp(casenumber, 'c3'))
    
    % Setup
    
    StochGenEc3 = matload([file,'_Esim', casename, '_c3.mat']);
    StochGenTc3 = matload([file,'_Tsim', casename, '_c3.mat']);
    StochGenPc3 = matload([file,'_Psim', casename, 'daily_c3.mat']);
    
    StochGenTc3All = NaN(length(StochGenTc3{1}), 400);
    StochGenEc3All = NaN(length(StochGenTc3{1}), 8000);
    
    % Ordering
    
    t=0;
    y=0;
    for k = 1:20
        for i = 1:20
            y=y+1;
            StochGenTc3All(:,y) = StochGenTc3{k}(:,i+3);
            
            for j = 1:20
                t = t+1;
                StochGenEc3All(:,t) = StochGenEc3{k,i}(:,j+3);
            end
        end
    end
    
end


%% Plots

figure_handle = figure(fignumber);

if any(strcmp(casenumber, 'c1'))
    % C1 - E
    StochGenEc1select=zeros(length(StochGenEc1),20);
    subplot(3,3,1)
    hold on
    for i = 1:20
        StochGenEc1select(:,i) = StochGenEc1(:, i+3);
        histogram(StochGenEc1select(:,i), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', 'red');
    end
    histogram(e_orig, 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', 'blue');
    title('E', 'FontWeight','Normal','FontSize',15)
    ylabel('Case 1','Fontsize', 15)
    
    % C1 - T: is not simulated
    subplot(3,3,2)
    axis off
    title('T', 'FontWeight','Normal','FontSize',15)
    
    % C1 - P: is not simulated
    subplot(3,3,3)
    axis off
    title('P', 'FontWeight','Normal','FontSize',15)
    
end

if any(strcmp(casenumber, 'c2'))
    
    % C2 - E: 20 values randomly selected from 400
    StochGenEc2select=zeros(length(StochGenEc2All),20);
    r = randi([1 400], [1 20]);
    subplot(3,3,4)
    hold on
    for i = 1:20
        StochGenEc2select(:,i) = StochGenEc2All(:, r(i));
        histogram(StochGenEc2select(:,i), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', 'red');
    end
    histogram(e_orig, 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', 'blue');
    ylabel('Case 2','Fontsize', 15)
    
    % C2 - T
    StochGenTc2select=zeros(length(StochGenTc2),20);
    subplot(3,3,5)
    hold on
    for i = 1:20
        StochGenTc2select(:,i) = StochGenTc2(:, i+3);
        histogram(StochGenTc2(:, i+3), 'BinWidth', 5, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', 'red')
    end
    histogram(t_orig, 'BinWidth', 5, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', 'blue');
    
    % C2 - P: not simulated
    subplot(3,3,6)
    axis off
    
end

if any(strcmp(casenumber, 'c3'))
    
    % C3 - E: 20 values randomly selected from 800
    StochGenEc3select=zeros(length(StochGenEc3All),20);
    r = randi([1 8000], [1 20]);
    subplot(3,3,7)
    title('E', 'FontWeight','Normal','FontSize',15)
    hold on
    for i = 1:20
        StochGenEc3select(:,i) = StochGenEc3All(:, r(i));
        histogram(StochGenEc3select(:,i), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', 'red')
    end
    histogram(e_orig, 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', 'blue')
    ylabel('Case 3','Fontsize', 15)
    
    % C3 - T: 20 values randomly selected from 400
    StochGenTc3select=zeros(length(StochGenTc3All),20);
    r = randi([1 400], [1 20]);
    subplot(3,3,8)
    title('T', 'FontWeight','Normal','FontSize',15)
    hold on
    for i = 1:20
        StochGenTc3select(:,i) = StochGenTc3All(:, r(i));
        histogram(StochGenTc3select(:,i), 'BinWidth', 5, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', 'red')
    end
    histogram(t_orig, 'BinWidth', 5, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', 'blue')
    
    % C3 - P
    StochGenPc3select=zeros(length(StochGenPc3),20);
    subplot(3,3,9)
    title('P', 'FontWeight','Normal','FontSize',15)
    hold on
    for i = 1:20
        StochGenPc3select(:,i) = StochGenPc3(:, i+3);
        histogram(StochGenPc3select(:,i), 'BinWidth', 5, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', 'red')
    end
    histogram(StochGenPc3(:, i+3), 'BinWidth', 5, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', 'blue')

    score(1) = Perkins(StochGenEc3select, e_orig, 20);
    score(2) = Perkins(StochGenTc3select, t_orig, 20);
    score(3) = Perkins(StochGenPc3select, p_orig, 20);
    
end

fullfig(figure_handle)

end