function stats = MomentsE(file, c, rep, timefile, casename, fignumber) 
%% MomentsE Calculates model performance E-generator
%   This function calculates the performance of the E generator for a case
%   of choice.
%
%   Inputs:
%       file: string with the basis file name (refers to both the original
%       file and the simulations based on this file)
%       c: case that has to be loaded
%       rep: number of repeats
%       timedata: time file necessary
%       casename: indicating the specific case that has to be loaded for
%       case 3
%       fignumber: number of the figure, to allow for multiple figures to
%       be plotted in the same script.
%
%   Outputs:
%       stats: statistics
%
    %   Last updated by J. Van de Velde on 10/06/'21: fignumber added

%% Setup

addpath(genpath('D:\Users\jpvdveld\Documents\PhD\Code\StochasticModelling'), genpath('E:\Users\jpvdveld\Onderzoek\Data')) %Both Code and Data paths need to be added with their subfolders.

if iscell(c) == 1
    c = c{1}; % Only works if there is only one assigned case
end

%% Loading Data

time = matload(timefile);

if any(strcmp(c, 'c1'))
        StochGen = matload([file, '_Esim_c1.mat']);
        [~, C] = size(StochGen);
        Ecase = cell(1, C-3);
        for n = 4:rep+3
            Ecase{n-3} = StochGen(:,n);
        end

        
end

if any(strcmp(c, 'c2'))
    StochGen = matload([file, '_Esim_c2.mat']);
    Ecase = cell(1, rep);
    rand1 = randi(rep,[1,rep]);
    rand2 = randi(rep,[1,rep]);
    for j = 1:rep
        Ecase{j} = StochGen{rand1(j)}(:,rand2(j)+3);
    end

end

if any(strcmp(c, 'c3'))
    StochGen = matload([file, '_Esim', casename, '_c3.mat']);
    Ecase = cell(1, rep);
    rand1 = randi(rep,[1,rep]);
    rand2 = randi(rep,[1,rep]);
    rand3 = randi(rep,[1,rep]);
    for k = 1:rep
        Ecase{k} = StochGen{rand1(k),rand2(k)}(:,rand3(k)+3);
    end
end

ClimSim = matload([file, '.mat']);

if isnumeric(ClimSim) == 1
    e_orig = rcm2uccle(ClimSim(:,4), time, '24');
else
    e_orig = rcm2uccle(ClimSim{1}(:,4),time, '24');
end

%% Calculation

e1full = [time cell2mat(Ecase)];
stats = CalMomentsE(e1full,rep,24,e_orig,1, fignumber);

end



 