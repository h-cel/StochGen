function stats =  MomentsT(file, c, rep, timefile, casename, fignumber) 
%% MomentsT Calculates model performance T-generator
%   This function calculates the performance of the T generator for a case
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
%       stats: statistics (4D-matrix)
%
%   Last updated by J. Van de Velde on 10/06/'21: fignumber added

%% Setup
 
addpath(genpath('D:\Users\jpvdveld\Documents\PhD\Code\StochasticModelling'), genpath('E:\Users\jpvdveld\Onderzoek\Data')) %Both Code and Data paths need to be added with their subfolders.

%% Data inladen

time = matload(timefile);

if any(strcmp(c, 'c2'))
    StochGen = matload([file,'_Tsim_c2.mat']);
    Tcase = cell(1, rep);
    for n = 4:rep+3
        Tcase{n-3} = StochGen(:,n);
    end

end

if any(strcmp(c, 'c3'))
    StochGen = matload([file,'_Tsim' casename, '_c3.mat']);
    Tcase = cell(1, rep);
    rand1 = randi(rep,[1,rep]);
    rand2 = randi(rep,[1,rep]);
    for j = 1:rep
        Tcase{j} = StochGen{rand1(j)}(:,rand2(j)+3);
    end
end

ClimSim = matload([file, '.mat']);
if isnumeric(ClimSim) == 1
    t_orig = rcm2uccle(ClimSim(:,5), time, '24');
else
    t_orig = rcm2uccle(ClimSim{1}(:,5),time, '24');
end

%% Calculations

t1full = [time cell2mat(Tcase)];
stats = CalMomentsT(t1full,rep,24,t_orig,1, fignumber);

end


 