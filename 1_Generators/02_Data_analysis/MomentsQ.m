function stats =  MomentsQ(names, c, rep, timefile, casename, fignumber) 
%% MomentsQ Calculates model performance Q
%   This function calculates the performance of the Q for a case of choice
%
%   Inputs:
%       names: strings with the basis file names (refers to both the original
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

%% Set-up
 
addpath(genpath('D:\Users\jpvdveld\Documents\PhD\Code\StochasticModelling'), genpath('E:\Users\jpvdveld\Onderzoek\Data')) %Both Code and Data paths need to be added with their subfolders.

%% Loading data

% basic files

PDM_par = matload('paramPDM.mat');
inputs.A = 385;

time = matload(timefile);

% Original Q

ClimSim = matload([names{1}, '.mat']);
if isnumeric(ClimSim) == 1
    e_orig = ClimSim(:,4);
    p_orig = ClimSim(:,6);
else
    e_orig = ClimSim{1}(:,4);
    p_orig = ClimSim{1}(:,6);
    p_orig = RemovePOutliers(p_orig,names);
end 

ndaysOrig = length(p_orig);
inputs.P = reshape(kron(p_orig/24, ones(1,24))', ndaysOrig*24,1); %Hourly precipitation input data
inputs.E = reshape(kron(e_orig/24, ones(1,24))', ndaysOrig*24,1); %Hourly evaporation input data

Q = PDMPieter(inputs, PDM_par);
q_orig = rcm2uccle(Q, time, '24');

% Case 1

if any(strcmp(c, 'c1'))
    StochGen = matload([names{1}, '_Esim_c1.mat']);
    Qcase = cell(1, rep);
    for n = 4:rep+3
        E = StochGen(:,n);
        inputs.P = reshape(kron(p_orig/24, ones(1,24))', ndaysOrig*24,1); %Hourly precipitation input data
        inputs.E = reshape(kron(E/24, ones(1,24))', ndaysOrig*24,1); %Hourly evaporation input data
        Qcase{n-3} = PDMPieter(inputs, PDM_par);
    end

end

if any(strcmp(c, 'c2'))
    StochGen = matload([names{1}, '_Esim_c2.mat']);
    Qcase = cell(1, rep);
    rand1 = randi(rep,[1,rep]);
    rand2 = randi(rep,[1,rep]);
    for j = 1:rep
        E = StochGen{rand1(j)}(:,rand2(n)+3);
        inputs.P = reshape(kron(p_orig/24, ones(1,24))', ndaysOrig*24,1); %Hourly precipitation input data
        inputs.E = reshape(kron(E/24, ones(1,24))', ndaysOrig*24,1); %Hourly evaporation input data
        Qcase{j} = PDMPieter(inputs, PDM_par);
    end
end

if any(strcmp(c, 'c3'))
    StochGenE = matload([names{1}, '_Esim', casename, '_c3.mat']);
    StochGenP =  matload([names{1}, '_Psim100daily_c3.mat']);
    Qcase = cell(1, rep);
    ndays = length(StochGenP);
    rand1 = randi(rep,[1,rep]);
    rand2 = randi(rep,[1,rep]);
    rand3 = randi(rep,[1,rep]);
    for k = 1:rep
        E = StochGenE{rand1(k),rand2(k)}(:,rand3(k));
        P = StochGenP(:,rand3(k));
        inputs.P = reshape(kron(P/24, ones(1,24))', ndays*24,1); %Hourly precipitation input data
        inputs.E = reshape(kron(E/24, ones(1,24))', ndays*24,1); %Hourly evaporation input data
        Qcase{k} = PDMPieter(inputs, PDM_par);
    end
end



%% Calculations

Q1full = [time cell2mat(Qcase)];
stats = CalMomentsQ(Q1full,rep,24,q_orig,1, fignumber);

end


 