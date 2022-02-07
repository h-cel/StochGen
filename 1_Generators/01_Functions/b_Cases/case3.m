%% Case 3    P via MBL toolbox -> T simulation -> E simulation -> Q simulation
%   Loads Case 3 of the stochastic generator chain. All data are loaded in
%   the LaunchCases function
%
%   Last updated by J. Van de Velde on 04/05/'21: added RemovePOutliers

%% Set-up

%Locations
loc_load = 'E:\Users\jpvdveld\Onderzoek\Data\1_biascorrection\';
loc_save = 'E:\Users\jpvdveld\Onderzoek\Data\3_cases\';

%% Initialization

outfile = [names{1} '.mat'];
data_tmp = matload(strcat(loc_load, outfile));
if isa(data_tmp, 'cell') == 1
    data_tmp = data_tmp{1};
end

data_tmp = RemovePOutliers(data_tmp, names);

time = CreateTimeMatrix(startdate, enddate);

vars = vars(1:3);

%% Simulations of P

[Psim, CalStr] = simBLBOX('RBL2', 7, data_tmp, vars, startdate, enddate, repeat, agg, calstr); %Outputs: simulated P and a calibration string

[Psimtmp,PnoiseMBL] = postprocessBLBOX(time, repeat, Psim, agg);

save([loc_save names{1} '_' names{6} '_c3.mat'], 'Psim', '-v7.3');
save([loc_save names{1} '_' names{7} '_c3.mat'], 'Psimtmp', '-v7.3');
save([loc_save names{1} '_', names{9} '.mat'], 'CalStr');
clear Psimtmp Psim CalStr

%% Simulations of T using stochastic P
%Initialization

Tsim = cell(1,repeat);
Tini = T(1);
%Loop over repeats of T
for i = 1:repeat
    Prep = PnoiseMBL(:,i+3);
    Tsim{i}(:,1:3) = time;
    %Loop over repeats of P (20*20 -> 400 repeats)
    for k = 1:repeat
        Tsim{i}(:,k+3) = generateT(Prep,RT,thetasT,familyTpPT,time,Tini,dataMonthT, copMatrixT);
        fprintf('T %d/%d of repetition %d/%d simulated.\n', k, repeat,i,repeat);
    end
end
disp('T case 3 simulated.')

%% Simulations of E using stochastic P and stochastic T
%Initialization
Eini = E(1);
Esim = cell(repeat,repeat);
%Loop over repeats of P
for i = 1:repeat
    Prep = PnoiseMBL(:,i+3);
    %Loop over repeats of T
    for k = 1:repeat
        Trep = Tsim{i}(:,k+3);
        Esim{i,k}(:,1:3) = time;
        %loop over repeats of E (20*20*20 -> 8000)
        for j = 1:repeat
            Esim{i,k}(:,j+3) =  generateE(Trep,Prep,RE,thetasE,familyTPEpE,time,Eini,dataMonthE, copMatrixE);
            fprintf('E repetition %d/%d of %d/%d of %d/%d simulated.\n', j, repeat, k, repeat,i,repeat);
         end
    end
end
disp('E case 3 simulated.')

%% Destandardize E and T
%Destandardization of E and T is seperated because these differ in repeat depth 
Tsim2 = cell(size(Tsim));
Esim2 = cell(size(Esim));
for i = 1:repeat
    for k = 1:repeat
        % Destandardize T
        for n = 1:length(vars)
            if strcmp(vars{n}, 'T') == 1
                for d = 1:length(time)
                    %Selection of time
                    md = time(d,2);
                    dd = time(d,3);
                    pos = find((dataStand(:,2) == md & dataStand(:,3) == dd));
                    Tsim2{i}(d,k+3) = Tsim{i}(d,k+3).*dataStand(pos(1),n+6) + dataStand(pos(1),n+3);
                    
                end
            end
        end
        
        % Destandardize E
        for j = 1:repeat
            for n = 1:length(vars)
                if strcmp(vars{n}, 'E') == 1
                    for d = 1:length(time)
                        %Selection of time
                        md = time(d,2);
                        dd = time(d,3);
                        pos = find((dataStand(:,2) == md & dataStand(:,3) == dd));
                        Esim2{i,k}(d,j+3) = Esim{i,k}(d,j+3).*dataStand(pos(1),n+6) + dataStand(pos(1),n+3);
                        
                    end
                end
            end
        end
        
    end
    Tsim2{i}(:,1:3) = time;
    Esim2{i,k}(:,1:3) = time;
end



% Save

save([loc_save names{1} '_' names{4} '_c3.mat'], 'Tsim2', '-v7.3');
clear Tsim2

%% Simulation of Q using stochastic P and stochastic E

%Initialization
ndays = length(time(:,1));
Qsim = cell(repeat,repeat);
inputs.A = 385;
%Loop over repeats of P
for i = 1:repeat
    Prep = PnoiseMBL(:,i+3);
    inputs.P = reshape(kron(Prep/24, ones(1,24))', ndays*24,1);
    %Loop over repeats of E
    for k = 1:repeat
        Qsim{i,k} = nan(ndays,repeat+3);
        Qsim{i,k}(:,1:3) = time;
        %20 repeats of Qsim
        for j = 4:repeat+3
            inputs.E = reshape(kron(Esim2{i,k}(:,j)/24, ones(1,24))', ndays*24,1);
            Qsim{i,k}(:,j) = PDMPieter(inputs, PDM_par); %Each repeat as a column in the cell of the (i,k)th repeat
        end
    end
end
disp('Q case 3 simulated.')

%% Save remaining output

save([loc_save names{1} '_' names{3} '_c3.mat'], 'Esim2', '-v7.3');
save([loc_save names{1} '_' names{5} '_c3.mat'], 'Qsim', '-v7.3');
save([loc_save names{1} '_' names{8} '_c3.mat'], 'time');
