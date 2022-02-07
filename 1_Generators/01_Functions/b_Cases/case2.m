%% Case 2    P measurements -> T simulation -> E simulation -> Q simulation

%Locations
loc_save = 'E:\Users\jpvdveld\Onderzoek\Data\3_cases\';

%% Simulations of T using original P
%Initialization
ndays = length(time(:, 1));
Tsim = nan(ndays,repeat+3);
Tsim(:,1:3) = time;
Tini = T(1); %T: created in launchCases
for i = 1:repeat %Number of repeats
    Tsim(:,i+3) = generateT(Pnoise,RT,thetasT,familyTpPT,time,Tini,dataMonthT, copMatrixT);
    fprintf('T %d of %d simulated.\n', i,repeat);
end
disp('T case 2 simulated.')

%% Simulations of E using original P and stochastic T
%Initialization
Esim = cell(1,repeat);
%Loop over repeats of T
for i = 1:repeat 
    Trep = Tsim(:,i+3);
    Esim{i}(:,1:3) = time;
    Eini = E(1); %E: created in launchCases; this wasn't implemented, but needed. Not sure if this is the right way.
    %Loop over repeats of E (-> 20*20 = 400 repeats created!)
    for j = 1:repeat
        Esim{i}(:,j+3) =  generateE(Trep,Pnoise,RE,thetasE,familyTPEpE,time,Eini,dataMonthE, copMatrixE);
        fprintf('E %d/%d of repetition %d/%d simulated.\n', j, repeat,i,repeat);
        for n = 1:length(vars)
            if strcmp(vars{n}, 'E') == 1
                Esim{i}(:,j+3) = Esim{i}(:,j+3).*dataStand(:,n+6) + dataStand(:,n+3); %Remove standardization
            end
        end
    end
end
disp('E case 2 simulated.')

%% Simulations of Q using original P and destandardized stochastic E

%Initialization
Qsim = cell(1,repeat);
inputs.P = reshape(kron(P/24, ones(1,24))', ndays*24,1);
inputs.A = 385;
%Loop over repeats
for i = 1:repeat
    Qsim{i} = nan(ndays,repeat+3);
    Qsim{i}(:,1:3) = time;
    %Loop over 20 E repeats
    for j = 4:repeat+3
        inputs.E = reshape(kron(Esim{i}(:,j)/24, ones(1,24))', ndays*24,1);
        Qsim{i}(:,j) = PDMPieter(inputs, PDM_par); %Each repeat as a column in the cell of the ith repeat
    end
end
disp('Q case 2 simulated.')

%% Destandardize T

for i = 1:repeat
    for n = 1:length(vars)
        if strcmp(vars{n}, 'T') == 1 
            Tsim(:,i+3) = Tsim(:,i+3).*dataStand(:,n+6) + dataStand(:,n+3);
        end
    end
end

%% Save output

save([loc_save names{1} '_' names{3} '_c2.mat'], 'Esim');
save([loc_save names{1} '_' names{4} '_c2.mat'], 'Tsim');
save([loc_save names{1} '_' names{5} '_c2.mat'], 'Qsim');
save([loc_save names{1} '_time_c2.mat'], 'time');


