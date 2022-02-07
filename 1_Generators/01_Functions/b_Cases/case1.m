%% Case 1    P&T measurements -> E simulation -> Q simulation

%% Setup

loc_save = 'E:\Users\jpvdveld\Onderzoek\Data\3_cases\';
           
%% Simulations of E using original P and original T

%Initialization
ndays = length(time(:,1));
Esim = nan(ndays,repeat+3); %Columns: time + each repeat
Esim(:,1:3) = time;
Eini = E(1); %E: created in launchCases
%Loop over repeats
for i = 1:repeat 
    Esim(:,i+3) =  generateE(T,Pnoise,RE,thetasE,familyTPEpE,time,Eini,dataMonthE, copMatrixE);
    fprintf('E %d of %d simulated.\n',i,repeat);
    for n = 1:length(vars)
        if strcmp(vars{n}, 'E') == 1
            Esim(:,i+3) = Esim(:,i+3).*dataStand(:,n+6) + dataStand(:,n+3); %Remove standardization
        end
    end
end
disp('E case 1 simulated.')

%% Simulation of Q using original P and stochastic E 

%Initialization
Qsim = nan(ndays,repeat+3);
Qsim(:,1:3) = time;
inputs.P = reshape(kron(P/24, ones(1,24))', ndays*24,1);
inputs.A = 385;
%Simulation of Q
for i = 4:repeat+3 %For each repeat
    inputs.E = reshape(kron(Esim(:,i)/24, ones(1,24))', ndays*24,1);
    Qsim(:,i) = PDMPieter(inputs, PDM_par); %Each repeat as a column, like E.
end
disp('Q case 1 simulated.')

%% Save output

save([loc_save names{1} '_' names{3} '_c1.mat'], 'Esim');
save([loc_save names{1} '_' names{5} '_c1.mat'], 'Qsim');
save([loc_save names{1} '_time_c1.mat'], 'time');
