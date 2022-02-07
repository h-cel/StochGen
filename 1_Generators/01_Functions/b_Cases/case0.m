%% Case 0    P&E measurements -> Q simulation

%% Simulations of Q using original P and original E 
%Initialization
ndays = size(time,1);
Qsim = nan(ndays,4);
Qsim(:,1:3) = time;
%Input preparation
inputs.P = reshape(kron(P/24, ones(1,24))', ndays*24,1); %Hourly precipitation input data
inputs.A = 385;
inputs.E = reshape(kron(E/24, ones(1,24))', ndays*24,1); %Hourly evaporation input data
%Simulation
Qsim(:,4) = PDMPieter(inputs, PDM_par);
disp('Q case 0 simulated.')

%% Save output

save(['E:\Users\jpvdveld\Onderzoek\Data\3_cases\' names{1} '_' names{5} '_c0.mat'], 'Qsim');
save([loc_save names{1} '_time_c0.mat'], 'time');
