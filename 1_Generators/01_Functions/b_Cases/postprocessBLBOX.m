function [Psimtmp,PnoiseMBL] = postprocessBLBOX(time, repeat, Psim, agg)
%POSTPROCESSBLBOX Postprocesses BLBOX output
%   This function postprocesses BLBOX output for use in the next steps of
%   the stochastic weather generator chain.

%   Last update by J. Van de Velde on 16/07/'21: made function by
%   refactoring

%% Aggregate to daily level

ndays= length(time);
Psimtmp = nan(ndays,repeat + 3);
for j=1:repeat % Aggregating to 24h for subsequent steps
    switch agg
        case '24'
            Psimtmp = Psim;
        case '1'
            row = size(Psim,1);
            step = 0:24:row;
            for i = 1:(row/24)
                Psimtmp(i,j+3) = sum(Psim((step(i)+1):step(i+1),j+4), 'omitnan');
            end
        case '1/6'
            row = size(Psim,1);
            step = 0:24:row;
            for i = 1:(row/24)
                Psimtmp(i,j+3) = sum(sum(Psim{j}((step(i)+1):step(i+1),5:10)), 'omitnan');
            end
    end
end

%% Add noise

vars2    = cell(1,size(Psimtmp,2)-3);
vars2(:) = {'P'};
PnoiseMBL = addNoise(Psimtmp,vars2,0.001); %Adds noise to the simulated values of P
PnoiseMBL(:, 1:3) = time;

end

