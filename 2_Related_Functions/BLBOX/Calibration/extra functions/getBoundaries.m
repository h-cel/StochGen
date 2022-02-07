function BOUND = getBoundaries(model)
%% getBoundaries This function retrieves the boundaries for a specified Bartlett-Lewis model
% BOUND : Specify search range per model
% For the modified models (i.e. everything except OBL) it is theoretically
% most correct to set the lower boundary of parameter alpha to 4. If alpha
% is lower than 4, the theoretical skewness is infinite, if alpha is lower
% than 3, the variance is infinite. 
%
% However, this claim on alpha was disputed in Onof and Wang (2020). Hence,
% for the boundaries of RBL2, a lower alpha value is chosen
%
%   Last updated by J. Van de Velde on 17/03/'21

%% Getboundaries

switch upper(model)
    case 'OBL'
        BOUND=[0,      0,     0,   0,    0; ...
        .1,   10,     1,    10,    10]';
    case 'MBL'
        % :,1: lambda|:,2: kappa|:, 3: phi|:,4: mux|:,5: alpha|:,6: nu (=1/theta); %
        % zoals in artikel Willem-Jan (Vanhaute et al, 2012)
        % lambda ondergrens:0.01 voor klimaatmodellen 
        BOUND=[0.01,      0,     0,   0,    4,    0; ...
        0.1,   20,    1,  15,   20, 20]';
    case 'RBL2'
        % :,1: lambda; :,2: nu; :,3: alpha; :,4: iota, :5: phi; :,6:
        % kappa; ;,7: omega
        BOUND =[0.01, 0, -2, 0, 0, 0, 0;...
            0.1, 20, 20, 20, 1, 20, 20]';
    case 'MBLG'
        BOUND=[0.01,0.05;0,2;0,1;4,100;0,50;0,10;0,5];
    case 'TBL'
        BOUND=[0.01,.1;0,1;0,.1;0,5;0,20;0,10;0,2];
    case 'TBLG'
        BOUND=[0.01,.1;0.1,2;0,.1;0,100;0,50;0,3;0,10;0,5;];        
    case 'TBLP'
        BOUND=[0.01,0.5;0,5;0,1;2,15;0,2;0,5;0,1;0,1;];      
    case 'MBLP'
        BOUND=[0,0.1;0,20;0,1;0,10;0,10;0,1/3;0.1,1];  
    case 'MBLA'
        BOUND=[0.01,0.05;0000001,2;0000001,5;0,8;0000001,5;0.000001,2;0.000001,1;0.000001,5];
    case 'MBLGFIT'
        BOUND=[0,2;0,1;2,8;0,1;0,3;0,2];  
end

