function E = RBL2_avg(theta,h)
%% RBL2_avg This function calculates the average precipitation for the RBL2 model.
%
%   This function calculates the average precipitation for the RBL2 model
%   (Kaczmarska et al., 2014) based upon the equations derived in Onof and
%   Wang (2020).
%   
%   Inputs:
%       theta: parameters for this model
%       h: aggregation level
%   Output:
%       E: mean precipitation for every month at aggregation level h
%
% Last updated by J. Van de Velde on 16/03/'21

%% Loading parameters
lambda=theta(:,1);
%nu=theta(:,2);
%alpha=theta(:,3);
iota=theta(:,4);
phi=theta(:,5);
kappa=theta(:,6);
%omega = theta(:,7);

%% Calculations

muc=1+kappa./phi;

E = lambda.*muc.*h.*iota;

end