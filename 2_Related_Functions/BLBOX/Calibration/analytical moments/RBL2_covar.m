function C = RBL2_covar(theta, h, k)
%% RBL2_covar This function calculates the covariance for the RBL2 model.
%
%   This function calculates the lag-k covariance in precipitation for the RBL2 model
%   (Kaczmarska et al., 2014) based upon the equations derived in Onof and
%   Wang (2020).
%   
%   Inputs:
%       theta: parameters for this model
%       h: aggregation level
%       k: lag used
%   Output:
%       C: covariance in precipitation at lag k for every month at
%       aggregation level h
%
% Last updated by J. Van de Velde on 16/03/'21

%% Loading parameters
lambda=theta(:,1);
nu=theta(:,2);
alpha=theta(:,3);
iota=theta(:,4);
phi=theta(:,5);
kappa=theta(:,6);
omega = theta(:,7);

%% Calculation
% Extra parameters etc.
eta=gamrnd(alpha,1/nu);
muc=1+kappa./phi;
[mu_x, var_x] = gamstat(omega, iota*eta/omega);
sigmax = sqrt(var_x);
mux2 = sigmax.^2+mu_x.^2;
f1 = mux2/(mu_x^2);
eta_0 = 0.001; % Onof and Wang (2020), p 2796

% Calculation

if alpha > 1
    C = lambda*muc*iota^2*((f1+(kappa*phi)/(phi^2-1))*(T_integral(1,(k-1)*h,0, nu, alpha) -2*T_integral(1,k*h,0, nu, alpha)+T_integral(1,(k+1)*h,0, nu, alpha)) ...
        -(kappa/(phi^2*(phi^2-1)))*(T_integral(1,phi*(k-1)*h,0, nu, alpha)-2*T_integral(1,phi*k*h,0, nu, alpha) +T_integral(1,phi*(k+1)*h,0, nu, alpha)));
elseif alpha>-1 && alpha <= 1
    C = lambda*muc*iota^2*((eta_0^(alpha+1)*h^2*nu^alpha)/(gamma(alpha)*(alpha+1))*(f1+kappa/(phi+1))+(f1+(kappa*phi)/(phi^2-1))*...
        (T_integral(1,(k-1)*h,eta_0, nu, alpha)-2*T_integral(1,k*h, eta_0,nu, alpha)+T_integral(1,(k+1)*h,eta_0, nu, alpha)) ...
        - (kappa/(phi^2*(phi^2-1)))*(T_integral(1,phi*(k-1)*h, eta_0, nu, alpha)-2*T_integral(1,phi*k*h,eta_0, nu, alpha)+T_integral(1,phi*(k+1)*h,eta_0,nu, alpha)));
else
    C= NaN;
end
end

