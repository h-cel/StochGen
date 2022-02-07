function S = RBL2_cmom3e(theta,h)
%% RBL2_cmom3e This function calculates the third central moment (skewness) for the RBL2 model.
%
%   This function calculates the third central moment of precipitation for the RBL2 model
%   (Kaczmarska et al., 2014) based upon the equations derived in Onof and
%   Wang (2020).
%   
%   Inputs:
%       theta: parameters for this model
%       h: aggregation level
%   Output:
%       S: skewness of precipitation for every month at
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
mux3 = 2.* sigmax.^3 + 3.*mu_x.*sigmax.^2 + mu_x.^3;
f1 = mux2/(mu_x^2);
f2 = mux3/(mu_x^3);
eta_0 = 0.001; % Onof and Wang (2020), p 2796

% Calculation

P1 = 6*T_integral(1,h,eta_0,nu, alpha)*phi^2*(phi*kappa^2*(2*phi^2-7*phi^2-3*phi+2)+2*phi*f2*(phi^6-6*phi^4+9*phi^2-4)+kappa*f1*(4*phi^6-22*phi^4-phi^3+25*phi^2+4*phi-4));

P2 = 6*T_integral(0,h,eta_0,nu,alpha)*phi^3*h*(f2*(phi^6-6*phi^4+9*phi^2-4)+phi*kappa*f1*(phi^2-1)*(phi^2-4));

P3 = 6*T_integral(1,phi*h,eta_0, nu, alpha)*kappa*(f1*(-phi^5+phi^4+6*phi^3-4*phi^2-8*phi)+kappa*(phi^5-3*phi^4+2*phi^3+14*phi^2-8));

P4 =6*T_integral(0,phi*h,eta_0,nu,alpha)*h*kappa^2*(phi^3*(5-phi^2)-4*phi);

P5 = T_integral(1,0,eta_0, nu, alpha)*(-12*phi^3*f2*(phi^6-6*phi^4+9*phi^2-4)+kappa^2*(-9*phi^7+39*phi^5+18*phi^4-12*phi^3-84*phi^2+48)-3*phi*kappa*f1*(7*phi^7-39*phi^5-2*phi^4+46*phi^3+12*phi^2-8*phi-16));

P6 = T_integral(0,0,eta_0,nu,alpha)*((6*h*phi^3*f2+12*h*phi^2*kappa*f1+6*h*phi*kappa^2)*(phi^6-6*phi^4+9*phi^2-4));

P7 = 3*T_integral(1,2*h,eta_0, nu, alpha)*phi^4*(1-phi^2)*iota^3*(phi*kappa^2+kappa*f1*(phi^2-4));

P8 = 6*T_integral(1,(1+phi)*h,eta_0, nu, alpha)*kappa*phi^2*(phi-2)*(phi-1)*iota^3*(f1*(phi+2)-phi*kappa);

if alpha > 1
    S = (lambda*muc*iota^3*(P1+P2+P3+P4+P5+P6+P7+P8))/((1+2*phi+phi^2)*(phi^4-2*phi^3-3*phi^2+8*phi-4)*phi^3);
elseif alpha>-2 && alpha <= 1
    S = (lambda*muc*iota^3)/((1+2*phi+phi^2)*(phi^4-2*phi^3-3*phi^2+8*phi-4)*phi^3)*((nu^alpha*eta_0^(alpha+2)*h^3)/(gamma(alpha)*(alpha+2))*(2*kappa^2*(phi^7-3*phi^6+phi^5+3*phi^4-2*phi^3)+f2*(phi^9-6*phi^7+9*phi^5-4*phi^3)+3*kappa*f1*(phi^8-phi^7-5*phi^6+5*phi^5+4*phi^4-4*phi^3))+(P1+P2+P3+P4+P4+P5+P6+P7+P8));
end
end

