function PH=RBL2_zdp(theta,h)
%% RBL2_zdp This function calculates the covariance for the RBL2 model.
%
%   This function calculates the zero depth probability at aggregation level h for the RBL2 model
%   (Kaczmarska et al., 2014), used as derived for the RBL1 model in
%   Rodriguez-Iturbe et al. (1988) (as was done in Kim and Onof (2020)
%   Although note that the exact origin of this function is unclear.
%   
%   Inputs:
%       theta: parameters for this model
%       h: aggregation level
%   Output:
%       PH: proportion dry for each month at aggregation level h
%
% Last updated by J. Van de Velde on 16/03/'21

%% Parameters and setup
[r,c]=size(theta);
lambda=theta(:,1);
nu=theta(:,2);
alpha=theta(:,3);
iota=theta(:,4);
phi=theta(:,5);
kappa=theta(:,6);
omega = theta(:,7);

eta=gamrnd(alpha,1/nu);
[mu_x, var_x] = gamstat(omega, iota*eta/omega);
sigmax = sqrt(var_x);
eps=0.001;   

%% Calculation

mux3 = 2.* sigmax.^3 + 3.*mu_x.*sigmax.^2 + mu_x.^3;
muc = 1+kappa./phi;
mux2 = sigmax.^2+mu_x.^2;
mut=( phi^(-1) + phi*inter2(kappa, phi) ) * Exp_fn_e(1, 0, eps, nu, alpha);

PH = exp( -lambda*(h+mut) + lambda*exp(-kappa)*(kappa+phi)^(-1)* inter(kappa, phi) * ...
    (phi*Exp_fn_e(1,0,eps,nu,alpha) + kappa*Exp_fn_e(1,(kappa+phi)*h,eps,nu,alpha)));
end

function interpdry=inter(kappa,phi)
 M = 5;
 j = 0:M;
 interpdry =sum(kappa.^j./(gamma(j+1).*(1+j+phi).*(j+phi)));
end

function mut=inter2(kappa,phi)
 M = 30;
 j =1:M;
 mut = 1/phi + sum((-kappa).^(j-1).*(kappa-j.^2-j).* ...
                                    exp(log(gamma(phi))-log(gamma(j+1+phi)))./ ...
					                     (j.*(j+1)));
end