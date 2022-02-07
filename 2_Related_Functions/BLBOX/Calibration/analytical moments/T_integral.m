function T_value = T_integral(k, u, l, nu, alpha)
%T_INTEGRAL Calculates the function T as used in Onof and Wang (2020)
%   
%   This function calculates T as given in Onof and Wang (2020). This
%   function is used in the analytical expressions of the moments for the
%   RBL2 model.
%   
%   Inputs:
%       k, u, l: parameters for the T function
%       nu, alpha: parameters of the RBL2 model
%   Outputs:
%       T_value: value of the function
%
%   Last updated by J. Van de Velde on 12/03/'21

%% Calculation

T_value = (nu.^alpha)/(nu+u).^(alpha-k)*gammainc(alpha-k, l*(nu+u))/gamma(alpha);
end

