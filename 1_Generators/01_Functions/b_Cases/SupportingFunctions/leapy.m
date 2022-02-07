function [I] = leapy(year)
%% Leapy determines if given year is a leap year
%
%   Input:
%       year: double, year
%
%   Output: 
%       I: Boolean, indicating whether or not 'year' is a leap year
%
%   Last updated by J. Van de Velde on 24/09/'20

%% Calculcation

        I = rem(year,4)==0 & (rem(year,100)~=0|rem(year,400)==0);% leap year
end