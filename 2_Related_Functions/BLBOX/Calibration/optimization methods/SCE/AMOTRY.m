% AMOTRY FUNCTION
% ---------------

function [YTRY,PTRY] = AMOTRY(FUN,P,FAC,LB,UB,varargin)
% Extrapolates by a factor FAC through the face of the simplex across from 
% the high point, tries it, and replaces the high point if the new point is 
% better.

global NDIM

% calculate coordinates of new vertex
PSUM = sum(P(1:NDIM,:))/NDIM;
PTRY = PSUM*(1-FAC)+P(end,:)*FAC;

% evaluate the function at the trial point.
YTRY = CALCULATE_COST(FUN,PTRY,LB,UB,varargin{:});

return
