% structminer will enable the user to mine the results of repeated
% calibrations of a certain Bartlett-Lewis model for its lowest value. The
% input is the structure resulting from calibration.
% 
% [param,z,details]=structminer(calstruct) will return the parameter set
% with the lowest function value z, and details about the calibration
% proces such as used routine, boundaries, number of function evaluations,
% calibration duration, etc.

function [param,zmin,details]=structminer(calstruct)
% unpack the structure 
x = calstruct.x;
z = calstruct.z;
rou = calstruct.ROUTINE;
% of  = calstruct.OBJFN;
fe  = calstruct.FEs;
bnd = calstruct.BOUND;
opt = calstruct.OPTIONS;
exi = calstruct.EXITFLAG;
tim = calstruct.TIME;
sc  = size(calstruct.z);
% tot = calstruct.TOTALTIME;

% find the minimum z value
z(z==0)=9999999;
param=zeros(sc(2),length(x(1,:,1)));
fe2=zeros(sc(2),1);
tim2=zeros(sc(2),1);
[zmin pos]=min(z);
for i = 1:sc(2)
    param(i,:)=x(pos(i),:,i);
    fe2(i)=fe(pos(i),i);
    tim2(i)=tim(pos(i),i);
end
zmin=zmin';
% details=struct('param',param,'z',zmin,'ROUTINE',rou,'OBJFN',of,'FEs',fe2,'BOUND',bnd, ...
%    'OPTIONS',opt,'EXITFLAG',exi,'TIME',tim2);
% 

    
    