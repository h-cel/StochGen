function [fitness stat]=OBJ_fn(par,varargin)
% This objective function corresponds to the objective function used in
% Verhoest et al. (1997)
% Used properties are specified by getRules.m 
%
% Load matrices created by OBJ_fn_init.m
model=varargin{1};
r=varargin{2};
o=varargin{3};
w=varargin{4};
month=varargin{6};

% Load aggregation levels
agg=[1/6 1/2 1 6 12 24 48 72];
lagg=length(agg);

% h=r.*repmat(agg,5,1);
stat=[];

for j=1:lagg+1
    if r(1,j)~=0 && j~=9
        stat=[stat ((AVG(par,agg(j),model)./o(1,j,month)-1).^2.*w(1,j,month))];
    elseif r(1,j)~=0 && j==9
        stat=[stat ((AVG(par,24*nrdiam(month),model)./o(1,j,month)-1).^2.*w(1,j,month))];        
    end
end


for j=1:lagg+1
    if r(2,j)~=0 && j~=9
        stat=[stat ((VAR(par,agg(j),model)./o(2,j,month)-1).^2.*w(2,j,month))];
    elseif r(2,j)~=0 && j==9
        stat=[stat ((VAR(par,24*nrdiam(month),model)./o(2,j,month)-1).^2.*w(2,j,month))];
    end
end

for j=1:lagg
    if r(3,j)~=0
    stat=[stat ((COVAR(par,agg(j),1,model)./o(3,j,month)-1).^2.*w(3,j,month))];
    end
end

for j=1:lagg
    if r(4,j)~=0
    stat=[stat ((ZDP(par,agg(j),model)./o(4,j,month)-1).^2.*w(4,j,month))];
    end
end

for j=1:lagg
    if r(5,j)~=0
    stat=[stat ((CMOM3(par,agg(j),model)./o(5,j,month)-1).^2.*w(5,j,month))];
    end
end


fitness=sum(stat,2);
        
    