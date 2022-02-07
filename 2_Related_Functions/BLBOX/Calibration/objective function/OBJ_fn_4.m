function [fitness, stat, statobs]=OBJ_fn_4(par,varargin)
% This objective function corresponds to a simplification of the objective function suggested
% by the theory of Hansen (1982), using the empirical variance matrix (diagonal matrix) as
% weighting matrix (see Kaczmarska (2011))
%  
% Used properties are specified by getRules.m 
%
% Load matrices created by OBJ_fn_init.m
model=varargin{1};
r=varargin{2};
o=varargin{3};
covgmm=varargin{5};
month=varargin{6};
vargmm=diag(diag(1./covgmm(:,:,month)));
% Load aggregation levels
agg=[1/6 1/2 1 6 12 24 48 72];
lagg=length(agg);


stat=[];
statobs=[];
for j=1:lagg
    if r(1,j)~=0 
        stat=[stat AVG(par,agg(j),model)];
        statobs=[statobs o(1,j,month)];
    end
    
    if r(2,j)~=0 
        stat=[stat VAR(par,agg(j),model)];
        statobs=[statobs o(2,j,month)];
    end
    
    if r(3,j)~=0
        stat=[stat COVAR(par,agg(j),1,model)];
        statobs=[statobs o(3,j,month)];
    end
    
    
    if r(4,j)~=0
        stat=[stat ZDP(par,agg(j),model)];
        statobs=[statobs o(4,j,month)];
    end


    if r(5,j)~=0
        stat=[stat CMOM3(par,agg(j),model)];
        statobs=[statobs o(5,j,month)];
    end  
end

if r(1,lagg+1)~=0 
    stat=[stat AVG(par,24*nrdiam(month),model)];     
    statobs=[statobs o(1,j,month)];
end

if r(2,lagg+1)~=0 
    stat=[stat VAR(par,24*nrdiam(month),model)];
    statobs=[statobs o(2,j,month)];
end

[rp cp]=size(par);
if cp==1
    par=par';
    [rp,cp]=size(par);
end

fitness=diag((stat-repmat(statobs,rp,1))*vargmm*(stat-repmat(statobs,rp,1))');

end
        
    