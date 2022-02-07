function out=OBJ_fn_init(model,rules,data,agg,zdpflag)
% This function initializes the fitting procedure. Model specifies which
% model is being fitted, rules -> see getRules.m, data corresponds to
% OBSERVATIONS in momFit.m, agg is a vector with aggregation levels. 

[rr rc]=size(rules);
wei=StatCalcm3(data,agg,zdpflag);
sw=size(wei);
if rc~=sw(2)+1
    error('Rule matrix column size must be equal to length of agg.');
end
x=monthlymean(data);
x=reshape(x',2,1,12);
x=[x;zeros(rr-2,1,12)];
wei=[wei,x];
w=zeros(rr,rc,12);
denom=zeros(rr,rc);
for i=1:numel(rules)
    if rules(i)>0
       denom(i)=1;
    end
end
o=repmat(denom,[1 1 12]).*wei;
w=repmat(denom,[1 1 12]).*repmat(rules,[1,1,12]);
covgmm=CovGMOM(data,agg,rules,zdpflag);


out=cell(5,1);
out{1}=model;
out{2}=rules;
out{3}=o;
out{4}=w;
out{5}=covgmm;
% save('OBJ_fn_covgmm','covgmm');
