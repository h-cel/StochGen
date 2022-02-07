function out=Ucc2row(data,agg)
n=length(data);
data=reshape(data(:,8:13)',6*agg,n/agg)';
out=sum(data,2);
end