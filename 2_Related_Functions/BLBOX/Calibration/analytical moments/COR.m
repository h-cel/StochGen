function corr=COR(par,h,lag,model)
% Obsolete? Seems to be replaced by the COVAR.m function
switch lower(model)
    case 'obl'
        corr=OBL_covar(par,h,lag);
    case 'mbl'
        corr=MBL_covar(par,h,lag);
    case 'mblg'
        corr=MBLG_covar(par,h,lag);
    case 'tbl'
        corr=TBL_covar(par,h,lag);
    case 'tblg'
        corr=TBLG_covar(par,h,lag);
    case 'tblp'
        corr=TBLP_covar(par,h);
    case 'mblp'
        corr=MBLP_covar(par,h);
end
corr=corr/VAR(par,h,model);
end