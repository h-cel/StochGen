function c=COVAR(par,h,lag,model)
%% AVG This function selects the analytical moment corresponding to the given model
%
%   This function selects the analytical moment for the lag-k covariance for the calibration of
%   Bartlett-Lewis models. 
%   
%   Inputs:
%       par: parameter set
%       h: aggregation level
%       lag: lag used for the calculation
%       model: model used
%   Output:
%       c: value for the covariance
%
%   Last updated by J. Van de Velde on 17/03/'21

%% Selection
switch lower(model)
    case 'obl'
        c=OBL_covar(par,h,lag);
    case 'mbl'
        c=MBL_covar(par,h,lag);
    case 'rbl2'
        c = RBL2_covar(par, h, lag);
    case 'mblg'
        c=MBLG_covar(par,h,lag);
    case 'tbl'
        c=TBL_covar(par,h,lag);
    case 'tblg'
        c=TBLG_covar(par,h,lag);
    case 'tblp'
        c=TBLP_covar(par,h,lag);
    case 'mblp'
        c=MBLP_covar(par,h,lag);
end
end
