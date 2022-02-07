function c=VAR(par,h,model)
%% AVG This function selects the analytical moment corresponding to the given model
%
%   This function selects the analytical moment function for the variance for the calibration of
%   Bartlett-Lewis models. 
%   
%   Inputs:
%       par: parameter set
%       h: aggregation level
%       model: model used
%   Output:
%       c: value for the variance
%
%   Last updated by J. Van de Velde on 17/03/'21

%% Selection

switch lower(model)
    case 'obl'
        c=OBL_var(par,h);
    case 'mbl'
        c=MBL_var(par,h);
    case 'rbl2'
        c=RBL2_var(par,h);
    case 'mblg'
        c=MBLG_var(par,h);
    case 'tbl'
        c=TBL_var(par,h);
    case 'tblg'
        c=TBLG_var(par,h);
    case 'tblp'
        c=TBLP_var(par,h);
    case 'mblp'
        c=MBLP_var(par,h);

end
end
