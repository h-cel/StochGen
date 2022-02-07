function dp=ZDP(par,h,model)
%% AVG This function selects the analytical moment corresponding to the given model
%
%   This function selects the analytical moment for the zero depth probability for the calibration of
%   Bartlett-Lewis models. 
%   
%   Inputs:
%       par: parameter set
%       h: aggregation level
%       model: model used
%   Output:
%       dp: value for the depth probability
%
%   Last updated by J. Van de Velde on 17/03/'21

%% Selection

switch lower(model)
    case 'obl'
        dp=OBL_zdp(par,h);
    case 'mbl'
        dp=MBLg_zdp(par,h,model);
    case 'rbl2'
        dp = RBL2_zdp(par, h);
    case 'mblg'
        dp=MBLg_zdp(par,h,model);
    case 'tbl'
        dp=TBL_zdp(par,h);
    case 'tblg'
        dp=TBLG_zdp(par,h);
    case 'tblp'
        dp=TBLP_zdp(par,h);
    case 'mblp'
        dp=MBLP_zdp(par,h);

end
end