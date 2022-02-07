function m=AVG(par,h,model)
%% AVG This function selects the analytical moment corresponding to the given model
%
%   This function selects the analytical moment for the calibration of
%   Bartlett-Lewis models. 
%   
%   Inputs:
%       par: parameter set
%       h: aggregation level
%       model: model used
%   Output:
%       m: value for the average
%
%   Last updated by J. Van de Velde on 17/03/'21

%% Selection
switch lower(model)
    case 'obl'
        m=OBL_avg(par,h);
    case 'mbl'
        m=MBL_avg(par,h);
    case 'rbl2'
        m =RBL2_avg(par,h);
    case 'mblg'
        m=MBLG_avg(par,h);
    case 'tbl'
        m=TBL_avg(par,h);
    case 'tblg'
        m=TBLG_avg(par,h);
    case 'tblp'
        m=TBLP_avg(par,h);
    case 'mblp'
        m=MBLP_avg(par,h);

end
end
