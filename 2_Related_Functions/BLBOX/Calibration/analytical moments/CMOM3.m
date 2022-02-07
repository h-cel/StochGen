function m3=CMOM3(par,h,model)
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
%       m3: value for the third central moment
%
%   Last updated by J. Van de Velde on 17/03/'21

%% Selection
switch lower(model)
    case 'mbl'
        m3=blr_cmom3e(par,h);
    case 'rbl2'
        m3 = RBL2_cmom3e(par,h);
    case 'tbl'
        m3=tbl_cmom3e(par,h);
    case 'obl'
        m3=blo_cmom3e(par,h);
    case 'mblg'
        m3=blrg_cmom3e(par,h);
    case 'tblg'
        m3=tblg_cmom3e(par,h);
    case 'tblp'
        m3=tblp_cmom3e(par,h);
    case 'mblp'
        m3=blrp_cmom3e(par,h);

end
end
