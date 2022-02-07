% COST FUNCTION EVALUATION
% ------------------------

function YTRY = CALCULATE_COST(FUN,PTRY,LB,UB,varargin)

global NDIM nFUN_EVALS

% add one to number of function evaluations
nFUN_EVALS = nFUN_EVALS + 1;

for i = 1:NDIM,
    % check lower bounds
    if PTRY(i) < LB(i),
        YTRY = 1e12+(LB(i)-PTRY(i))*1e6;
        return
    end
    % check upper bounds
    if PTRY(i) > UB(i),
        YTRY = 1e12+(PTRY(i)-UB(i))*1e6;
        return
    end
end

% calculate cost associated with PTRY
YTRY = feval(FUN,PTRY,varargin{:});

return