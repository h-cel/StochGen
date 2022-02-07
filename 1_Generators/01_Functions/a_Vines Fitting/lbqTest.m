function [h,pValue,stat,cValue] = lbqTest(res,varargin)
%LBQTEST Ljung-Box Q-test for residual autocorrelation
%
% Syntax:
%
%   [h,pValue,stat,cValue] = lbqtest(res)
%   [h,pValue,stat,cValue] = lbqtest(res,param1,val1,param2,val2,...)
%
% Description:
%
%   The "portmanteau" test of Ljung and Box assesses the null hypothesis
%   that a series of residuals exhibits no autocorrelation for a fixed
%   number of lags L, against the alternative that some autocorrelation
%   coefficient rho(k), k = 1, ..., L, is nonzero. The test statistic is
%
%                  L
%   	Q = T(T+2)Sum(rho(k)^2/(T-k)),
%                 k=1
%
%   where T is the sample size, L is the number of autocorrelation lags,
%   and rho(k) is the sample autocorrelation at lag k. Under the null, the
%   asymptotic distribution of Q is chi-square with L degrees of freedom.
%
% Input Arguments:
%
%   res - Vector of residuals for which the test statistic is computed. The
%   	last element corresponds to the most recent observation. Typically,
%       res contains (standardized) residuals obtained by fitting a model
%       to an observed time series.
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME        VALUE
%
%   'lags'      Scalar or vector of positive integers indicating the number
%               of lags L used to compute the test statistic. Each element
%               must be less than the length of res. The default value is
%               min[20,T-1].
%
%   'alpha'     Scalar or vector of nominal significance levels for the
%               tests. Elements must be greater than zero and less than
%               one. The default value is 0.05.
%
%   'dof'       Scalar or vector of degree-of-freedom parameters for the
%               asymptotic chi-square distributions of the test statistics.
%               Elements must be positive integers less than the
%               corresponding element of lags. The default value is the
%               value of lags.
%
%   Scalar parameter values are expanded to the length of any vector value
%   (the number of tests). Vector values must have equal length. If any
%   value is a row vector, all outputs are row vectors.
%
% Output Arguments:
%
%   h - Vector of Boolean decisions for the tests, with length equal to the
%   	number of tests. Values of h equal to 1 indicate rejection of the
%       null of no autocorrelation in favor of the alternative. Values of h
%       equal to 0 indicate a failure to reject the null.
%
%   pValue - Vector of p-values of the test statistics, with length equal
%       to the number of tests.
%
%   stat - Vector of test statistics, with length equal to the number of
%       tests.
%
%   cValue - Vector of critical values for the tests, determined by alpha,
%       with length equal to the number of tests.
%
% Notes:
%
%   o The input lags affects the power of the test. If L is too small, the
%     test will not detect high-order autocorrelations; if it is too large,
%     the test will lose power when significant correlation at one lag is
%     washed out by insignificant correlations at other lags. The default
%     value of min[20,T-1] is suggested by Box, Jenkins, and Reinsel [1].
%     Tsay [4] cites simulation evidence that a value approximating log(T)
%     provides better power performance.
%
%   o When res is obtained by fitting a model to data, the degrees of
%     freedom are reduced by the number of estimated coefficients,
%     excluding constants. For example, if res is obtained by fitting an
%     ARMA(p,q) model, dof should be L-p-q.
%
%   o LBQTEST does not test directly for serial dependencies other than
%     autocorrelation, but it can be used to identify conditional
%     heteroscedasticity (ARCH effects) by testing squared residuals. See,
%     e.g., McLeod and Li [3]. Engle's test, implemented by ARCHTEST, tests
%     for ARCH effects directly.
%
% Example:
%
%   % Test exchange rates for autocorrelation, ARCH effects:
%
%   load Data_MarkPound
%   returns = price2ret(Data);
%   residuals = returns-mean(returns);
%   h1 = lbqtest(residuals)
%   h2 = lbqtest(residuals.^2)
%
% References:
%
%   [1] Box, G.E.P., G.M. Jenkins, and G.C. Reinsel. Time Series Analysis:
%       Forecasting and Control. 3rd ed. Upper Saddle River, NJ:
%       Prentice-Hall, 1994.
% 
%   [2] Gourieroux, C. ARCH Models and Financial Applications. New York:
%       Springer-Verlag, 1997.
%
%   [3] McLeod, A.I. and W.K. Li. "Diagnostic Checking ARMA Time Series
%       Models Using Squared-Residual Autocorrelations." Journal of Time
%       Series Analysis. Vol. 4, 1983, pp. 269-273.
%
%   [4] Tsay,R.S. Analysis of Financial Time Series. Hoboken, NJ: John
%       Wiley & Sons, Inc., 2005.
%
% See also AUTOCORR, ARCHTEST.

% Copyright 1999-2010 The MathWorks, Inc.   

% Parse inputs and set defaults. This parse supports the deprecated syntax
% in which lags, alpha, and dof are input as the second, third, and fourth
% arguments, respectively. However, parameter/value pairs for lags, alpha,
% or dof override any value in the deprecated syntax.

T = length(res); % Effective sample size
defaultLags = min(20,T-1); % Recommended in [1]
defaultAlpha = 0.05;
defaultDof = defaultLags;

parseObj = inputParser;
parseObj.addRequired('res',@resCheck);
parseObj.addOptional('valueBasedLags',defaultLags,@lagsCheck);
parseObj.addOptional('valueBasedAlpha',defaultAlpha,@alphaCheck);
parseObj.addOptional('valueBasedDof',defaultDof,@dofCheck);
parseObj.addParameter('lags',defaultLags,@lagsCheck);
parseObj.addParameter('alpha',defaultAlpha,@alphaCheck);
parseObj.addParameter('dof',defaultDof,@dofCheck);

parseObj.parse(res,varargin{:});

res = parseObj.Results.res;
valueBasedLags = parseObj.Results.valueBasedLags;   % Deprecated
valueBasedAlpha = parseObj.Results.valueBasedAlpha; % Deprecated
valueBasedDof = parseObj.Results.valueBasedDof;     % Deprecated
lags = parseObj.Results.lags;
alpha = parseObj.Results.alpha;
dof = parseObj.Results.dof;

if ~isequal(valueBasedLags,defaultLags) && ~sum(strcmpi('lags',varargin))
    lags = valueBasedLags; % Accept lags in deprecated syntax
end

if ~isequal(valueBasedAlpha,defaultAlpha) && ~sum(strcmpi('alpha',varargin))
    alpha = valueBasedAlpha; % Accept alpha in deprecated syntax
end

if ~isequal(valueBasedDof,defaultDof) && ~sum(strcmpi('dof',varargin))
    dof = valueBasedDof; % Accept dof in deprecated syntax
end

% Check for non-default lags and unspecified dof:
if ~isequal(lags,defaultLags) && ... % Default dof differs from default lags
   ~sum(strcmpi('dof',varargin)) && ... % Dof unspecified as parameter/value pair
   (nargin < 4 || ischar(varargin{2}) || ischar(varargin{3})) % Dof unspecified as value-based argument

    dof = lags; % Reset dof default
    
end

% Check parameter values for commensurate lengths, expand scalars, and
% convert all variables to columns:

[rowOutput,lags,alpha,dof] = sizeCheck(lags,alpha,dof);
res = res(:);

if any(dof > lags)
       
	error(message('DofTooLarge'))
        
end

% Compute the sample ACF out to the largest lag:

maxLag = max(lags);
ACF = autocorr(res,maxLag); % Lags 0, 1, ..., maxLag
%ACF = ACF(2:end);           % Strip off ACF at lag 0

% Compute Q-statistics to the largest lag; keep only those requested:

idx = (T-(1:maxLag))';
stat = T*(T+2)*cumsum((ACF.^2)./idx);
stat = stat(lags);

% Compute p-values:

pValue = 1-chi2cdf(stat,dof);

% Compute critical values, if requested:

if nargout >= 4
    
   cValue = chi2inv(1-alpha,dof);
   
else
    
   cValue = [];
   
end

% Perform the test:

h = (alpha >= pValue);

% Display outputs as row vectors if lags is a row vector:

if rowOutput
    
   h = h';
   pValue = pValue';
   stat = stat';

   if ~isempty(cValue)
       
      cValue = cValue';
      
   end
   
end

%-------------------------------------------------------------------------
% Check input res
function OK = resCheck(res)

    if isempty(res)
        
        error(message('DataUnspecified'))
          
    elseif ~isnumeric(res)
        
        error(message('DataNonNumeric'))
          
    elseif ~isvector(res)
        
        error(message('DataNonVector'))
          
    else
        
        OK = true;
        
    end

end % resCheck

%-------------------------------------------------------------------------
% Check value of 'lags' parameter (or value-based lags)
function OK = lagsCheck(lags)
    
    if ~isnumeric(lags)
        
        error(message('LagsNonNumeric'))
          
    elseif ~isvector(lags)
        
        error(message('LagsNonVector'))
          
    elseif any(mod(lags,1) ~= 0) || any(lags <= 0)
        
        error(message('LagsOutOfRange'))
          
    elseif any(lags > (T-1))
       
        error(message('LagsTooLarge'))
          
    else
        
        OK = true;
        
    end

end % lagsCheck

%-------------------------------------------------------------------------
% Check value of 'alpha' parameter (or value-based alpha)
function OK = alphaCheck(alpha)
    
    if ~isnumeric(alpha)
        
        error(message('AlphaNonNumeric'))
          
    elseif ~isvector(alpha)
        
        error(message('AlphaNonVector'))
          
	elseif any(alpha <= 0) || any(alpha >= 1)
        
        error(message('AlphaOutOfRange'))
          
    else
        
        OK = true;
        
    end

end % alphaCheck

%-------------------------------------------------------------------------
% Check value of 'dof' parameter (or value-based dof)
function OK = dofCheck(dof)
    
    if ~isnumeric(dof)
        
        error(message('DofNonNumeric'))
          
    elseif ~isvector(dof)
        
        error(message('DofNonVector'))
          
    elseif any(mod(dof,1) ~= 0) || any(dof <= 0)
        
        error(message('DofOutOfRange'))
          
    else
        
        OK = true;
        
    end
 
end % dofCheck

%-------------------------------------------------------------------------
% Check parameter values for commensurate lengths, expand scalars, and
% convert all variables to columns
function [rowOutput,varargout] = sizeCheck(varargin)

% Initialize outputs:

numTests = 1;
rowOutput = false;

% Determine vector lengths, row output flag:

for i = 1:nargin
        
    ivar = varargin{i};
    iname = inputname(i);
    
    paramLength.(iname) = length(ivar);
    numTests = max(numTests,paramLength.(iname));
    
    if ~isscalar(ivar)
        rowOutput = rowOutput || (size(ivar,1) == 1);
    end    
    
end

% Check for commensurate vector lengths:

for i = 1:(nargin-1)
    iname = inputname(i);
    for j = (i+1):nargin
        jname = inputname(j);
        if (paramLength.(iname) > 1) && (paramLength.(jname) > 1) ...
            && (paramLength.(iname) ~= paramLength.(jname))
        
            error(message('ParameterSizeMismatch', iname, jname))
              
        end        
    end
end

% Expand scalars:

for i = 1:nargin
    
    ivar = varargin{i};
    if paramLength.(inputname(i)) == 1
        varargout{i} = ivar(ones(numTests,1)); %#ok
    else
        varargout{i} = ivar(:);  %#ok Column output 
    end
    
end

end % sizeCheck

end % LBQTEST