function [CalStr output] = momFit(MODEL,OBSERVATIONS,OBJFN,ROUTINE,MONTH,REPEAT,RULES,BOUND,SAVE,LOAD,zdpflag,varargin)
% This function is used to fit a Bartlett-Lewis models (MODEL) to a
% set of OBSERVATIONS, using the generalized method of moments. Different
% objective functions are possible (OBJFN), which differ in configuration
% and choice of weights. The fitting is conducted by an optimization
% algorithm of choice (ROUTINE), for a given set of months (MONTH), and is
% repeated REPEAT times. The resulting structure, containing the results is
% saved to SAVE. LOAD, loads an existing (unfinished) structure.
%
%
%
% MODEL     = model of choice to be fitted, string
%             - 'obl'
%             - 'mbl'
%             - 'rbl2'
%             - 'mblg'
%             - 'mblp'
%             - 'tbl'  (see Verhoest et al. (1997) and Onof (2011)) 
%             - 'tblg
%          or - 'tblp' 
%
% OBSERVATIONS = Uccle style matrix of observations
%
% OBJFN = name of the objective function, preceded by function handle or as
% string
%
% ROUTINE = optimization routine, string
%           - 'pso' (Particle Swarm Optimization)
%           - 'simpsa' (Simplex-Simulated Annealing)
%           - 'neldermead' (Simplex method by Nelder & Mead)
%           - 'sce' (Shuffled Complex Evolution)
%
% MONTH = months to be calibrated, integer vector (1-12)
%
% REPEAT = integer, number of repetitions
%
% SAVE = string, location to which the resulting structure will be saved
%
% LOAD = string, location from where initial structure is to be loaded
%
% Willem Jan Vanhaute 25/08/2011
%
%   Last updated by J. Van de Velde on 17/03/'21

%% Fit

% check if OBSERVATIONS is string or matrix, load matrix if string.
% Function matload must be added to path
if ischar(OBSERVATIONS)
    OBSERVATIONS=matload(OBSERVATIONS);
end

if ~isempty(varargin)
    W=varargin{1};
else
    W=[];
end
% Define potential aggregation levels used in fitting
agg     = [1/6 1/2 1 6 12 24 48 72];
while agg(1)<1/length(OBSERVATIONS(1,8:end))
    agg(1)=[];
end

% change boundaries to a more easy to handle formation
LB      =BOUND(:,1);
UB      =BOUND(:,2);
global dim
dim = length(BOUND);

% initialize the objective function. OBJ_fn_init is independent of the used
% OBJFN
initout=OBJ_fn_init(MODEL,RULES,OBSERVATIONS,agg,zdpflag);

%Initialize the matrices for data-storage
if ~isempty(LOAD)
    loadstruct=matload(LOAD);
    x           = loadstruct.x;
    z           = loadstruct.z;
    nFUN_EVAL   = loadstruct.FEs;
    EXITFLAG    = loadstruct.EXITFLAG;
    OPTIONS     = loadstruct.OPTIONS;
    TIME        = loadstruct.TIME; 
else
    x           = zeros(REPEAT,length(BOUND),length(MONTH));
    z           = zeros(REPEAT,length(MONTH));
    nFUN_EVAL   = zeros(REPEAT,length(MONTH));
    EXITFLAG    = zeros(REPEAT,length(MONTH));
    OPTIONS     = getOptions(MODEL,ROUTINE);
    TIME        = zeros(REPEAT,length(MONTH)); 
end

% Start of fitting loop
for i=1:length(MONTH)    
    for j=1:REPEAT
        mo=MONTH(i);        
        if max(x(j,:,i))==0 && z(j,i)==0 && nFUN_EVAL(j,i)==0 
            disp(['model: ' MODEL ', routine: ' ROUTINE ', month: ' num2str(mo) ', repetition: ' num2str(j)]);        
            tic
            switch lower(ROUTINE)
                case 'pso' 
                    rout='pso';
                    [x(j,:,i),z(j,i),nFUN_EVAL(j,i),EXITFLAG(j,i)]=PSO(OBJFN,BOUND,[],initout{:},mo,W);
                case 'simpsa'
                    rout='simpsa';
                    initpoint=rand(1,length(BOUND)).*repmat((UB-LB)',1,1)+repmat(LB',1,1); 
                    [x(j,:,i),z(j,i),EXITFLAG(j,i),out]=SIMPSA(OBJFN,initpoint,LB,UB,OPTIONS,initout{:},mo,W);
                    nFUN_EVAL(j,i)=out.nFUNEVALS;
                case 'neldermead'
                    rout='neldermead';
                    initpoint=(rand(1,length(BOUND)).*repmat((UB-LB)',1,1)+repmat(LB',1,1))'; 
                    [x(j,:,i),z(j,i),EXITFLAG(j,i),out]=fminsearchbnd(OBJFN,initpoint,LB,UB,OPTIONS,initout{:},mo,W);
                    nFUN_EVAL(j,i)=out.funcCount;
                    for r=1:29 %% repeat the simplex algorithm 30 times to increase accuracy
                        initpoint=(rand(1,length(BOUND)).*repmat((UB-LB)',1,1)+repmat(LB',1,1))'; 
                        [xr,zr,exi,out2]=fminsearchbnd(OBJFN,initpoint,LB,UB,OPTIONS,initout{:},mo,W);                
                        if zr<z(j,i)
                            z(j,i)=zr;
                            x(j,:,i)=xr;
                            EXITFLAG(j,i)=exi;
                            out=out2;
                            nFUN_EVAL(j,i)=out.funcCount;
                        end
                    end
                case 'sce'
                    rout='sce';
                    initpoint=(rand(1,length(BOUND)).*repmat((UB-LB)',1,1)+repmat(LB',1,1))'; 
                    [x(j,:,i),z(j,i),EXITFLAG(j,i),out]=SCE(OBJFN,initpoint,LB,UB,OPTIONS,initout{:},mo,W);
                    nFUN_EVAL(j,i)=out.nFUN_EVALS;
                    output=out;
            end
            TIME(j,i) = toc;
        end
    end
end
CalStr=struct('x',x,'z',z,'ROUTINE',rout,'OBJFN',OBJFN,'FEs',nFUN_EVAL,'RULES',RULES,'BOUND',BOUND,'OPTIONS',OPTIONS,'EXITFLAG',EXITFLAG,'TIME',TIME);
if ~isempty(SAVE)
save(SAVE,'CalStr');
end
% The output of momFit is the structure CalStr. CalStr consist of the
% caliberation results as well as information about the calibration in
% itself. 
% x         parameters
% z         corresponding objective function values 
% ROUTINE   used optimization method
% OBJFN     used objective function
% FEs       number of function evaluations for any given monthly repetition
% RULES     rule matrix (see getRules)
% BOUND     boundaries to the parameter space (see getBoundaries)
% OPTIONS   options of the optimization method (see getOptions)
% EXITFLAG  exitflags for the optimization (see optimization method script
% in question)
% TIME      time needed to complete calibration




