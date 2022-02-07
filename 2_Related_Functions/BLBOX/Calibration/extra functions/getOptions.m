function options=getOptions(model,routine)

%   This function returns a structure containing the options for the chosen
%   optimization routine. Options will differ in function of the model
%   under calibration.
popgrootte=30; % size of the population (for PSO)
tolx=10^(-6);  % stopping criterium: if the euclidian distance between subsequent parameter estimations is smaller than tolx, the algorithm will stop.
tolfun=10^(-6);% stopping criterium: if the difference between subsequent objective function values of parameter estimations is smaller than tolfun, the algorithm will stop. 
MaxIter=5000;  % stopping criterium: maximum number of iterations
global dim
MaxFunEvals=round(MaxIter*popgrootte/dim/(dim+1)); % stopping criterium: maximum number of function evaluations
change=10^10;  % optional criterium for PSO: if global best does not change for 'change' number of steps, the algorithm will stop 
display='yes'; % change to yes if detailed output of calibration is required
%% 
switch lower(model) % Switch to model (results as in Vanhaute et al. 2012)
    case 'obl'
        switch lower(routine)
            case 'pso'
                
                c1=1.9;
                c2=1.8;
                w=0.3;
            case 'simpsa'
                coolrate=0.8;
            case 'sce'
                % parameters obtained by exhaustive parameter search
                % (August, 2011) 
                nC = 5;
                nI = 20;
        end
    case 'mbl'
       switch lower(routine)
            case 'pso' 
                % parameters obtained by exhaustive parameter search
                % (August, 2011) 
                c1=1.5;
                c2=2.5;
                w=0.4;
            case 'simpsa'
                % parameters obtained by exhaustive parameter search
                % (September, 2011) 
                coolrate=0.5;
                tempend=1;
           case 'sce'
                % parameters obtained by exhaustive parameter search
                % (August, 2011) 
                nC = 5;
                nI = 20;
        end
    case 'mblg'
       switch lower(routine)
            case 'pso' 
                
                c1=1.9;
                c2=1.4;
                w=0.4;
            case 'simpsa'
                coolrate=0.6;
            case 'sce'
                % parameters obtained by exhaustive parameter search
                % (August, 2011) 
                nC = 5;
                nI = 20;
        end
    case 'tbl'
        switch lower(routine)
            case 'pso'
                
                c1=1.9;
                c2=1.4;
                w=0.4;
            case 'simpsa'
                coolrate=0.3;
            case 'sce'
               
                nC = 10;
                nI = 20;
        end
    otherwise % For models not included in Vanhaute et al. 2012
       switch lower(routine)
            case 'pso'
                
                c1=1.9;
                c2=1.4;
                w=0.4;
            case 'simpsa'
                coolrate=0.3;
            case 'sce'
               
                nC = 10;
                nI = 20;
       end
        
end
%% ordering options in structure
switch lower(routine)
    case 'pso'
        if exist('c1')
            options=struct('popgrootte',popgrootte,'c1',c1,'c2',c2,'w',w,'MaxIter',MaxIter,'vmax',1,'tolpso',tolx,'change',change,'tolfun',tolfun,'display',display);
        else
            options=[];
        end
    case 'simpsa'
        tempend=1;
        if strcmp(display,'yes')
            options=SIMPSASET('COOL_RATE',coolrate,'TEMP_END',tempend, ... 
                     'TOLFUN',tolfun,'TOLX',tolx,'MAX_FUN_EVALS',MaxFunEvals,'DISPLAY','iter');
        else
            options=SIMPSASET('COOL_RATE',coolrate,'TEMP_END',tempend, ... 
                     'TOLFUN',tolfun,'TOLX',tolx,'MAX_FUN_EVALS',MaxFunEvals);
        end
    case 'neldermead'
        if strcmp(display,'yes')
            MaxFunEvals=30*MaxIter;
            options=OPTIMSET('TolFun',tolfun,'TolX',tolx,'MaxFunEvals',MaxFunEvals,'Display','iter');
        else
            options=OPTIMSET('TolFun',tolfun,'TolX',tolx,'MaxFunEvals',MaxFunEvals);
        end
    case 'sce'
        if strcmp(display,'yes')
            options=SCESET('nCOMPLEXES',nC,'nITER_INNER_LOOP',nI,'MAX_ITER', ...
                MaxIter,'MAX_FUN_EVALS',MaxFunEvals,'TOLX',tolx,'TOLFUN',tolfun,'DISPLAY','iter');
                     
        else
            options=SCESET('nCOMPLEXES',nC,'nITER_INNER_LOOP',nI,'MAX_ITER', ...
                MaxIter,'MAX_FUN_EVALS',MaxFunEvals,'TOLX',tolx,'TOLFUN',tolfun);
        end 
end
end

