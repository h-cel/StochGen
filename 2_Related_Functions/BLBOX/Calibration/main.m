% This is a step-by-step example of how to use the Bartlett-Lewis Optimization toolbox
%%
% 1. Add the Bartlett-Lewis Optimization toolbox and its subfolders to
%    Matlab's path

 addpath(genpath('D:\Users\kdroos\Documents\MATLAB\Copulas\2_Related_Functions\BLBOX'));

%% 2. Specify which model you want to calibrate
% choices are : - OBL (original Bartlett-Lewis model)
%               - MBL (modified Bartlett-Lewis model)
%               - MBLG (modified Bartlett-Lewis Gamma model)
%               - MBLP (modified Bartlett-Lewis Pareto model)
%               - TBL (Truncated Bartlett-Lewis model)
%               - TBLG (Truncated Bartlett-Lewis Gamma model)
%               - TBLP (Truncated Bartlett-Lewis Pareto model)

MODEL = 'MBL';

%% 3. Choose observation data set (make sure it is in the style of the 105-year Uccle timeseries)

OBSERVATIONS = 'ukkel72_10min.mat';

%% 4. Specify objective function
% choices are : - OBJ_fn (standard objective function)
%               - OBJ_fn_2 (implicit weights)
%               - OBJ_fn_3 (Covariance matrix as weigths)
%               - OBJ_fn_4 (Variance as weights)

OBJFN = 'OBJ_fn';

%% 5. Choose which months to calibrate, and how many times to repeat this
% calibration (0 is not an option)
% i.e. if you want to calibrate months 4 and 7, and repeat each calibration 5
% times, then MONTH = [4,7]; REPEAT = 5;.

MONTH = 1:12;
REPEAT = 20;

%% 6. Choose an optimization method
% choices are : - 'neldermead' (simplex method)
%               - 'pso' (particle swarm optimization)
%               - 'sce' (shuffled complex evolution)
%               - 'simpsa' (simplex simulated annealing)

ROUTINE = 'SCE'; 

% if you want to CHANGE the OPTIMIZATION METHOD'S SETTINGS, 
% open getOptions.m and edit that function.

%% 7. Set boundaries of the parameter space by opening getBoundaries.m and making changes

%% 8. Set properties to be used by the generalized method of moments by opening getRules and making changes

%% 9. Start the fitting procedure 

% 9b. (optional) specify where to save output
SAVE = 'CalStr_Boundaries_ukkel_10min';

zdpflag = 0;
output=momFit(MODEL,OBSERVATIONS,OBJFN,ROUTINE,MONTH,REPEAT,getRules,getBoundaries(MODEL),SAVE,[],zdpflag);



