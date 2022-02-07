function out = sampleCvine(X,W,family_or,theta, copMatrix)
%% sampleCvine Sample from a C-vine
%   Implements a sampling algorithm for C-vines. Based on, but slightly
%   different than the sampling algorithm in Aas et al. (2009).
%
%   Inputs:
%       X: matrix with the known tree values
%       W: random generator
%       family: cell with the used families
%       theta: cell with the used families
%   Output:
%       Out: sample of the C-vine
%
%   Originally implemented by M. T. Pham and H. Vernieuwe
%   Last update by J. Van de Velde on 14/07/'21: fixed last bugs

%% Check X

for i=1:length(X)
    if X(i) == 1
        X(i) = 1-rand(1)*0.0001;
    elseif X(i) == 0
        X(i) = 0 + rand(1)*0.0001;
    end
end  

%% Reconstruction of matrices

if size(family_or,2) == 3
    theta_tmp = [theta(2), theta(3); nan theta(6)]';
    X(3) = NaN;
    Xnew = [X(copMatrix(1,3,1)), X(copMatrix(1,2,1)), X(copMatrix(1,1,1))];
    family = {family_or{:,2,3} []; family_or{:,1,3} family_or{:,1,2} };
    theta = [theta_tmp(2,2) nan; theta_tmp(2,1) theta_tmp(1,1)];
else
    theta_tmp = [theta(2), theta(3) theta(4); nan theta(7) theta(8); nan nan theta(12)]';
    X(4) = NaN;
    Xnew = [X(copMatrix(1,4,1)),  X(copMatrix(1,3,1)), X(copMatrix(1,2,1)), X(copMatrix(1,1,1))];
    family = {family_or{:,3,4} [] []; family_or{:,2,4} family_or{:,2,3}  [];family_or{:,1,4}  family_or{:,1,3} family_or{:,1,2}   };
    theta = [theta_tmp(3,3) nan nan; theta_tmp(3,2) theta_tmp(2,2)  nan; theta_tmp(3,1) theta_tmp(2,1) theta_tmp(1,1) ];
end

%% Simulation algorithm
n = length(W);
out = nan(n,1);
Vine = nan(size(family_or,2));
locNaN= nan(size(Xnew,2)-1,1);
% Loop over vine copula
Vine(:,1) = Xnew(1,:);
for nVars = 2:size(Xnew,2)
    for nTree = 1:nVars-1
        if isnan(Vine(nTree,nTree))
            locNaN(nVars-1)=nVars;
            break
        elseif isnan(Vine(nVars,nTree))
            locNaN(nTree)=nVars;
            if nTree == size(Xnew,2)-1 %Should be nVars-1 if sampling flexibly
                Vine(nTree+1,nTree+1)=W(1);
                for k=nVars-1:-1:1 % Separate function?
                    Vine(locNaN(k),k) = h(family{locNaN(k)-1,k},Vine(locNaN(k),k+1),Vine(k,k),theta(locNaN(k)-1,k),'U','inv');
                    out = Vine(locNaN(k),k);
                end
            end
        else
            Vine(nVars,nTree+1) = h(family{nVars-1,nTree},Vine(nTree,nTree),Vine(nVars,nTree),theta(nVars-1,nTree),'U','norm');
        end
    end
    if ~isnan(out)
        break
    end
end

end
