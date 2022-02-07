function prepareTPEpE(data,modelname)
%% prepareTPEpE This function builds the monthly C-vines for the configuration TPEpE
%   Set the data in the correct form to build the monthly C-vines
%   Monthly datasets will be built as follows: 3-tuples (T(j+1), P(j+1), E(j+1), E(j))
%
%   Format of observation dataset:
%   year month day T P E
%   Noise is already added
%
%   Inputs:
%       data: date used for the fit
%       modelname: a base file name used for saves
% 
%   Originally implemented by M. T. Pham and H. Vernieuwe (?)
%   Last updated by J. Van de Velde on 08/04/'21: documentation

%% Initialization
DOM = [31 28 31 30 31 30 31 31 30 31 30 31; % Normal year
    31 29 31 30 31 30 31 31 30 31 30 31]; %Leap year
startyear = min(data(:,1));
endyear = max(data(:,1));
nyears = endyear-startyear+1;
YI = 1:1:nyears; %Year index
YIfeb = YI;
yearnrs = [startyear YI+startyear-1];
I = rem(yearnrs,4)==0 & (rem(yearnrs,100)~=0|rem(yearnrs,400)==0); % leap year

%% Preparation loop

for i = 1:12
    dataMonth = data(data(:,2)==i,4:6); %Selection of data for the current month (out of all years)
    
    %For each month, an Xst and Xend (for Ep) and Yst and Yend (for the
    %other variables) are calculated.
    if i == 2
        %February index changes depending on the leap year
        YIfeb(I) = DOM(2,i);
        YIfeb(~I) = DOM(1,i);
        %Further setup
        Xst=[1 cumsum(YIfeb(1:end-1))+1]; %Creates vector with size nyr
        Xend = Xst(2:end)-2; %Vector
        Xst = Xst(1:end-1);
        Yst = Xst+1;
        Yend = Xend+1;
    else
        Xst = (DOM(1,i)*YI-DOM(1,i))+1;
        Xend= DOM(1,i)*YI-1;
        Xend(end)=Xend(end)-1;
        Yst= Xst+1;
        Yend = Xend+1;
    end
    
    X=[];
    Y=[];
    Z=[];
    W=[];
    for k =  1:length(Xst)
        %Making a monthly vector for each year, puts the vectors under each
        %other
        X = [X; dataMonth(Yst(k):Yend(k),1)];       
        Y = [Y; dataMonth(Yst(k):Yend(k),2)];
        Z = [Z; dataMonth(Yst(k):Yend(k),3)];
        W = [W; dataMonth(Xst(k):Xend(k),3)]; 
    end
    filename = sprintf('TPEpE_%d.mat',i); %Different filename for each month
    D = [X Y W Z];
    save(['E:\Users\jpvdveld\Onderzoek\Data\2_detrended+vines\' modelname '_' filename],'D');
end

end



