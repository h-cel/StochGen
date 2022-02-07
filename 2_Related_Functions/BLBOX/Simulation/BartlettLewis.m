function synthrain=BartlettLewis(param,ybegin,yend,model,outputfile)
%% BartlettLewis Simulates the Bartlett-Lewis Rectangular Pulses models
% 
%   The Bartlett-Lewis Rectangular Pulses modelssimulate temporal rainfall as a series of clustered rectangular pulses. 
%   Each of these clusters represent a single rainfall event, whose
%   structure and intensity is characterized by a set of rectangular pulses. 
%
%   synthrain=BartlettLewis(param,ybegin,yend,model) generates a timeseries
%   aggregated over 10 min intervals between year ybegin and yend. ybegin
%   and yend are both scalar. param is the parameter matrix, consisting of
%   12 rows, representing the different months. model is a string which
%   specifies which type of model is to be used, this can be:
%
%       'OBL'   Original Bartlett-Lewis model           (param=12x5)
%       'MBL'   Modified Bartlett-Lewis model           (param=12x6)
%       'RBL2'   Randomized Bartlett-Lewis model 2      (param=12x7)
%       'MBLG'  Modified Bartlett-Lewis Gamma model     (param=12x7)
%       'TBL'   Truncated Bartlett-Lewis model          (param=12x7)
%       'TBLG'  Truncated Bartlett-Lewis Gamma model    (param=12x8)
%       'MBLP'  Modified Bartlett-Lewis Pareto model    (param=12x7)
%       'TBLP'  Modified Bartlett-Lewis Pareto model    (param=12x8)
%
%   synthrain=BartlettLewis(param,ybegin,yend,model,outputfile) does
%   exactly the same, and writes to output to the specified outputfile.
%   ouputfile is a string.
%
%   Last updated by J. Van de Velde on 12/03/'21

%% Controlled randomization

rand('state',sum(100*clock))
%% Model-input: parameters
if ischar(param)==1
    param=matload(param);
end
% ybegin = Begin year of simulation
% yend = End year of simulation

%% Simulation parameters
years=ybegin:1:yend;
nry=length(years); %number of years
aggl=1/6; %10-min aggregation level
simd=zeros(nry,1);
nti=zeros(nry,1);
days=zeros(nry,12);
months=zeros(nry,12);

for y=1:nry
    year=years(y);
    %Taking leap years into account during simulation
    if (rem(year,4)==0 && rem(year,100)~=0) || (rem(year,400)==0)
        months(y,:)=cumsum([31 29 31 30 31 30 31 31 30 31 30 31]*24);
        simd(y)=(29+4*30+7*31)*24;
        nti(y)=simd(y)/aggl;
        days(y,:)=[31 29 31 30 31 30 31 31 30 31 30 31];
    else
        months(y,:)=cumsum([31 28 31 30 31 30 31 31 30 31 30 31]*24);
        simd(y)=(28+4*30+7*31)*24;
        nti(y)=simd(y)/aggl;
        days(y,:)=[31 28 31 30 31 30 31 31 30 31 30 31];
    end
end
cumnti=cumsum(nti);

%% Output file
synthrain=zeros(sum(nti)*aggl,1/aggl); %10-min rainfall-amounts, 1 row=1hour

%Make matrix with headers for 10-min rainfall (Uccle-data style)
header=zeros(sum(nti)*aggl,7);
count=0;
row=1;
for y=1:nry
    year=years(y);
    for m=1:12
        for d=1:days(y,m)
            for i=1:24
            hour=0:1:23;
            header(row,:)=[6 6100 year m d hour(i) 1];
            row=row+1;
            end
            count=count+24;
        end
    end
end


%% OBL model
switch lower(model)
    case 'obl'
%% Simulation

intensnext=zeros(1,nti(1));
%rand('twister',-9999);
for y=1:nry
    
    year=years(y);
    % Initialisation each year
    bts=0; %Begin time storm [h]
    btc=0; %Begin time cell [h]
    lambda=param(1,1); %Lambda for January
    
    intens=intensnext;  %intensities/aggregationlevel/year
                        %When cell durations are larger than a year,
                        %intensities are stored for the next year
    
    if y<nry %intensity array for next year with appropriate length
    intensnext=zeros(1,nti(y+1)); 
    else
    intensnext=zeros(1,nti(y));
    end
    
    % First storm origin
    bts=bts+randexp(lambda);

    while bts < simd(y)

        month=find(months(y,:)>=bts,1,'first'); %month in which storm begins

        lambda=param(month,1);
        beta=param(month,2);
        gamma=param(month,3);
        muix=param(month,4);
        eta=param(month,5);

        %Storm duration
        ds=randexp(gamma); %duration of the storm
        ets=bts+ds; %end time the storm

        if ets>simd(y)
            ets=simd(y);
        end

        btc=bts; %first cell at storm origin

        while btc < ets %Generate cells as long as beginning is in storm

            %Cells

            dc=randexp(eta); %duration of cell
            etc=btc+dc; %end time of cell

            %Cell intensity
            intcell=randexp(1/muix);
            intcell=intcell*aggl; %cell intensity from mm/h --> mm/agglevel
            bir=btc/aggl; %cell begin in unit of agglevel
            eir=etc/aggl; %cell end in unit of agglevel
            a=1-(bir-floor(bir));
            b=eir-floor(eir);
            bi=floor(bir)+1;
            ei=floor(eir)+1;

            if floor(bir) == floor(eir)
                intens(bi)=intens(bi)+(eir-bir)*intcell;
            elseif ei >nti(y)
                intens(bi)=intens(bi)+a*intcell;
                intensnext(ei-nti(y)) = intensnext(ei-nti(y))+b*intcell;
                for ti=bi+1:nti(y)
                    intens(ti)=intens(ti)+intcell;
                end
                for ti=1:ei-nti(y)-1
                    intensnext(ti) = intensnext(ti)+intcell;
                end
            else
                intens(bi)=intens(bi)+a*intcell;
                intens(ei)=intens(ei)+b*intcell;
                for ti=bi+1:ei-1
                    intens(ti)=intens(ti)+intcell;
                end
            end

            %Begin timing next cell
            btc=btc+randexp(beta);
        end
        
        bts=bts+randexp(lambda); %begin time of a storm
    end

    
    intens=reshape(intens',1/aggl,nti(y)*aggl)';
   
    synthrain((cumnti(y)*aggl-nti(y)*aggl+1):cumnti(y)*aggl,:)=intens;
    

end


%Save output

synthrain=[header synthrain];
if nargin==5
    save(outputfile,'synthrain');
end
%% MBL model
    case 'mbl'
%% Simulation

intensnext=zeros(1,nti(1));
%rand('twister',-9999);
for y=1:nry
    year=years(y);
    % Initialisation each year
    bts=0; %Begin time storm [h]
    btc=0; %Begin time cell [h]
    lambda=param(1,1); %Lambda for January
    
    intens=intensnext;  %intensities/aggregationlevel/year
                        %When cell durations are larger than a year,
                        %intensities are stored for the next year
    
    if y<nry %intensity array for next year with appropriate length
    intensnext=zeros(1,nti(y+1)); 
    else
    intensnext=zeros(1,nti(y));
    end
    
    % First storm origin
    bts=bts+randexp(lambda);

    while bts < simd(y)

        month=find(months(y,:)>=bts,1,'first'); %month in which storm begins

        lambda=param(month,1);
        kappa=param(month,2);
        phi=param(month,3);
        muix=param(month,4);
        alpha=param(month,5);
        theta=1/param(month,6);

        %Storm duration
        eta=gamrnd(alpha,theta);
        gamma=phi*eta;
        
        ds=randexp(gamma); %duration of the storm
        ets=bts+ds; %end time the storm

        if ets>simd(y)
            ets=simd(y);
        end

        btc=bts; %first cell at storm origin

        while btc < ets %Generate cells as long as beginning is in storm

            %Cells
           
            dc=randexp(eta); %duration of cell
            etc=btc+dc; %end time of cell

            %Cell intensity
            intcell=randexp(1/muix);
 
            intcell=intcell*aggl; %cell intensity from mm/h --> mm/agglevel
            bir=btc/aggl; %cell begin in unit of agglevel
            eir=etc/aggl; %cell end in unit of agglevel
            a=1-(bir-floor(bir));
            b=eir-floor(eir);
            bi=floor(bir)+1;
            ei=floor(eir)+1;

            if floor(bir) == floor(eir)
                intens(bi)=intens(bi)+(eir-bir)*intcell;
            elseif ei >nti(y)
                intens(bi)=intens(bi)+a*intcell;
                intensnext(ei-nti(y)) = intensnext(ei-nti(y))+b*intcell;
                for ti=bi+1:nti(y)
                    intens(ti)=intens(ti)+intcell;
                end
                for ti=1:ei-nti(y)-1
                    intensnext(ti) = intensnext(ti)+intcell;
                end
            else
                intens(bi)=intens(bi)+a*intcell;
                intens(ei)=intens(ei)+b*intcell;
                for ti=bi+1:ei-1
                    intens(ti)=intens(ti)+intcell;
                end
            end

            %Begin timing next cell
            beta=kappa*eta;
            btc=btc+randexp(beta);
        end
        
        bts=bts+randexp(lambda); %begin time of a storm
    end

    intens=reshape(intens',1/aggl,nti(y)*aggl)';
   
    synthrain((cumnti(y)*aggl-nti(y)*aggl+1):cumnti(y)*aggl,:)=intens;
    

end


%Save output

synthrain=[header synthrain];
if nargin==5
    save(outputfile,'synthrain');
end

%% RBL2 model
    case 'rbl2'
%% Simulation

intensnext=zeros(1,nti(1));
%rand('twister',-9999);
for y=1:nry
    year=years(y);
    % Initialisation each year
    bts=0; %Begin time storm [h]
    btc=0; %Begin time cell [h]
    lambda=param(1,1); %Lambda for January
    
    intens=intensnext;  %intensities/aggregationlevel/year
                        %When cell durations are larger than a year,
                        %intensities are stored for the next year
    
    if y<nry %intensity array for next year with appropriate length
    intensnext=zeros(1,nti(y+1)); 
    else
    intensnext=zeros(1,nti(y));
    end
    
    % First storm origin
    bts=bts+randexp(lambda);

    while bts < simd(y)

        month=find(months(y,:)>=bts,1,'first'); %month in which storm begins

        lambda=param(month,1);
        nu=param(month,2);
        alpha=param(month,3);
        iota=param(month,4);
        phi=param(month,5);
        kappa=param(month,6);
        omega = param(month,7);

        %Storm duration
        eta=gamrnd(alpha,1/nu);
        gamma=phi*eta;
        
        ds=randexp(gamma); %duration of the storm
        ets=bts+ds; %end time the storm

        if ets>simd(y)
            ets=simd(y);
        end

        btc=bts; %first cell at storm origin

        while btc < ets %Generate cells as long as beginning is in storm

            %Cells
           
            
            ei = length(intensnext)*10;
            while (ei-nti(y))> length(intensnext)% Sometimes, the simulated value of dc is too large
                dc=randexp(eta); %duration of cell
                etc=btc+dc; %end time of cell
                eir=etc/aggl; %cell end in unit of agglevel
                ei=floor(eir)+1;
            end

            %Cell intensity
            intcell=gamrnd(omega, iota*eta/omega);
 
            intcell=intcell*aggl; %cell intensity from mm/h --> mm/agglevel
            bir=btc/aggl; %cell begin in unit of agglevel
            
            a=1-(bir-floor(bir));
            b=eir-floor(eir);
            if bir <0
                bi=ceil(bir)+1;
            else
                bi=floor(bir)+1;
            end
            

            if floor(bir) == floor(eir)
                intens(bi)=intens(bi)+(eir-bir)*intcell;
            elseif ei >nti(y)
                intens(bi)=intens(bi)+a*intcell;
                intensnext(ei-nti(y)) = intensnext(ei-nti(y))+b*intcell;
                for ti=bi+1:nti(y)
                    intens(ti)=intens(ti)+intcell;
                end
                for ti=1:ei-nti(y)-1
                    intensnext(ti) = intensnext(ti)+intcell;
                end
            else
                intens(bi)=intens(bi)+a*intcell;
                intens(ei)=intens(ei)+b*intcell;
                for ti=bi+1:ei-1
                    intens(ti)=intens(ti)+intcell;
                end
            end

            %Begin timing next cell
            beta=kappa*eta;
            btc=btc+randexp(beta);
        end
        
        bts=bts+randexp(lambda); %begin time of a storm
    end

    intens=reshape(intens',1/aggl,nti(y)*aggl)';
   
    synthrain((cumnti(y)*aggl-nti(y)*aggl+1):cumnti(y)*aggl,:)=intens;
    

end


%Save output

synthrain=[header synthrain];
if nargin==5
    save(outputfile,'synthrain');
end
%% MBLG model
     case 'mblg'
%% Simulation

intensnext=zeros(1,sum(nti));
%rand('twister',-9999);
for y=1:nry
    year=years(y);
    % Initialisation each year
    bts=0; %Begin time storm [h]
    btc=0; %Begin time cell [h]
    lambda=param(1,1); %Lambda for January
    
    intens=intensnext(1:nti(y));  %intensities/aggregationlevel/year
                        %When cell durations are larger than a year,
                        %intensities are stored for the next year
   
    if y<nry %intensity array for next year with appropriate length
        intensnext=intensnext(nti(y)+1:end); 
    else
        intensnext=zeros(1,nti(y));
    end
    
    % First storm origin
    bts=bts+randexp(lambda);

    while bts < simd(y)

        month=find(months(y,:)>=bts,1,'first'); %month in which storm begins

        lambda=param(month,1);
        kappa=param(month,2);
        phi=param(month,3);
        alpha=param(month,4);
        theta=1/param(month,5);
        p=param(month,6);
        delta=1/param(month,7);

        %Storm duration
        eta=gamrnd(alpha,theta);
        gamma=phi*eta;
        
        ds=randexp(gamma); %duration of the storm
        ets=bts+ds; %end time the storm

        if ets>simd(y)
            ets=simd(y);
        end

        btc=bts; %first cell at storm origin

        while btc < ets %Generate cells as long as beginning is in storm

            %Cells

            dc=randexp(eta); %duration of cell
            etc=btc+dc; %end time of cell

            %Cell intensity
            intcell=gamrnd(p,delta);
            intcell=intcell*aggl; %cell intensity from mm/h --> mm/agglevel
            bir=btc/aggl; %cell begin in unit of agglevel
            eir=etc/aggl; %cell end in unit of agglevel
            a=1-(bir-floor(bir));
            b=eir-floor(eir);
            bi=floor(bir)+1;
            ei=floor(eir)+1;

            if floor(bir) == floor(eir)
                intens(bi)=intens(bi)+(eir-bir)*intcell;
            elseif ei >nti(y)
                intens(bi)=intens(bi)+a*intcell;
                if ei-nti(y)>length(intensnext)
                    intensnext(end)=intensnext(end)+b*intcell;
                    for ti=1:length(intensnext)-1
                        intensnext(ti) = intensnext(ti)+intcell;
                    end
                else
                    intensnext(ei-nti(y)) = intensnext(ei-nti(y))+b*intcell;
                    for ti=1:ei-nti(y)-1
                        intensnext(ti) = intensnext(ti)+intcell;
                    end
                end
                for ti=bi+1:nti(y)
                    intens(ti)=intens(ti)+intcell;
                end                
            else
                intens(bi)=intens(bi)+a*intcell;
                intens(ei)=intens(ei)+b*intcell;
                for ti=bi+1:ei-1
                    intens(ti)=intens(ti)+intcell;
                end
            end


            %Begin timing next cell
            beta=kappa*eta;
            btc=btc+randexp(beta);
        end
        
        bts=bts+randexp(lambda); %begin time of a storm
    end

    
    intens=reshape(intens',1/aggl,nti(y)*aggl)';
   
    synthrain((cumnti(y)*aggl-nti(y)*aggl+1):cumnti(y)*aggl,:)=intens;
    

end



%Save output

synthrain=[header synthrain];
if nargin==5
    save(outputfile,'synthrain');
end

    case 'tbl'
%% TBL model
%% Simulation

intensnext=zeros(1,nti(1));
%rand('twister',-9999);
% h = waitbar(0,'Please wait...simulation in progress...');
    
for y=1:nry,
    disp(y);
% waitbar(y/nry)


    year=years(y);
    % Initialisation each year
    bts=0; %Begin time storm [h]
    btc=0; %Begin time cell [h]
    lambda=param(1,1); %Lambda for January
    
    intens=intensnext;  %intensities/aggregationlevel/year
                        %When cell durations are larger than a year,
                        %intensities are stored for the next year
    
    if y<nry %intensity array for next year with appropriate length
    intensnext=zeros(1,nti(y+1)); 
    else
    intensnext=zeros(1,nti(y));
    end
    
    % First storm origin
    bts=bts+randexp(lambda);

    while bts < simd(y)

        month=find(months(y,:)>=bts,1,'first'); %month in which storm begins

        lambda=param(month,1);
        kappa=param(month,2);
        phi=param(month,3);
        muix=param(month,4);
        alpha=param(month,5);
        theta=1/param(month,6);
        eps=param(month,7);
        %Storm duration
        
        eta=gamrnd(alpha,theta);
        while eta<eps
            eta=gamrnd(alpha,theta);
        end
        gamma=phi*eta;
        
        ds=randexp(gamma); %duration of the storm
        ets=bts+ds; %end time the storm

        if ets>simd(y)
            ets=simd(y);
        end

        btc=bts; %first cell at storm origin

        while btc < ets %Generate cells as long as beginning is in storm

            %Cells
           
            dc=randexp(eta); %duration of cell
            etc=btc+dc; %end time of cell

            %Cell intensity
            intcell=randexp(1/muix);
 
            intcell=intcell*aggl; %cell intensity from mm/h --> mm/agglevel
            bir=btc/aggl; %cell begin in unit of agglevel
            eir=etc/aggl; %cell end in unit of agglevel
            a=1-(bir-floor(bir));
            b=eir-floor(eir);
            bi=floor(bir)+1;
            ei=floor(eir)+1;

            if floor(bir) == floor(eir)
                intens(bi)=intens(bi)+(eir-bir)*intcell;
            elseif ei >nti(y)
                intens(bi)=intens(bi)+a*intcell;
                intensnext(ei-nti(y)) = intensnext(ei-nti(y))+b*intcell;
                for ti=bi+1:nti(y)
                    intens(ti)=intens(ti)+intcell;
                end
                for ti=1:ei-nti(y)-1
                    intensnext(ti) = intensnext(ti)+intcell;
                end
            else
                intens(bi)=intens(bi)+a*intcell;
                intens(ei)=intens(ei)+b*intcell;
                for ti=bi+1:ei-1
                    intens(ti)=intens(ti)+intcell;
                end
            end

            %Begin timing next cell
            beta=kappa*eta;
            btc=btc+randexp(beta);
        end
        
        bts=bts+randexp(lambda); %begin time of a storm
    end

    
    intens=reshape(intens',1/aggl,nti(y)*aggl)';
   
    synthrain((cumnti(y)*aggl-nti(y)*aggl+1):cumnti(y)*aggl,:)=intens;
    

end


%Save output

synthrain=[header synthrain];
if nargin==5
    save(outputfile,'synthrain');
end
% close(h) 




%% TBLG model
    case 'tblg'

     %% Simulation

intensnext=zeros(1,nti(1));
%rand('twister',-9999);
for y=1:nry
    year=years(y);
    disp(y);
    % Initialisation each year
    bts=0; %Begin time storm [h]
    btc=0; %Begin time cell [h]
    lambda=param(1,1); %Lambda for January
    
    intens=intensnext;  %intensities/aggregationlevel/year
                        %When cell durations are larger than a year,
                        %intensities are stored for the next year
    
    if y<nry %intensity array for next year with appropriate length
    intensnext=zeros(1,nti(y+1)); 
    else
    intensnext=zeros(1,nti(y));
    end
    
    % First storm origin
    bts=bts+randexp(lambda);

    while bts < simd(y)

        month=find(months(y,:)>=bts,1,'first'); %month in which storm begins

        lambda=param(month,1);
        kappa=param(month,2);
        phi=param(month,3);
        alpha=param(month,4);
        theta=1/param(month,5);
        p=param(month,6);
        delta=1/param(month,7);
        eps=param(month,8);

        %Storm duration
        eta=gamrnd(alpha,theta);
        while eta<eps
            eta=gamrnd(alpha,theta);
        end
        gamma=phi*eta;
        
        ds=randexp(gamma); %duration of the storm
        ets=bts+ds; %end time the storm

        if ets>simd(y)
            ets=simd(y);
        end

        btc=bts; %first cell at storm origin

        while btc < ets %Generate cells as long as beginning is in storm

            %Cells

            dc=randexp(eta); %duration of cell
            etc=btc+dc; %end time of cell

            %Cell intensity
            intcell=gamrnd(p,delta);
            intcell=intcell*aggl; %cell intensity from mm/h --> mm/agglevel
            bir=btc/aggl; %cell begin in unit of agglevel
            eir=etc/aggl; %cell end in unit of agglevel
            a=1-(bir-floor(bir));
            b=eir-floor(eir);
            bi=floor(bir)+1;
            ei=floor(eir)+1;

            if floor(bir) == floor(eir)
                intens(bi)=intens(bi)+(eir-bir)*intcell;
            elseif ei >nti(y)
                intens(bi)=intens(bi)+a*intcell;
                intensnext(ei-nti(y)) = intensnext(ei-nti(y))+b*intcell;
                for ti=bi+1:nti(y)
                    intens(ti)=intens(ti)+intcell;
                end
                for ti=1:ei-nti(y)-1
                    intensnext(ti) = intensnext(ti)+intcell;
                end
            else
                intens(bi)=intens(bi)+a*intcell;
                intens(ei)=intens(ei)+b*intcell;
                for ti=bi+1:ei-1
                    intens(ti)=intens(ti)+intcell;
                end
            end

            %Begin timing next cell
            beta=kappa*eta;
            btc=btc+randexp(beta);
        end
        
        bts=bts+randexp(lambda); %begin time of a storm
    end

    
    intens=reshape(intens',1/aggl,nti(y)*aggl)';
   
    synthrain((cumnti(y)*aggl-nti(y)*aggl+1):cumnti(y)*aggl,:)=intens;
    

end


%Save output

synthrain=[header synthrain];
if nargin==5
    save(outputfile,'synthrain');
end
%% TBLP model
    case 'tblp'
%% Simulation

intensnext=zeros(1,nti(1));
%rand('twister',-9999);
for y=1:nry
    disp(y);
    year=years(y);
    % Initialisation each year
    bts=0; %Begin time storm [h]
    btc=0; %Begin time cell [h]
    lambda=param(1,1); %Lambda for January
    
    intens=intensnext;  %intensities/aggregationlevel/year
                        %When cell durations are larger than a year,
                        %intensities are stored for the next year
    
    if y<nry %intensity array for next year with appropriate length
    intensnext=zeros(1,nti(y+1)); 
    else
    intensnext=zeros(1,nti(y));
    end
    
    % First storm origin
    bts=bts+randexp(lambda);

    while bts < simd(y)

        month=find(months(y,:)>=bts,1,'first'); %month in which storm begins

        lambda=param(month,1);
        kappa=param(month,2);
        phi=param(month,3);
        alpha=param(month,4);
        theta=1/param(month,5);
        p=param(month,6);
        delta=param(month,7);
        eps=param(month,8);

        %Storm duration
        eta=gamrnd(alpha,theta);
        while eta<eps
            eta=gamrnd(alpha,theta);
        end
        gamma=phi*eta;
        
        ds=randexp(gamma); %duration of the storm
        ets=bts+ds; %end time the storm

        if ets>simd(y)
            ets=simd(y);
        end

        btc=bts; %first cell at storm origin

        while btc < ets %Generate cells as long as beginning is in storm

            %Cells

            dc=randexp(eta); %duration of cell
            etc=btc+dc; %end time of cell

            %Cell intensity
            intcell=gprnd(p,delta,delta/p);
            intcell=intcell*aggl; %cell intensity from mm/h --> mm/agglevel
            bir=btc/aggl; %cell begin in unit of agglevel
            eir=etc/aggl; %cell end in unit of agglevel
            a=1-(bir-floor(bir));
            b=eir-floor(eir);
            bi=floor(bir)+1;
            ei=floor(eir)+1;

            if floor(bir) == floor(eir)
                intens(bi)=intens(bi)+(eir-bir)*intcell;
            elseif ei >nti(y)
                intens(bi)=intens(bi)+a*intcell;
                intensnext(ei-nti(y)) = intensnext(ei-nti(y))+b*intcell;
                for ti=bi+1:nti(y)
                    intens(ti)=intens(ti)+intcell;
                end
                for ti=1:ei-nti(y)-1
                    intensnext(ti) = intensnext(ti)+intcell;
                end
            else
                intens(bi)=intens(bi)+a*intcell;
                intens(ei)=intens(ei)+b*intcell;
                for ti=bi+1:ei-1
                    intens(ti)=intens(ti)+intcell;
                end
            end

            %Begin timing next cell
            beta=kappa*eta;
            btc=btc+randexp(beta);
        end
        
        bts=bts+randexp(lambda); %begin time of a storm
    end

    
    intens=reshape(intens',1/aggl,nti(y)*aggl)';
   
    synthrain((cumnti(y)*aggl-nti(y)*aggl+1):cumnti(y)*aggl,:)=intens;
    

end


%Save output

synthrain=[header synthrain];
if nargin==5
    save(outputfile,'synthrain');
end
%% MBLP model
    case 'mblp'
%% Simulation

intensnext=zeros(1,nti(1));
%rand('twister',-9999);
for y=1:nry
    disp(y);
    year=years(y);
    % Initialisation each year
    bts=0; %Begin time storm [h]
    btc=0; %Begin time cell [h]
    lambda=param(1,1); %Lambda for January
    
    intens=intensnext;  %intensities/aggregationlevel/year
                        %When cell durations are larger than a year,
                        %intensities are stored for the next year
    
    if y<nry %intensity array for next year with appropriate length
    intensnext=zeros(1,nti(y+1)); 
    else
    intensnext=zeros(1,nti(y));
    end
    
    % First storm origin
    bts=bts+randexp(lambda);

    while bts < simd(y)

        month=find(months(y,:)>=bts,1,'first'); %month in which storm begins

        lambda=param(month,1);
        kappa=param(month,2);
        phi=param(month,3);
        alpha=param(month,4);
        theta=1/param(month,5);
        p=param(month,6);
        delta=param(month,7);
        

        %Storm duration
        eta=gamrnd(alpha,theta);
        gamma=phi*eta;
        
        ds=randexp(gamma); %duration of the storm
        ets=bts+ds; %end time the storm

        if ets>simd(y)
            ets=simd(y);
        end

        btc=bts; %first cell at storm origin

        while btc < ets %Generate cells as long as beginning is in storm

            %Cells

            dc=randexp(eta); %duration of cell
            etc=btc+dc; %end time of cell

            %Cell intensity
            intcell=gprnd(p,delta,delta/p);
            intcell=intcell*aggl; %cell intensity from mm/h --> mm/agglevel
            bir=btc/aggl; %cell begin in unit of agglevel
            eir=etc/aggl; %cell end in unit of agglevel
            a=1-(bir-floor(bir));
            b=eir-floor(eir);
            bi=floor(bir)+1;
            ei=floor(eir)+1;

            if floor(bir) == floor(eir)
                intens(bi)=intens(bi)+(eir-bir)*intcell;
            elseif ei >nti(y)
                intens(bi)=intens(bi)+a*intcell;
                intensnext(ei-nti(y)) = intensnext(ei-nti(y))+b*intcell;
                for ti=bi+1:nti(y)
                    intens(ti)=intens(ti)+intcell;
                end
                for ti=1:ei-nti(y)-1
                    intensnext(ti) = intensnext(ti)+intcell;
                end
            else
                intens(bi)=intens(bi)+a*intcell;
                intens(ei)=intens(ei)+b*intcell;
                for ti=bi+1:ei-1
                    intens(ti)=intens(ti)+intcell;
                end
            end

            %Begin timing next cell
            beta=kappa*eta;
            btc=btc+randexp(beta);
        end
        
        bts=bts+randexp(lambda); %begin time of a storm
    end

    
    intens=reshape(intens',1/aggl,nti(y)*aggl)';
   
    synthrain((cumnti(y)*aggl-nti(y)*aggl+1):cumnti(y)*aggl,:)=intens;
    

end


%Save output

synthrain=[header synthrain];
% if nargin==5
%     save(outputfile,'synthrain');
% end

%% TBLG model
    case 'gpdbl'

     %% Simulation

intensnext=zeros(1,nti(1));
%rand('twister',-9999);
for y=1:nry
    year=years(y);

    % Initialisation each year
    bts=0; %Begin time storm [h]
    btc=0; %Begin time cell [h]
    gpdpar1=param(1,1);
    gpdpar2=param(1,2);%Lambda for January

    intens=intensnext;  %intensities/aggregationlevel/year
                        %When cell durations are larger than a year,
                        %intensities are stored for the next year
    
    if y<nry %intensity array for next year with appropriate length
    intensnext=zeros(1,nti(y+1)); 
    else
    intensnext=zeros(1,nti(y));
    end
    
    % First storm origin
    bts=bts+gprnd(gpdpar1,gpdpar2,0);

    while bts < simd(y)

        month=find(months(y,:)>=bts,1,'first'); %month in which storm begins

        gpdpar1=param(month,1);
        gpdpar2=param(month,2);
        kappa=param(month,3);
        phi=param(month,4);
        alpha=param(month,5);
        theta=1/param(month,6);
        p=param(month,7);
        delta=1/param(month,8);
        eps=param(month,9);

        %Storm duration
        eta=gamrnd(alpha,theta);
        while eta<eps
            eta=gamrnd(alpha,theta);
        end
        gamma=phi*eta;
        
        ds=randexp(gamma); %duration of the storm
        ets=bts+ds; %end time the storm

        if ets>simd(y)
            ets=simd(y);
        end

        btc=bts; %first cell at storm origin

        while btc < ets %Generate cells as long as beginning is in storm

            %Cells

            dc=randexp(eta); %duration of cell
            etc=btc+dc; %end time of cell

            %Cell intensity
            intcell=gamrnd(p,delta);
            intcell=intcell*aggl; %cell intensity from mm/h --> mm/agglevel
            bir=btc/aggl; %cell begin in unit of agglevel
            eir=etc/aggl; %cell end in unit of agglevel
            a=1-(bir-floor(bir));
            b=eir-floor(eir);
            bi=floor(bir)+1;
            ei=floor(eir)+1;

            if floor(bir) == floor(eir)
                intens(bi)=intens(bi)+(eir-bir)*intcell;
            elseif ei >nti(y)
                intens(bi)=intens(bi)+a*intcell;
                intensnext(ei-nti(y)) = intensnext(ei-nti(y))+b*intcell;
                for ti=bi+1:nti(y)
                    intens(ti)=intens(ti)+intcell;
                end
                for ti=1:ei-nti(y)-1
                    intensnext(ti) = intensnext(ti)+intcell;
                end
            else
                intens(bi)=intens(bi)+a*intcell;
                intens(ei)=intens(ei)+b*intcell;
                for ti=bi+1:ei-1
                    intens(ti)=intens(ti)+intcell;
                end
            end

            %Begin timing next cell
            beta=kappa*eta;
            btc=btc+randexp(beta);
        end
        
        bts=bts+gprnd(gpdpar1,gpdpar2,0); %begin time of a storm
    end

    
    intens=reshape(intens',1/aggl,nti(y)*aggl)';
   
    synthrain((cumnti(y)*aggl-nti(y)*aggl+1):cumnti(y)*aggl,:)=intens;
    

end


%Save output

synthrain=[header synthrain];
if nargin==5
    save(outputfile,'synthrain');
end

%% TBLG model
    case 'gambl'

     %% Simulation

intensnext=zeros(1,nti(1));
%rand('twister',-9999);
for y=1:nry
    year=years(y);

    % Initialisation each year
    bts=0; %Begin time storm [h]
    btc=0; %Begin time cell [h]
    gampar1=param(1,1);
    gampar2=param(1,2);%Lambda for January

    intens=intensnext;  %intensities/aggregationlevel/year
                        %When cell durations are larger than a year,
                        %intensities are stored for the next year
    
    if y<nry %intensity array for next year with appropriate length
    intensnext=zeros(1,nti(y+1)); 
    else
    intensnext=zeros(1,nti(y));
    end
    
    % First storm origin
    bts=bts+gamrnd(gampar1,gampar2);

    while bts < simd(y)

        month=find(months(y,:)>=bts,1,'first'); %month in which storm begins

        gampar1=param(month,1);
        gampar2=param(month,2);
        kappa=param(month,3);
        phi=param(month,4);
        alpha=param(month,5);
        theta=1/param(month,6);
        p=param(month,7);
        delta=1/param(month,8);
        eps=param(month,9);

        %Storm duration
        eta=gamrnd(alpha,theta);
        while eta<eps
            eta=gamrnd(alpha,theta);
        end
        gamma=phi*eta;
        
        ds=randexp(gamma); %duration of the storm
        ets=bts+ds; %end time the storm

        if ets>simd(y)
            ets=simd(y);
        end

        btc=bts; %first cell at storm origin

        while btc < ets %Generate cells as long as beginning is in storm

            %Cells

            dc=randexp(eta); %duration of cell
            etc=btc+dc; %end time of cell

            %Cell intensity
            intcell=gamrnd(p,delta);
            intcell=intcell*aggl; %cell intensity from mm/h --> mm/agglevel
            bir=btc/aggl; %cell begin in unit of agglevel
            eir=etc/aggl; %cell end in unit of agglevel
            a=1-(bir-floor(bir));
            b=eir-floor(eir);
            bi=floor(bir)+1;
            ei=floor(eir)+1;

            if floor(bir) == floor(eir)
                intens(bi)=intens(bi)+(eir-bir)*intcell;
            elseif ei >nti(y)
                intens(bi)=intens(bi)+a*intcell;
                intensnext(ei-nti(y)) = intensnext(ei-nti(y))+b*intcell;
                for ti=bi+1:nti(y)
                    intens(ti)=intens(ti)+intcell;
                end
                for ti=1:ei-nti(y)-1
                    intensnext(ti) = intensnext(ti)+intcell;
                end
            else
                intens(bi)=intens(bi)+a*intcell;
                intens(ei)=intens(ei)+b*intcell;
                for ti=bi+1:ei-1
                    intens(ti)=intens(ti)+intcell;
                end
            end

            %Begin timing next cell
            beta=kappa*eta;
            btc=btc+randexp(beta);
        end
        
        bts=bts+gamrnd(gampar1,gampar2); %begin time of a storm
    end

    
    intens=reshape(intens',1/aggl,nti(y)*aggl)';
   
    synthrain((cumnti(y)*aggl-nti(y)*aggl+1):cumnti(y)*aggl,:)=intens;
    

end


%Save output

synthrain=[header synthrain];
if nargin==5
    save(outputfile,'synthrain');
end

%% cop3 model
    case 'cop3'

%% Simulation

intensnext=zeros(1,sum(nti));
%rand('twister',-9999);
for y=1:nry
    year=years(y);
    % Initialisation each year
    bts=0; %Begin time storm [h]
    btc=0; %Begin time cell [h]
    lambda=param(1,1); %Lambda for January
    
    intens=intensnext(1:nti(y));  %intensities/aggregationlevel/year
                        %When cell durations are larger than a year,
                        %intensities are stored for the next year
   
    if y<nry %intensity array for next year with appropriate length
        intensnext=intensnext(nti(y)+1:end); 
    else
        intensnext=zeros(1,nti(y));
    end
    
    % First storm origin
    bts=bts+randexp(lambda);

    while bts < simd(y)

        month=find(months(y,:)>=bts,1,'first'); %month in which storm begins

        lambda=param(month,1);
        kappa=param(month,2);
        phi=param(month,3);
        alpha=param(month,4);
        theta=1/param(month,5);
        p=param(month,6);
        delta=1/param(month,7);
        copparam=param(month,8);
        eps=param(month,9);

        %Storm duration
        a = gamcdf(eps,alpha,theta);
        r = a + (1-a).*rand;
        eta = gaminv(r,alpha,theta);
        gamma=phi*eta;
        
        u1=rand;        
        ds=expinv(u1,1/gamma); %duration of the storm
        ets=bts+ds; %end time the storm

        if ets>simd(y)
            ets=simd(y);
        end

        btc=bts; %first cell at storm origin

        while btc < ets %Generate cells as long as beginning is in storm

            %Cells
            u2=cond_frank(u1,copparam);
            while u2==1
                u2=cond_frank(u1,copparam);
            end
            dc=expinv(u2,1/eta); %duration of cell
            etc=btc+dc; %end time of cell

            %Cell intensity
            intcell=gamrnd(p,delta);
            intcell=intcell*aggl; %cell intensity from mm/h --> mm/agglevel
            bir=btc/aggl; %cell begin in unit of agglevel
            eir=etc/aggl; %cell end in unit of agglevel
            a=1-(bir-floor(bir));
            b=eir-floor(eir);
            bi=floor(bir)+1;
            ei=floor(eir)+1;

            if floor(bir) == floor(eir)
                intens(bi)=intens(bi)+(eir-bir)*intcell;
            elseif ei >nti(y)
                intens(bi)=intens(bi)+a*intcell;
                if ei-nti(y)>length(intensnext)
                    intensnext(end)=intensnext(end)+b*intcell;
                    for ti=1:length(intensnext)-1
                        intensnext(ti) = intensnext(ti)+intcell;
                    end
                else
                    intensnext(ei-nti(y)) = intensnext(ei-nti(y))+b*intcell;
                    for ti=1:ei-nti(y)-1
                        intensnext(ti) = intensnext(ti)+intcell;
                    end
                end
                for ti=bi+1:nti(y)
                    intens(ti)=intens(ti)+intcell;
                end                
            else
                intens(bi)=intens(bi)+a*intcell;
                intens(ei)=intens(ei)+b*intcell;
                for ti=bi+1:ei-1
                    intens(ti)=intens(ti)+intcell;
                end
            end


            %Begin timing next cell
            beta=kappa*eta;
            btc=btc+randexp(beta);
        end
        
        bts=bts+randexp(lambda); %begin time of a storm
    end

    
    intens=reshape(intens',1/aggl,nti(y)*aggl)';
   
    synthrain((cumnti(y)*aggl-nti(y)*aggl+1):cumnti(y)*aggl,:)=intens;
    

end


%Save output

synthrain=[header synthrain];
if nargin==6
    save([outputfile num2str(herhaal) '.mat'],'synthrain');
end

end
end

function u2=cond_frank(u1,alpha)
p=rand;
if abs(alpha) > log(realmax)
    u2 = (u1 < 0) + sign(alpha).*u1; % u1 or 1-u1
elseif abs(alpha) > sqrt(eps)
    u2 = -log((exp(-alpha.*u1).*(1-p)./p + exp(-alpha))./(1 + exp(-alpha.*u1).*(1-p)./p))./alpha;
%             u2 = -log(1 + (1-exp(alpha))./(exp(alpha) + exp(alpha.*(1-u1)).*(1-p)./p))./alpha;
else
    u2 = p;
end
end










    
    
    




    
