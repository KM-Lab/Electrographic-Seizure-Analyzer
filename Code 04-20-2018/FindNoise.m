function noise=FindNoise(data,posspikes,negspikes,fs,Settings,DoDisplay)
%Look through the data, using the parameters set in 'Settings', and return all noise indices.
%1. Use a sliding window to determine various metrics of the signal
%2. Compare metrics of interest to thresholds, and flag as noise if they exceed set levels
%
%Input: 'data'                Row vector containing data
%       'posspikes'           [location, height, width] info for all positive spikes, where location is index of data
%       'negspikes'           [location, height, width] info for all negative spikes, where location is index of data
%       'fs'                  Number of samples per second
%       'Settings'            Contains parameter settings. Relevant are:
%       'DoDisplay'           1 = show popup of noise identification
%                             0 = don't show popup
%
%       'noise'               Matrix with [start,end] of all noise events as indices of 'data'
%
%Calls:         none
%
%Called by:     'MenuAnalyze_Callback' in 'findallofthem.m'
%

%================================================
%Step 1: Calculate metrics on all sliding windows
%================================================
totallength=length(data);
ThrUpper=Settings.BaselineUpper; %???
ThrLower=Settings.BaselineLower; %???
WindowFactor=1;  %(Not necessary if you are using parameters in 'Settings' - Will scale sliding window extra)

%Initialize sliding window at the start of the data file
start=1;finish=Settings.movingwin(1)*fs*WindowFactor;
i=1;

%Initialize counters for efficiently calculating spike metrics in each window
filteri=1;filterj=1;filterdone=false;
filteri2=1;filterj2=1;filter2done=false;

%Intitialize data structure used to analyze various metrics 
metrics=struct;
noise=[];


%Loop through all windows, calculating and storing metrics on each window
%Note that not all metrics are used for noise calculation. Rather, we
%can consider this a place to examine the metrics and tune them with the
%observed noise patterns in the data in step 2.
while finish<totallength 
    %Extract data in current sliding window
    currentdata=data(start:finish);
    currentlength=length(currentdata);
    
    %-------------------
    %Calculate statistics on the current frame
    %-------------------
        
    metrics(i).average=mean(currentdata);       %mean of data
    metrics(i).minimum=min(currentdata);        %minimum of data
    metrics(i).maximum=max(currentdata);        %maximum of data
    metrics(i).st=std(currentdata);             %standard deviation of data
    metrics(i).absmean=mean(abs(currentdata));  %mean of the absolute value of data
    
    %Calculate the 'up-down' length of the trace (sum of vertical distances of consecutive datapoints)
    metrics(i).tracelength=sum(abs(currentdata(2:currentlength)-currentdata(1:currentlength-1)));
    
    %Normalized number of times the window data 'crosses' 0
    metrics(i).zerox=sum(currentdata(2:currentlength).*currentdata(1:currentlength-1)<0)/currentlength;

    %Normalized number of times the window data 'crosses' the window average
    metrics(i).meanx=sum((currentdata(2:currentlength)-metrics(i).average).*(currentdata(1:currentlength-1)-metrics(i).average)<0)/currentlength;
    
    %Normalized number of direction changes (number of peaks) {this metric seems effective at finding noise edges}
    rshift=currentdata(2:currentlength-1)-currentdata(1:currentlength-2); %+ for increase, - for decrease
    lshift=currentdata(2:currentlength-1)-currentdata(3:currentlength); %- for increase, + for decrease
    metrics(i).nrchanges=sum(rshift.*lshift>0)/currentlength; %count #positive products, which gives us number of peaks
      
    %Mean and Standard Deviation of positive data exceeding the set data thresholds
    ah=mean(currentdata(currentdata>ThrUpper));
    if ~isnan(ah)
        metrics(i).average_high=ah;
        metrics(i).st_high=std(currentdata(currentdata>ThrUpper));
    else
        metrics(i).average_high=0;
        metrics(i).st_high=0;
    end
    
    %Mean and Standard Deviation of negative data exceeding the data thresholds
    al=mean(currentdata(currentdata<ThrLower));
    if ~isnan(al)
        metrics(i).average_low=al;
        metrics(i).st_low=std(currentdata(currentdata<ThrLower));
    else
        metrics(i).average_low=0;
        metrics(i).st_low=0;
    end

    %The portion of the spikes that are high (measured as exceeding twice the thresholds
    pks=peakseek(currentdata, 50 ,ThrUpper);
    metrics(i).highpeakperc=(sum(currentdata(pks)>2*ThrUpper)+sum(currentdata(pks)<2*ThrLower))/currentlength;
    
    
    %mean and standard deviation of spike width and height for negative spikes
    length_negspikes=size(negspikes,1);
    negfilter=[];
    swn=0;swnst=0;shn=0;shnst=0;
    %Find first negative spike in window
    while filteri<=length_negspikes && negspikes(filteri,1)<start 
        filteri=filteri+1;
    end
    if filteri>length_negspikes
        filterdone=true;
    else
        %Find last negative spike in window
        filterj=filteri; 
        while filterj<=length_negspikes && negspikes(filterj,1)<finish
            filterj=filterj+1;
        end
        
        %Calculate metrics (mean, std) or negative spike width and height
        negfilter=negspikes(filteri:filterj-1,3);
        if sum(negfilter)==0
            swn=0;swnst=0;shn=0;shnst=0;
        else
            swn=mean(negfilter);
            swnst=std(negfilter);
            shn=mean(negspikes(filteri:filterj-1,2));
            shnst=std(negspikes(filteri:filterj-1,2));
        end
        
    end
    metrics(i).spikewidth_neg_st=swnst;     %standard deviation of neg spike width in window
    metrics(i).spikewidth_neg=swn;          %mean of neg spike width in window
    metrics(i).spikeheight_neg_st=shnst;    %standard deviation of neg spike height in window
    metrics(i).spikeheight_neg=shn;         %mean of neg spike height in window
    
    %Spike rate of positive spikes
    %Find the segment of interest in posspikes
    length_posspikes=size(posspikes,1);
    swp=0;swpst=0;shp=0;shpst=0;
    posfilter=[];
    %Find first positive spike in window
    while filteri2<=length_posspikes && posspikes(filteri2,1)<start
        filteri2=filteri2+1;
    end
    if filteri2>length_posspikes
        filter2done=true;
    else
        %Find last positive spike in window
        filterj2=filteri2; 
        while filterj2<=length_posspikes && posspikes(filterj2,1)<finish
            filterj2=filterj2+1;
        end
        
        %mean and standard deviation of spike width and height for negative spikes
        posfilter=posspikes(filteri2:filterj2-1,3);
        if sum(posfilter)==0
            swp=0;swpst=0;shp=0;shpst=0;
        else
            swp=mean(posfilter);
            swpst=std(posfilter);
            shp=mean(posspikes(filteri2:filterj2-1,2));
            shpst=std(posspikes(filteri2:filterj2-1,2));
        end
        
    end
    metrics(i).spikewidth_pos_st=swpst;     %standard deviation of pos spike width in window
    metrics(i).spikewidth_pos=swp;          %mean of pos spike width in window
    metrics(i).spikeheight_pos_st=shpst;    %standard deviation of neg spike width in window
    metrics(i).spikeheight_pos=shp;         %mean of neg spike width in window
    
    %Calculate spike rates in window
    metrics(i).spikerate_neg=length(negfilter)/Settings.movingwin(1);
    metrics(i).spikerate_pos=length(posfilter)/Settings.movingwin(1);
    metrics(i).spikerate_total=metrics(i).spikerate_pos+metrics(i).spikerate_neg;
        
    
    %-------------------
    %Slide window to next location
    start=start+Settings.movingwin(2)*fs*WindowFactor;
    finish=finish+Settings.movingwin(2)*fs*WindowFactor;
    i=i+1;
end


nrmn=mean([metrics.nrchanges]);
for i=1:length([metrics.nrchanges])
    metrics(i).nrchanges=metrics(i).nrchanges-nrmn;
end

xrange=Settings.movingwin(1)*fs*WindowFactor/2:Settings.movingwin(2)*fs*WindowFactor:Settings.movingwin(1)*fs*WindowFactor/2+(length([metrics.average])-1)*Settings.movingwin(2)*fs*WindowFactor;


%================================================
%Step 2: Compare relevant metrics to thresholds to identify noise
%================================================
% If a combination of the metrics triggers a noise-alert, this window is
% labeled as "noisy" in the vector 'noisecriteria'
% We then use the metric 'nrchanges' to find the edges of that noise event
% Finally, the vector 'noisecriteria' is used to create a matrix of noise
% events [start, end], both as indices of data.

%Specify Noise Criteria Components
%1. Number of direction changes 
noise1=[metrics.nrchanges]<1.5*ThrLower; %blue

%2. Is average data too large or too little
noise2=([metrics.average]>0.5*ThrUpper)+([metrics.average]<0.5*ThrLower); %red

%3.Is average of positive peaks too high, or average of negative peaks too low
noise3=max(([metrics.average_high]>Settings.NoiseUpper),([metrics.average_low]<Settings.NoiseLower)); %magenta

%4. Does data fluctuate too much
noise4=[metrics.st]>1*Settings.NoiseUpper; %green

%5. Are the positive data or negative peak heights fluctuating too much?
noise5=max(([metrics.spikeheight_pos_st]>ThrUpper),([metrics.spikeheight_neg_st]>-1*ThrLower)); %cyan

%Combine Noise Criteria of Interest - This is where you ultimately
%determine whether or not to flag a window as noise!
totalnoise=noise2+noise3+0*noise4+0*noise5;
noisecriteria=[totalnoise>1]; %yellow

%Fill noise array if there was noise detected
if sum(noisecriteria)>0
    
    %Find boundaries of noise event, using nrchanges
    %For each window, 'noisecriteria' will flag it as noisy (1) or not (0)
    i=1;
    while i<length(noisecriteria)
        if and(noisecriteria(i)==1,metrics(i).nrchanges<0)
            ileft=i;iright=i;
            while (metrics(ileft).nrchanges<0) && (metrics(ileft).maximum-metrics(ileft).minimum>ThrUpper-ThrLower) && (ileft>1)
                noisecriteria(ileft)=1;
                ileft=ileft-1;
            end
            while (metrics(iright).nrchanges<0) && (metrics(iright).maximum-metrics(iright).minimum>ThrUpper-ThrLower) && (iright<length(noisecriteria))
                noisecriteria(iright)=1;
                iright=iright+1;
            end
            i=max(i+1,iright);
        else
            i=i+1;
        end
    end
    
    %Use 'noisecriteria' to determine a matrix of [start, end] for each
    %noise event, as indices of 'data'
    
    %Check wether we start in noise or out of noise
    if (noisecriteria(1)==0) 
        InNoise=false; 
    else
        InNoise=true;
        noisestart=1;
    end
        
    %Continue search
    for i=1:length(noisecriteria)     
         if InNoise         %in noise event
             if noisecriteria(i)
                 %keep looking for the end of the current noise event
             else
                 %add finished event and start looking for start of next
                 noise=[noise ; [xrange(noisestart) xrange(i)]];
                 InNoise=false;
             end
         else               %not in noise event
             if noisecriteria(i)
                 %found start of next event -> start looking for the end
                 noisestart=i;
                 InNoise=true;
             else
                 %keep looking for the start of the next noise event
             end
         end
               
    end
       
    %Check wether we end in noise or out of noise
    if InNoise
        %file is over while we were in noise. Add it to the noise list
        noise=[noise ; [xrange(noisestart) xrange(i)]];
    end

    %Glue together noise events that are close together
    %Note: Noise must be separated by more than 'glue_noise' to be counted as separate
    gluednoise=[];
    i=1;
    if ~isempty(noise)
        while i<length(noise(:,1))
            glued=noise(i,:);
            while and(i<length(noise),noise(i,1)-glued(2)<Settings.glue_noise*fs); 
                glued=[glued(1) noise(i,2)];
                i=i+1;
            end
            gluednoise=[gluednoise; glued];
        end
    end
    noise=gluednoise;
end
 
%================================================
%Step 3: Optionally display window metrics and noise identification in popup window
%================================================
if (DoDisplay)
    FigureNoise=figure('units','normalized','outerposition',[0 0 1 1],'Name','Find Noise');
    subplot(2,1,1)
    plot(1:length(data),data);                                                      %plot data
    hold on
    plot(xrange,[metrics.average],'r');
    plot(xrange,[metrics.minimum],'Color',[0.5,0.5,0.5]);
    plot(xrange,[metrics.maximum],'Color',[0.5,0.5,0.5]);
    %plot(xrange,[metrics.spikewidth_pos],'k');
    %plot(xrange,[metrics.spikewidth_pos_st],'--k');
    %plot(xrange,[metrics.spikeheight_pos],'b');
    %plot(xrange,[metrics.spikeheight_pos_st],'--b');
    %plot(xrange,-1*[metrics.spikeheight_neg_st],'--r');
    %plot(xrange,-1*[metrics.spikewidth_neg_st],'r');
    
    plot([xrange(1) xrange(end)],[ThrUpper ThrUpper],'m');
    plot([xrange(1) xrange(end)],[ThrLower ThrLower],'m');
    plot([xrange(1) xrange(end)],[Settings.NoiseUpper Settings.NoiseUpper],'--r');
    plot([xrange(1) xrange(end)],[Settings.NoiseLower Settings.NoiseLower],'--r');
    plot([xrange(1) xrange(end)],[Settings.SpikesUpper Settings.SpikesUpper],'y');
    plot([xrange(1) xrange(end)],[Settings.SpikesLower Settings.SpikesLower],'y');
    minperc=num2str(round(Settings.cl*100,0));
    maxperc=num2str(round(Settings.ch*100,0));
    posperc=num2str(round(Settings.pch*100,0));
    negperc=num2str(round(Settings.ncl*100,0));
    legend('data','average','minimum','maximum',[maxperc '% data threshold (baseline)'],[minperc '% data threshold (baseline)'],[negperc '% neg spikes (noise cutoff)'],[posperc '% pos spikes (noise cutoff)'],'SD spike cutoff','SD spike cutoff');
    hold off

    
    %Second Plot is used to show only those metrics used in noise
    %detection, as well as the picked up noise for those metrics
    subplot(2,1,2)
    pldata=plot(1:length(data),data);
    axis([0,inf,-1,2]);
    hold on
    
    plot(xrange,[metrics.nrchanges],'b');           %source of noise 1
    plot(xrange,[metrics.average],'r');             %source of noise 2
    plot(xrange,[metrics.average_high],'m');        %source of noise 3
    plot(xrange,[metrics.average_low],'m');         %source of noise 3
    plot(xrange,[metrics.st],'g');                  %source of noise 4
    plot(xrange,[metrics.spikeheight_pos_st],'c');  %source of noise 5
    plot(xrange,[metrics.spikeheight_neg_st],'c');  %source of noise 5
       
    pln1=plot(xrange,1.5*(3*noise1-2),'.b');                   %noise 1 triggered
    pln2=plot(xrange,1.55*(3*noise2-2),'.r');                  %noise 2 triggered
    pln3=plot(xrange,1.6*(3*noise3-2),'.m');                   %noise 3 triggered
    pln4=plot(xrange,1.65*(3*noise4-2),'.g');                  %noise 4 triggered
    pln5=plot(xrange,1.7*(3*noise5-2),'.c');                   %noise 5 triggered
    plnn=plot(xrange,1.75*(3*noisecriteria-2),'.y','Linewidth',1); %NOISE triggered
    
    if ~isempty(noise)
       for i=1:length(noise(:,1))
           plot([noise(i,1) noise(i,2)],[1.8 1.8],'y','LineWidth',2);
       end
       %lastnoise=plot([noise(i,1) noise(i,2)],[1.8 1.8],'y','LineWidth',2);
    end
    
    plot([xrange(1) xrange(end)],0.5*[ThrUpper ThrUpper],'--r');
    plot([xrange(1) xrange(end)],0.5*[ThrLower ThrLower],'--r');
    minperc=num2str(round(Settings.cl*100,0));
    maxperc=num2str(round(Settings.ch*100,0));
    posperc=num2str(round(Settings.pch*100,0));
    negperc=num2str(round(Settings.ncl*100,0));

    legend([pldata,pln1,pln2,pln3,pln4,pln5,plnn],{'Data','Noise 1 triggered','Noise 2 triggered', 'Noise 3 triggered', 'Noise 4 triggered','Noise 5 triggered','Identified Noise'});
    
    %UNUSED METRIC PLOTS
    %plot(xrange,[metrics.average_high],'r');
    %plot(xrange,[metrics.average_low],'c');
    %plot(xrange,[metrics.highpeakperc],'b');
    %plot(xrange,[metrics.average],'m');
    %plot(xrange,[metrics.minimum],'k');
    %plot(xrange,[metrics.maximum],'y');
    %plot([xrange(1) xrange(end)],[ThrUpper ThrUpper],'--b');
    %plot([xrange(1) xrange(end)],[ThrLower ThrLower],'--b');
    %plot([xrange(1) xrange(end)],[Settings.NoiseUpper Settings.NoiseUpper],'--m');
    %plot([xrange(1) xrange(end)],[Settings.NoiseLower Settings.NoiseLower],'--m');
    %plot(xrange,0.1*[metrics.tracelength],'c');
    %plot(xrange,[metrics.spikewidth_pos],'c');
    %plot(xrange,[metrics.spikewidth_neg],'m');
    %plot(xrange,[metrics.spikewidth_total],'--k');
    %plot(xrange,[metrics.increaserate],'--y');
    %plot(xrange,[metrics.spikerate_pos],'--k');
    %plot(xrange,[metrics.spikerate_neg],'--y');
    %plot(xrange,[metrics.spikewidth_total],'--k');
    %plot(xrange,[metrics.increaserate],'--y');
    %plot(xrange,[metrics.absmean],'--y');
    %plot(xrange,[metrics.zerox],'r');
    %plot(xrange,[metrics.meanx],'--r');
    %plot(xrange,[metrics.st_high],'--k');
    %plot(xrange,[metrics.spikewidth_pos],'k');
    %plot(xrange,[metrics.highpeakperc],'b');

    hold off
end
