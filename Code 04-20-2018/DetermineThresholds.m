function [chline,clline,pchline,nclline,sdposline,sdnegline,sddataline]=DetermineThresholds(data,Settings,DoDisplay)
%Given data, shifted to recenter around 0, determine cutoffs of positive data points, negative datapoints, and all datapoints. 
%Used to set automatic thresholds in spike-determination. Optionally display the results in popout     
%Input:         'data'       contains data points of 1 file, as a row vector, recentered around 0
%               'Settings'   contains all settings parameters
%               'DoDisplay'  0=don't show, 1=show
%    
%Output:        'chline'     cutoff high of settings.ch% of all data
%               'clline'     cutoff low of settings.cl% of all data
%               'pchline'    cutoff high of settings.pch% of all data
%               'nclline'    cutoff low of settings.ncl% of all data
%
%Calls:         DetermineCutoffs
%
%Called by:     'CutoffsAutomic_Callback' in 'findallofthem.m'
%
    %determine ALL spikes, including the teeny ones
    spp = peakseek(data, Settings.all_dist_pos,Settings.all_thrsh_pos);   %positive peak loc & data
    spn = peakseek(-data, Settings.all_dist_neg,Settings.all_thrsh_neg);  %negative peak loc & data

    %determine standard deviations of data, and peaks
    sd_multiplier=Settings.number_of_sd;
    sd_data=std(data); sddataline=sd_multiplier*sd_data;
    sd_pos=std(data(spp)); sdposline=sd_multiplier*sd_pos;
    sd_neg=std(data(spn)); sdnegline=-1*sd_multiplier*sd_neg;
    
    
    %histogram count, and cumulative histogram. Also find indices for low and high cutoffs
    [N,Edges, cutoff_l,cutoff_h,clline,chline]=DetermineCutoffs(data,Settings.cl,Settings.ch);  

    %Find the low and high cutoff of the positive spikes
    pchline=0;cnlline=0;pN=0;nN=0; %Make sure program doesn't crash if spp,spn are empty
    if length(spp)>0
        ppeakdata=data(spp(1,:));
        [pN,pEdges, pcutoff_l,pcutoff_h,pclline,pchline]=DetermineCutoffs(ppeakdata,Settings.pcl,Settings.pch);
    end
    
    %Find the low and high cutoff of the negative spikes
    if length(spn)>0
        npeakdata=data(spn(1,:));
        [nN,nEdges, ncutoff_l,ncutoff_h,nclline,nchline]=DetermineCutoffs(npeakdata,Settings.ncl,Settings.nch);
    end

    %Display histograms on screen
    if DoDisplay
        FigureThresholds=figure('units','normalized','outerposition',[0 0 1 1],'Name','Determine Thresholds');
        %plot data
        subplot(3,3,[1 2 3]);
        plot(data);
        hold on
        plotcl=plot([0 length(data)-1],[clline clline],'m');                %baseline (% of ALL data)
        plotch=plot([0 length(data)-1],[chline chline],'m');                %baseline (% of ALL data)
        plotpch=plot([0 length(data)-1],[pchline pchline],'--r');           %noise cutoff (%of +data)
        plotncl=plot([0 length(data)-1],[nclline nclline],'--r');           %noise cutoff (%of +data)
        plotsddata=plot([0 length(data)-1],[-sddataline -sddataline],'b');  
        plotsdneg=plot([0 length(data)-1],[sdnegline sdnegline],'y');       %spikes (#SD of -data)
        plotsdpos=plot([0 length(data)-1],[sdposline sdposline],'y');       %spikes (#SD of +data
        plotsddata2=plot([0 length(data)-1],[sddataline sddataline],'b');
        hold off
        
        %plot histograms
        subplot(3,3,4);
        histogram(data);
        title('All datapoints');
        axis([2*clline 2*chline 0 0.5*max(N)]);
        hold on
        plot([clline clline],[0 0.5*max(N)],'m');
        plot([chline chline],[0 0.5*max(N)],'m');
        plot([-sddataline -sddataline],[0 0.5*max(N)],'b');
        plot([sddataline sddataline],[0 0.5*max(N)],'b');
        plot([sdnegline sdnegline],[0 max(nN)],'y');
        plot([sdposline sdposline],[0 max(pN)],'y');
        hold off
        
        subplot(3,3,5);
        histogram(data(spp(1,:)));
        title('Positive data');
        axis([0 2*pchline 0 inf]);
        hold on
        plot([pchline pchline],[0 max(pN)],'--r');
        plot([sdposline sdposline],[0 max(pN)],'y');
        hold off
        
        subplot(3,3,6);
        histogram(data(spn(1,:)));
        title('Negative data');
        axis([2*nclline 0 0 inf]);
        hold on
        plot([nclline nclline],[0 max(nN)],'--r');
        plot([sdnegline sdnegline],[0 max(nN)],'y');
        hold off

        %Add legend
        subplr=subplot(3,3,[7,8,9]);
        axis off
        legendbaseline=[num2str(100*(Settings.ch-Settings.cl)) '% of data (detect baseline)'];
        noisebaseline=[num2str(100*Settings.pch) '% of + or - data (used for noise)'];
        noisebaseline2=[num2str(100*(1-Settings.ncl)) '%'];
        legend(subplr,[plotcl,plotsddata,plotsdneg,plotsdpos,plotpch],{legendbaseline,[num2str(sd_multiplier) 'SD data'],[num2str(sd_multiplier) 'SD -Data (detect -spikes)'],[num2str(sd_multiplier) 'SD +Data (detect +spikes)'],noisebaseline})
        
    end
end


function [N,Edges,cutoff_low,cutoff_high,clline,chline]=DetermineCutoffs(data,cl,ch)
%Given a dataset, and a lower and upper threshold for data, return the histcount,
%the indices for the cutoff values specified and the values reached here.
%
%Input:         'data'          Row vector of numerical data
%               'cl'            Proportion of data - low threshold 
%               'ch'            Proportion of data - high threshold
%
%Output:        'N'             Number of data points in each bin
%               'Edges'         Edges of each bin
%               'cutoff_low'    Number of data points below 'cl'%
%               'cutoff_high'   Number of data points above 'ch'%
%               'clline'        Value of data at 'cl'%
%               'chline'        Value of data at 'ch'%
%
%Calls:         
%
%Called by:     DetermineThresholds
%       
    [N,Edges]=histcounts(data);
    cumN=zeros(1,length(N));
    cutoff_low=0; clline=0;
    cutoff_high=0; chline=0;
    for i=1:length(N)
        cumN(i)=sum(N(1:i))/length(data);
        if cutoff_low==0 && cumN(i)>cl
            cutoff_low=i;
            clline=Edges(cutoff_low);
        end
        if cutoff_high==0 && cumN(i)>ch
            cutoff_high=i;
            chline=Edges(cutoff_high);
        end
    end
end

