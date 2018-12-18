function Settings=DefaultSettings()
%Load in all the default settings for the automated temporal analysis
%This will create a Settings data structure, which is fed into the various
%functions, who will take the components they need.
% Input:        none
%
% Output:       'Settings'      Datastructure containing all settings
%                               More components are added throughout code
% Calls:  
%
% Called by:    'set_default_Callback' in 'findallofthem.m'
%               'findallofthem_OpeningFcn' in 'findallofthem.m'

%Create data structure
Settings=struct;

%User Interface & File Management
Settings.buffer = 5;             %Allow program to load upto # files simultaneously. These are stored in buffer for quick toggling between files.
Settings.eventpadding=1;         %When saving an event, save 1s before and after event for context (not currently using)
Settings.outputfolder='Output';  %Subfolder for all output events


%SETTINGS FOR PREANALYSIS (HISTOGRAM ANALYSIS OF PEAKS)
%Peekseek Settings for Preanalysis of ALL spikes. The cutoffs of 95% of data are
%pretty helpful in detecting all spikes. This is why the thresholds are set so small. 
Settings.all_dist_pos = 50;       %min distance between pos spikes
Settings.all_dist_neg = 50;       %min distance between neg spikes
Settings.all_thrsh_pos = 0.001;   %amplitude threshold for positive spikes 
Settings.all_thrsh_neg = 0.001;   %amplitude threshold for negative spikes
Settings.number_of_sd = 2;        %number of standard deviations away from mean to consider spike


%Threshold Settings for Amplitude Analysis
Settings.cl=0.05;   %low cutoff for all spikes
Settings.ch=0.95;   %high cutoff for all spikes
Settings.pcl=0.05;  %low cutoff for positive spikes (in preanalysis)
Settings.pch=0.95;  %high cutoff for positive spikes (in preanalysis)
Settings.ncl=0.05;  %low cutoff for negative spikes (in preanalysis)
Settings.nch=0.95;  %high cutoff for negative spikes (in preanalysis)


%SETTINGS FOR FINDING SPIKES (peekseak in main loop)
%Peaks are only considered to be 'spikes' if they exceed a multiplier of the threshold, are at least dist apart (in samples), 
%and the width is measured at the height prescribed, as a percentage of peak height. 
%Let's say 95% of ALL negative spikes lie below -0.6, then given a
%multiplier of 0.5 (50%) we only consider those spikes that exceed -0.3 (50% of -0.6)
Settings.dist_pos = 50;      %positive distance between spikes
Settings.dist_neg = 50;      %negative distance between spikes

Settings.thrsh_pos = 0.50;      %threshold for positive spikes 
Settings.thrsh_neg = 0.50;      %threshold for negative spikes 
   
Settings.ht_perc = 0.5;      %determines percent of height for each spike needed to det. width
Settings.req_width = 0.05;   %req width threshold for each spike (in seconds)

   
    
%Settings for finding electrographic seizures
Settings.artif_rmvl  = 0.1;  %threshold to consider a spike artificial or not
Settings.artif_limit = 0.8;  %evnt needs to have less than limit to not be considered artificial event
                             %not currently used

%SETTINGS FOR FINDING EVENTS
%An event is defined as a series of 'spikes' that:
%  * are dense enough  (maintains 'min_nr_spikes' / 'min_szre_windw')
%  * lasts long enough ('min_szre_lngth')
%If two events are closer than 'eventglue' (s), make them the same event
%If two events are only separated by few wide spikes that don't indicate noise, optionally merge them
Settings.min_szre_lngth = 2; %min length req for seiz evnt 
Settings.min_szre_windw = 2; %min length req for rolling time window 
Settings.min_nr_spikes = 4;  %minimum number of spikes to count as event
Settings.connect_events = 1; %attempt to connect events if separated by few dense wide spikes
                             %The wide spikes would have to adhere to spike conditions
                             %0 = don't merge
                             %1 = merge
Settings.eventglue=2;        %glue events together if they are less than #s apart


%SETTINGS FOR IDENTIFYING NOISE
%NOTE: More noise settings are automatically added after threshold analysis
Settings.noisecutoff=0.5*Settings.min_szre_windw; %delete noise shorter than this (s)
Settings.glue_noise=1;                            %glue noise events together if they are less than #s apart
Settings.movingwin=[1, 0.5];                      %Sliding window [size,overlap] in s
                                                  %Resolution of noise detection determined by 'size'



%SETTINGS FOR MERGING POSITIVE EVENTS, NEGATIVE EVENTS, NOISE
Settings.deal_with_noise=1; %0=ignore all noise, 1=truncate events at noise, 2=ignore noisy events altogether
Settings.deal_with_merge=1; %0=union, 1=intersect, 2=positive only, 3=negative only
end

