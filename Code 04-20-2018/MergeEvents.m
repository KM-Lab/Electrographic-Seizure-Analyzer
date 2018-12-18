%%MergeEvents
function [events,noise]=MergeEvents(posevents,negevents,noise,HowToMerge,DealWithNoise,EventCutoff,NoiseCutoff,EventGlue)
%Combines lists of positive, and negative events, as well as noise into one using either intersection, union, or another method
%
%Input: 'posevents'    matrix of all positive events [start,end] as indices of 'data'
%       'negevents'    matrix of all negative events [start,end] as indices of 'data'
%       'noise'        matrix of all noise [start, end] as indices of 'data'
%
%       'HowToMerge'    Specify how to merge positive and negative events
%                       0 - Union events      (events get longer)
%                       1 - Intersect events  (only overlapping sections count)
%                       2 - Only take positive events
%                       3 - Only take negative events
%
%       'DealWithNoise' Specify how to reduce events, based on noisy signal
%                       0 - Ignore all noise 
%                       1 - Split events up at noise locations, ignoring short noise (specified in 'NoiseCutoff')
%                       2 - Ignore events that contain noise, eliminating the entire events
%
%       'NoiseCutoff'   Ignore all noise events less long than this.
%       'EventCutoff'   If the merging process results in events less long than 'EventCutoff', remove them
%       'EventGlue'     Glue events that are closer to each other  than this (in s)
%
%Output: 'events'       Matrix of event locations, [start,end] as indices of 'data'
%        'noise'        Matrix of noise locations, [start,end] as indices of 'data'
%
%Calls:         none
%
%Called by:     'MenuAnalyze_Callback' in 'findallofthem.m'

%First merge all positive and negative events
switch HowToMerge
    case 3
        %only consider negative events
        events=negevents;
    case 2
        %only consider positive events
        events=posevents;
    case 1
        %intersect positive and negative events
        number_of_positives=size(posevents,1);
        number_of_negatives=size(negevents,1);
        
        %Initialize
        i_p=1; %index of positive event
        i_n=1; %index of negative event
        events=[];
        
        while and(i_p<=number_of_positives,i_n<=number_of_negatives)
            %Search for overlap, keep moving events until found
            foundoverlap=1;
            if posevents(i_p,2)<negevents(i_n,1)
                i_p=i_p+1;
                foundoverlap=0;
            elseif negevents(i_n,2)<posevents(i_p,1)
                i_n=i_n+1;
                foundoverlap=0;
            end
            
            %This part is only reached if there is overlap
            if foundoverlap
                %Add event, if still long enough
                events=[events; [max(posevents(i_p,1),negevents(i_n,1)),min(posevents(i_p,2),negevents(i_n,2))]];
                
                if posevents(i_p,2)<negevents(i_n,2)
                    i_p=i_p+1;
                else
                    i_n=i_n+1;
                end
            end
        end
        %After positives or negatives run out, there is no more overlap
        
    case 0
        %union positive and negative events
        number_of_positives=size(posevents,1);
        number_of_negatives=size(negevents,1);
        
        %Initialize
        i_p=1; %index of positive event
        i_n=1; %index of negative event
        events=[];
        InEvent=0; %0 is no event, 1 is in event
        
        %Start looping through the events 
        while and(i_p<=number_of_positives,i_n<=number_of_negatives)
            if ~InEvent
                %Find the start of the next event
                if posevents(i_p,1)<negevents(i_n,1)
                    currentevent='positive';
                    startevent=posevents(i_p,1);
                else
                    currentevent='negative';
                    startevent=negevents(i_n,1);
                end
                InEvent=1;
            else
                %See if the event continues
                if strcmp('positive',currentevent);
                    %check if the negative event starts before the current 
                    %positive event ends.
                    if  negevents(i_n,1)<=posevents(i_p,2)
                        %event continues
                        if negevents(i_n,2)<=posevents(i_p,2)
                            currentevent='positive';
                            i_n=i_n+1;
                        else
                            currentevent='negative';
                        end
                    else
                        %event is over, add event to union list
                        events=[events; [startevent posevents(i_p,2)]];
                        InEvent=0;
                        i_p=i_p+1;
                    end
                    
                else
                    %check if the next positive event starts before the
                    %current negative event ends.
                    if  posevents(i_p,1)<=negevents(i_n,2)
                        %event continues
                        if posevents(i_p,2)<=negevents(i_n,2)
                            currentevent='negative';
                            i_p=i_p+1;
                        else
                            currentevent='positive';
                        end
                    else
                        %event is over, add event to union list
                        events=[events; [startevent negevents(i_n,2)]];
                        InEvent=0;
                        i_n=i_n+1;
                    end
                    
                end
            end
        end
        %We ran out of either positive or negative events
        %Add the event we were in the middle off
        if InEvent
            if strcmp(currentevent,'positive')
                events=[events; [startevent posevents(i_p,2)]];
                i_p=i_p+1;
            else
                events=[events; [startevent negevents(i_n,2)]];
                i_n=i_n+1;
            end
        end
        
        %At most one of the next two while loops will run
        while i_p<=number_of_positives
            events=[events; [posevents(i_p,1) posevents(i_p,2)]];
            i_p=i_p+1;
        end
        while i_n<=number_of_negatives
            events=[events; [negevents(i_n,1) negevents(i_n,2)]];
            i_n=i_n+1;
        end
            
        
end



%if noise is empty, don't attempt to deal with it
if or(isempty(noise),isempty(events))
    DealWithNoise=0;
    afternoise=events;
end

%Deal with noise
switch DealWithNoise
   
    case 0
        %Ignore noise completely
        afternoise=events;
        
    case 1
        %Delete all noise from events, cutting events into smaller events
        %Delete those events if they become too short.
        %Also, ignore noise if it is shorter than given amount
        afternoise=[];  %Eventlist after noise
        newnoise=[];    %Noise list after ignoring short noise
        i_events=1;
        i_noise=1;
        eventstart=0;
        if size(events,1)>0
            eventstart=events(i_events,1);
        end
        while and(i_events<size(events,1),i_noise<=size(noise,1))
            if noise(i_noise,2)-noise(i_noise,1)<NoiseCutoff
                %Noise that is too short gets ignored
                i_noise=i_noise+1;
            else
                if (noise(i_noise,2)<eventstart)
                    %If current noise ends before current event starts, move to next noise
                    i_noise=i_noise+1;
                elseif (events(i_events,2)<noise(i_noise,1))
                    %If current event ends before current noise starts, add it unchanged
                    %Then move to the next event
                    afternoise=[afternoise; [eventstart events(i_events,2)]];
                    i_events=i_events+1;
                    if size(events,1)>=i_events
                        eventstart=events(i_events,1);
                    end
                else
                    if noise(i_noise,2)>=events(i_events,2)
                        %noise lasts longer than event
                        if noise(i_noise,1)>eventstart
                            %event is cut short by noise
                            cutevent=[eventstart noise(i_noise,1)];
                            afternoise=[afternoise; cutevent];    
                        else
                            %completely ignore the event, since it is smaller
                            %than the noise
                        end
                        i_events=i_events+1;
                        if size(events,1)>=i_events
                            eventstart=events(i_events,1);
                        end
                    else
                        %noise ends before event does. There may be another
                        %noise starting
                        if noise(i_noise,1)>eventstart
                            %event is cut short by noise
                            cutevent=[eventstart noise(i_noise,1)];
                            afternoise=[afternoise; cutevent];    
                        end
                        %still in the same event, but with different start
                        %point.
                        eventstart=noise(i_noise,2);
                        i_noise=i_noise+1;
                    end
                end
            end
        end
        while i_events<=size(events,1)
            %Add any events that happen after the last noise
            afternoise=[afternoise; [eventstart events(i_events,2)]];
            i_events=i_events+1;
            if i_events<=size(events,1)
                eventstart=events(i_events,1);
            end
        end
        
    case 2
        %Ignore all events that contain any noise, exceeding the NoiseCutoff
        i_noise=1;
        i_events=1;
        afternoise=[];
        newnoise=[];
        while and(i_events<size(events,1),i_noise<size(noise,1))
            if events(i_events,1)>noise(i_noise,2)
                %Current event starts after current noise. Move noise up
                i_noise=i_noise+1;
            elseif  events(i_events,2)<noise(i_noise,1)
                %Current event ends before current noise. Move event up
                afternoise=[afternoise; events(i_events,:)];
                i_events=i_events+1;
            else
                %There is Noise-Event overlap.
                if noise(i_noise,2)-noise(i_noise,1)>NoiseCutoff
                    %This noise event is too long -> Ignore the event
                    i_events=i_events+1;
                else
                    %This is a short noise event and does not warrant
                    %removing event. Consider the next noise
                    i_noise=i_noise+1;
                end
            end
        end
        
        while i_events<=size(events,1)
            %Add any events that happen after the last noise
            afternoise=[afternoise; events(i_events,:)];
            i_events=i_events+1;
        end
        
        %return events
end

%Remove the snort noise, that doesn't meet the cutoffs
oldnoise=noise;
noise=[];
if ~isempty(oldnoise)
    for i=1:length(oldnoise(:,1))
        if oldnoise(i,2)-oldnoise(i,1)>=NoiseCutoff
            noise=[noise; oldnoise(i,:)];
        end
    end
end

%Remove the short events, that don't meet the cutoffs
%This might have occured by chopping up events at noise locations
events=[];
if ~isempty(afternoise)
  for i=1:length(afternoise(:,1))
    if afternoise(i,2)-afternoise(i,1)>=EventCutoff
          events=[events; afternoise(i,:)];
    end
  end
end

%Join events that are close enough to each other (<EventGlue (in s))
%Make sure to stop them if noise is detected
joinedevents=[];
if ~isempty(events) %only run if events is nonempty
    nrevents=length(events(:,1));
    if nrevents>1 %only run if there is more than 1 event
        istart=1;icurrent=2;
        while icurrent<=length(events(:,1))
            if or(events(icurrent,1)-events(icurrent-1,2)>EventGlue,... %gap too long
                    (~isempty(noise) && sum((noise(:,1)>=events(icurrent-1,2)).*(noise(:,2)<=events(icurrent,1)))>0)) %noise in gap
                joinedevents=[joinedevents; events(istart,1) events(icurrent-1,2)];
                istart=icurrent;
                icurrent=istart+1;
            else
                icurrent=icurrent+1;
            end
        end
        joinedevents=[joinedevents; events(istart,1) events(icurrent-1,2)];
    else
        joinedevents=events;
    end
    
    events=joinedevents;
end



