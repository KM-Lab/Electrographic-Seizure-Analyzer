function varargout = findallofthem(varargin)
% FINDALLOFTHEM MATLAB code for findallofthem.fig
%      FINDALLOFTHEM, by itself, creates a new FINDALLOFTHEM or raises the existing
%      singleton*.
%
%      H = FINDALLOFTHEM returns the handle to a new FINDALLOFTHEM or the handle to
%      the existing singleton*.
%
%      FINDALLOFTHEM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FINDALLOFTHEM.M with the given input arguments.
%
%      FINDALLOFTHEM('Property','Value',...) creates a new FINDALLOFTHEM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before findallofthem_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to findallofthem_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help findallofthem

% Last Modified by GUIDE v2.5 03-Nov-2017 10:09:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @findallofthem_OpeningFcn, ...
                   'gui_OutputFcn',  @findallofthem_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before findallofthem is made visible.
function findallofthem_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to findallofthem (see VARARGIN)

% Choose default command line output for findallofthem
handles.output = hObject;

%Change visibility of various panels
set(handles.HidePanelLeft,'Visible','off'); %Don't hide left panel
set(handles.HidePanelRight,'Visible','on'); %Hide right panel
set(handles.ChangeStep2,'Visible','off');
set(handles.MenuShowIt,'Value',1);

%Load in default settings
AllInfo = struct;
AllInfo.Settings=DefaultSettings();
set(handles.MenuFiles,'UserData',AllInfo);

% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = findallofthem_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%% -----------------------------------------------------------------
%%% Step 1: Load in Files
%%% -----------------------------------------------------------------

% --- Executes on button press in MenuPreanalysis.
function MenuPreanalysis_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPreanalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Check if a directory and channel have been selected
CurrentDirectory = get(handles.MenuCurrentDir,'String');
if strcmp(CurrentDirectory,'No Directory Selected')
   %Show warning on screen
   set(handles.MenuWarning,'String','Please select a directory');
   return
end

set(handles.MenuWarning,'String','');
idx= get(handles.MenuChannel,'Value');   % get the index of the selected channel
CurrentChannel = idx-1;
if idx<2
    %Show warning on screen
    set(handles.MenuWarning,'String','Please select a channel');
    return 
end

AllInfo=get(handles.MenuFiles,'UserData');  %obtain the structure AllInfo
outputfolder=AllInfo.Settings.outputfolder  %obtain the ouput folder

%Analyze chosen directory; determine number of files
DoCheck=get(handles.MenuOverwriteMeta,'Value');
AllFiles=AnalyzeDirectory(CurrentDirectory,CurrentChannel,DoCheck,outputfolder);

%Show warning on screen if there are no files for the channel
if ~iscell(AllFiles)
    set(handles.MenuWarning,'String',['No files available for channel ' num2str(CurrentChannel)]);
    return
end

NumberOfFiles=length(AllFiles);
AllInfo=get(handles.MenuFiles,'UserData');  %obtain the structure AllInfo
Settings=AllInfo.Settings;                  %unpack the settings from AllInfo

if NumberOfFiles<1
    %show warning on screen, since directory is empty
    set(handles.MenuWarning,'String','Directory does not contain *.mat files');
    return
else
    %Create data structure with information on the directory
    AllInfo.CurrentDirectory=CurrentDirectory;
    AllInfo.AllFiles=AllFiles;
    AllInfo.NumberOfFiles=NumberOfFiles;
    AllInfo.ActiveFile=1;
    AllInfo.CurrentChannel=CurrentChannel;
    AllInfo.BufferSize=Settings.buffer;
    AllInfo.NextBuffer=1;
    
    %Initialize data structure to store information about each individual file
    DoLoadAnalysis=get(handles.MenuLoadAnalysis,'Value');
    Didit=0;    %We did not load in previous analysis yet
    
    if DoLoadAnalysis
        %Load analysis if it exists
        set(handles.MenuFiles,'UserData',AllInfo);
        Didit=LoadAnalysis(handles);
        if Didit
            set(handles.MenuWarningLoading,'String','Loaded Analysis Information From File');
        else
            set(handles.MenuWarningLoading,'String','Could Not Load Analysis From File');
        end
        AllInfo=get(handles.MenuFiles,'UserData');
    else
        set(handles.MenuWarningLoading,'String','Loaded Data (w/o possible previous analysis)');
    end
    if or(~DoLoadAnalysis,~Didit)
        %Create a new data-structure (did not load data, or loading failed)
        AllInfo.Settings.DoShowAnalysis=false;  
        AllInfo.FileInfo=struct;
        for i=1:NumberOfFiles
            AllInfo.FileInfo(i).Loaded=0;       %Mark that no fileinfo has been loaded yet for each file
            AllInfo.FileInfo(i).DataBuffered=0; %Mark the location in the buffer of this data (0 is not buffered)
        end
        
        
        %Initialize buffer structure
        %Buffer.File:   refers to filenumber stored in this buffer
        %Buffer.Data:   contains the data
        AllInfo.Buffer=struct;
        for i=1:Settings.buffer
            AllInfo.Buffer(i).File=0;
        end
        
        %Empty the settings regarding cutoffs
        chline=0;
        clline=0;
        pchline=0;
        nclline=0;
        AllInfo.Settings.BaselineUpper=chline;
        AllInfo.Settings.BaselineLower=clline;
        AllInfo.Settings.SpikesUpper=pchline;
        AllInfo.Settings.SpikesLower=nclline;
        set(handles.BaselineUpper,'String',num2str(chline));
        set(handles.BaselineLower,'String',num2str(clline));
        set(handles.SpikesUpper,'String',num2str(pchline));
        set(handles.SpikesLower,'String',num2str(nclline));
    end
    
    %Store this data in the MenuFiles part of the GUI
    set(handles.MenuFiles,'UserData',AllInfo);
    
    %Hide inputpanel, by overlaying summary panel
    set(handles.PanelFileName,'String',CurrentDirectory);
    set(handles.PanelChannel,'String',['Channel: ' num2str(CurrentChannel)]);
    set(handles.HidePanelLeft,'Visible','on');
    set(handles.HidePanelMiddle,'Visible','off');
    set(handles.ChangeStep2,'Visible','on');
    
    %Initialize screen for first file in directory
    DisplayFile(AllInfo.ActiveFile,handles);
end

% --- Executes on button press in MenuSelectDir.
function MenuSelectDir_Callback(hObject, eventdata, handles)
% hObject    handle to MenuSelectDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder_name = uigetdir();
set(handles.MenuCurrentDir, 'String', folder_name);


function data=LoadData(FileNumber,AddToBuffer,handles)
%Return data associated with given file number.
%If the data is in the buffer, obtain it from there.
%If the data is not in the buffer, then add it to the buffer if and only if
%AddToBuffer=1.
Info=get(handles.MenuFiles,'UserData');

if (Info.FileInfo(FileNumber).DataBuffered==0)
    %File is not in the buffer
    file_name = char(fullfile(Info.CurrentDirectory, Info.AllFiles(FileNumber)));  
    imported  = load(file_name); 
    data = imported.sbuf(1:end,Info.CurrentChannel);   %imports data values
    data=data-mean(data);                              %recenter around 0
    
    if (AddToBuffer)
        %Add data to buffer and update buffer reference
        if (Info.Buffer(Info.NextBuffer).File~=0)
            Info.FileInfo(Info.Buffer(Info.NextBuffer).File).DataBuffered=0; %Remove reference to buffer location from file
        end 
        Info.FileInfo(FileNumber).DataBuffered=Info.NextBuffer;          %Add reference to buffer location for new file
        Info.Buffer(Info.NextBuffer).File=FileNumber;                    %Buffer references current file number
        Info.Buffer(Info.NextBuffer).Data=data;                          %Load data into buffer
        Info.NextBuffer=1+mod(Info.NextBuffer,Info.BufferSize);          %Update the next buffer to change
        set(handles.MenuFiles,'UserData',Info);
    end
else
    %File is in the buffer already. Simply return that data
    data=Info.Buffer(Info.FileInfo(FileNumber).DataBuffered).Data;    
end

%If data has never been loaded, also load in other parameters
if not(Info.FileInfo(FileNumber).Loaded)
    Info.FileInfo(FileNumber).Fs=imported.fs; %load in the sampling frequency
    Info.Fs=imported.fs; %All files in the directory should share the same sampling frequency.  %%%BUILD IN CHECK%%%
    Info.FileInfo(FileNumber).upperrange=(length(data)-1)/imported.fs;
    Info.FileInfo(FileNumber).Loaded=1; %The file has now been loaded
    timestamps=imported.trdata.timestamp;
    Info.FileInfo(FileNumber).StartTime=timestamps(1);  %Start time of the file
    Info.FileInfo(FileNumber).EndTime=timestamps(end);  %End time of the file
    
    set(handles.MenuFiles,'UserData',Info);
end



%%% -----------------------------------------------------------------
%%% Step 2: Change Automatic Thresholds Panel
%%% -----------------------------------------------------------------

% --- Executes on button press in MenuCutoffs.
% User accepted the preanalysis cutoff values (automatic or manual)
function MenuCutoffs_Callback(hObject, eventdata, handles)
% hObject    handle to MenuCutoffs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Overlay the hiding panel, and fill it with information
set(handles.HidePanelMiddle,'Visible','on');
set(handles.ChangeStep2,'Visible','on');
set(handles.HidePanelRight,'Visible','off');
set(handles.CutoffChoice,'String',get(handles.CutoffsMethod,'String'));
set(handles.b1,'String',get(handles.BaselineUpper,'String'));
set(handles.b2,'String',get(handles.BaselineLower,'String'));
set(handles.s1,'String',get(handles.SpikesUpper,'String'));
set(handles.s2,'String',get(handles.SpikesLower,'String'));


% --- Apply manual changes to the thresholds
function CutoffsApplyChanges_Callback(hObject, eventdata, handles)
% hObject    handle to CutoffsApplyChanges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Info=get(handles.MenuFiles,'UserData');
Settings=Info.Settings;
set(handles.CutoffsMethod,'String','Current Values: Manual');

%Change UserData: load selected values into settings
Settings.BaselineUpper=str2num(get(handles.BaselineUpper,'String'));
Settings.BaselineLower=str2num(get(handles.BaselineLower,'String'));
Settings.SpikesUpper=str2num(get(handles.SpikesUpper,'String'));
Settings.SpikesLower=str2num(get(handles.SpikesLower,'String'));
Settings.CutoffsMethod=get(handles.CutoffsMethod,'String');

%Reset the analysis, because with changed values it is no longer valid
%%% ???this code appears twice
    for i=1:Info.NumberOfFiles
        Info.FileInfo(i).SpikesPositive=[];
        Info.FileInfo(i).SpikesNegative=[];
        Info.FileInfo(i).psn=[];
        Info.FileInfo(i).nsn=[];
        Info.FileInfo(i).psw=[];
        Info.FileInfo(i).nsw=[];
        Info.FileInfo(i).Noise=[];
        Info.FileInfo(i).EventsPos=[];
        Info.FileInfo(i).EventsNeg=[];
        Info.FileInfo(i).EventsCombined=[];
    end

Info.Settings=Settings;
set(handles.MenuFiles,'UserData',Info);

RefreshDisplay(handles);


% --- Determine the thresholds automatically, based on current settings parameters
function CutoffsAutomatic_Callback(hObject, eventdata, handles)
% hObject    handle to CutoffsAutomatic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Determine the cutoffs automatically based on current file.
%If settings haven't changed, just load them in from earlier this session.
Info=get(handles.MenuFiles,'UserData');
Settings=Info.Settings;
data=Info.Buffer(Info.FileInfo(Info.ActiveFile).DataBuffered).Data;

if not(isfield(Info.FileInfo(Info.ActiveFile),'AutoCutoffs'))
    Info.FileInfo(Info.ActiveFile).AutoCutoffs=[];
end

SettingsCutoffsChanged=true; %%%currently always rerun, but could optimize here later
MustCalculate=or(SettingsCutoffsChanged,...                                     %Setting changed
                 isempty(Info.FileInfo(Info.ActiveFile).AutoCutoffs));          %never calculated before
             
if MustCalculate %always runs
    ShowThresholdPopup=get(handles.PopupThreshold,'Value');
    %Calculate automatic threshold cutoffs based on histogram information
    %of positive and negative data peaks
    [chline, clline,pchline,nclline,sdposline,sdnegline,sddataline]=DetermineThresholds(data,Settings,ShowThresholdPopup);
    Info.FileInfo(Info.ActiveFile).AutoCutoffs=[chline,clline,pchline,nclline,sdposline,sdnegline,sddataline];
else %currently never runs
    cutoffs=Info.FileInfo(Info.ActiveFile).AutoCutoffs;
    chline=cutoffs(1);
    clline=cutoffs(2);
    pchline=cutoffs(3);
    nclline=cutoffs(4);
    sdposline=cutoffs(5);
    sdnegline=cutoffs(6);
    sddataline=cutoffs(7);
end

%Update GUI
set(handles.BaselineUpper,'String',num2str(chline));
set(handles.BaselineLower,'String',num2str(clline));
set(handles.SpikesUpper,'String',num2str(sdposline));
set(handles.SpikesLower,'String',num2str(sdnegline));
set(handles.CutoffsMethod,'String',['Current Values: Automatic - File ' num2str(Info.ActiveFile)]);    

%Change UserData: load selected values into settings
CutoffsSettingsChanged=false;          %if settings changed, analysis is invalid
if chline~=Settings.BaselineUpper
    CutoffsSettingsChanged=true;
    Settings.BaselineUpper=chline;
end
if clline~=Settings.BaselineLower
    CutoffsSettingsChanged=true;
    Settings.BaselineLower=clline;
end
if sdposline~=Settings.SpikesUpper
    CutoffsSettingsChanged=true;
    Settings.SpikesUpper=sdposline;
end
if sdnegline~=Settings.SpikesLower
    CutoffsSettingsChanged=true;
    Settings.SpikesLower=sdnegline;
end
Settings.NoiseUpper=pchline;
Settings.NoiseLower=nclline;
Settings.CutoffsMethod=get(handles.CutoffsMethod,'String');

%If any cutoff settings were changed, reset the analysis of all data files. 
%Analysis no longer valid with different parameters
if CutoffsSettingsChanged
    for i=1:Info.NumberOfFiles
        Info.FileInfo(i).SpikesPositive=[];
        Info.FileInfo(i).SpikesNegative=[];
        Info.FileInfo(i).psn=[];
        Info.FileInfo(i).nsn=[];
        Info.FileInfo(i).psw=[];
        Info.FileInfo(i).nsw=[];
        Info.FileInfo(i).Noise=[];
        Info.FileInfo(i).EventsPos=[];
        Info.FileInfo(i).EventsNeg=[];
        Info.FileInfo(i).EventsCombined=[];
    end
end

%Update the settings and the display
Info.Settings=Settings;
set(handles.MenuFiles,'UserData',Info);
RefreshDisplay(handles);



%%% -----------------------------------------------------------------
%%% Step 3: Analyze Events and Noise
%%% -----------------------------------------------------------------

function FileExists=LoadAnalysis(handles)
%Load existing analysis if present. Update values accordingly
    Info=get(handles.MenuFiles,'UserData');
    outputfolder=Info.Settings.outputfolder;
    inputfolder=Info.CurrentDirectory;
    channel=Info.CurrentChannel;
    AnalysisFileName = char(fullfile(inputfolder, outputfolder, ['Analysis_ch ' num2str(channel) '.mat']));  
        
    %Try to load in any existing analysisfile
    try
        AnalysisFile=load(AnalysisFileName);
        FileExists=1;
    catch
        FileExists=0;
    end
        
    if FileExists
        %Load in the information stored so far
        LoadedSettings=AnalysisFile.Info.Settings;
        Info=AnalysisFile.Info;
        set(handles.MenuFiles,'UserData',Info);
        
        %Update the GUI on cutoffs, using the chosen cutoff from settings,
        %rather than the one made for this trace
        set(handles.BaselineUpper,'String',num2str(Info.Settings.BaselineUpper));
        set(handles.BaselineLower,'String',num2str(Info.Settings.BaselineLower));
        set(handles.SpikesUpper,'String',num2str(Info.Settings.SpikesUpper));
        set(handles.SpikesLower,'String',num2str(Info.Settings.SpikesLower));
        try
            set(handles.CutoffsMethod,'String',Info.Settings.CutoffsMethod);
        catch
            set(handles.CutoffsMethod,'String','Method unknown');
        end
    end

    % --- Executes on button press in MenuAnalyze.
function MenuAnalyze_Callback(hObject, eventdata, handles)
% hObject    handle to MenuAnalyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Load in the analysis parameters chosen
Info=get(handles.MenuFiles,'UserData'); 
DoAllFiles=get(handles.MenuAllFiles,'Value');   %Do 1 file, or all files
DoNoise=get(handles.MenuFindNoise,'Value');     %Analyze trace to find noise?
DoEvents=get(handles.MenuFindSeizures,'Value'); %Analyze trace to find events?
DoGraph=1; %Always update the GUI display -> could change later to not do this
Info.Settings.DoShowAnalysis=1;                 %Changed internally to keep track of whether or not there is new analysis to display

%Create output directory if it does not yet exists
inputfolder=Info.CurrentDirectory;
outputfolder=Info.Settings.outputfolder;
try
    mkdir(char(fullfile(inputfolder,outputfolder)));   
    DirectoryExists=0;
catch
    DirectoryExists=1;
end
        

%Determine which files need to be analyzed
if (DoAllFiles)
    startfile=1;
    endfile=Info.NumberOfFiles;
else
    startfile=Info.ActiveFile;
    endfile=Info.ActiveFile;
end       

%Noise & Event collection
TotalNoiseList=[];
TotalEventList=[];

%Loop through all files that need to be analyzed
for FileNumber=startfile:endfile
    
  %analyze the filenumber, as long as it is in the inclusion list.  
  if 1 %TO DO: INTENDED TO GIVE OPT OUT OPTION FOR SPECIFIC FILES
    %Load in data
    set(handles.ProgressUpdate,'String',[num2str(FileNumber) ': Load Data']);
    data=LoadData(FileNumber,1,handles);                                                 %Load in data, without updating buffer
    Info=get(handles.MenuFiles,'UserData'); 
    Info.ActiveFile=FileNumber;
    
    %Find ALL spike information
    set(handles.ProgressUpdate,'String',[num2str(FileNumber) ': Find Spikes']);
    [posspikes,negspikes,psn,nsn,psw,nsw]=SpikeFinder(data,Info.FileInfo(FileNumber).Fs,Info.Settings,0);  %Find positive and negative spikes
    Info.FileInfo(FileNumber).SpikesPositive=posspikes;     %Store positive spikes
    Info.FileInfo(FileNumber).SpikesNegative=negspikes;     %Store negative spikes
    Info.FileInfo(FileNumber).psn=psn;
    Info.FileInfo(FileNumber).nsn=nsn;
    Info.FileInfo(FileNumber).psw=psw;
    Info.FileInfo(FileNumber).nsw=nsw;
    set(handles.MenuFiles,'UserData',Info);
    
    %Find all Noise
    if DoNoise
       set(handles.ProgressUpdate,'String',[num2str(FileNumber) ': Find Noise']);
       ShowPopupNoise=get(handles.PopupNoise,'Value');
       noise=FindNoise(data,posspikes,negspikes,Info.FileInfo(FileNumber).Fs,Info.Settings,ShowPopupNoise);
       Info.FileInfo(FileNumber).Noise=noise;              %Stores noise
       set(handles.MenuFiles,'UserData',Info);
    else
       noise=[];
    end
    
    %Find All Events
    evt_stats=[];
    evt_stats_with_header=[];
    if DoEvents
       set(handles.ProgressUpdate,'String',[num2str(FileNumber) ': Find Events']);
       [posevents,negevents]=FindEvents(data,posspikes,negspikes,psn,nsn,psw,nsw,Info.FileInfo(FileNumber).Fs,Info.Settings,0); 
       Info.FileInfo(FileNumber).EventsPos=posevents;              %Stores events
       Info.FileInfo(FileNumber).EventsNeg=negevents;              %Stores events

       %Combine positive and negative events,and incorporate noise
       set(handles.ProgressUpdate,'String',[num2str(FileNumber) ': Merge Events']);
       eventcutoff=Info.Settings.min_szre_windw*Info.FileInfo(FileNumber).Fs;
       noisecutoff=Info.Settings.noisecutoff*Info.FileInfo(FileNumber).Fs;
       eventglue=Info.Settings.eventglue*Info.FileInfo(FileNumber).Fs;
       dwn=Info.Settings.deal_with_noise;
       dwm=Info.Settings.deal_with_merge;
       [combinedevents,noise]=MergeEvents(posevents,negevents,noise,dwm,dwn,eventcutoff,noisecutoff,eventglue);
       Info.FileInfo(FileNumber).Noise=noise;                      %Update noise after removing short ones
       Info.FileInfo(FileNumber).EventsCombined=combinedevents;    %combine events     
       Info.Settings.DoShowAnalysis=true;                          %Analysis performed, display
       set(handles.MenuFiles,'UserData',Info);
       
       %Calculate statistics on the combined events
       nr_of_events=size(combinedevents,1);
       filestart=Info.FileInfo(FileNumber).StartTime;
       fileend=Info.FileInfo(FileNumber).EndTime;
       filelength=etime(datevec(fileend),datevec(filestart));
       fs=Info.FileInfo(FileNumber).Fs;

       statsheader={'Start Date',' ',' ','Start Time',' ',' ','Finish Date',' ',' ','Finish Time',' ',' ','Duration (s)',...
                    '#+ Spikes', '#- Spikes', '#+W Spikes','#-W Spikes',...
                    '#+N Spikes','#-N Spikes','Avg(+ Width)','Avg(- Width)',...
                    '%+N Spikes','%-N Spikes','+N Spikes/s','-N Spikes/s','+ Spikes/s','- Spikes/s',...
                    'File Number','FS','Date Last Analyzed'};
       
       %Loop through all combined events in the file
       for i=1:nr_of_events
           time_start=filestart+Info.FileInfo(FileNumber).EventsCombined(i,1)/(24*60*60*fs);   %start time
           [y_start,m_start,d_start,h_start,mn_start,s_start]=datevec(time_start);             %components
           time_finish=filestart+Info.FileInfo(FileNumber).EventsCombined(i,2)/(24*60*60*fs);  %end time
           [y_finish,m_finish,d_finish,h_finish,mn_finish,s_finish]=datevec(time_finish);      %components
           duration=etime(datevec(time_finish),datevec(time_start));                           %duration
           
           %Run spikefinder on the found event alone to get counts of spikes.
           [evt_p,evt_n,evt_psn,evt_nsn,evt_psw,evt_nsw]=SpikeFinder(data(combinedevents(i,1):combinedevents(i,2)),Info.FileInfo(FileNumber).Fs,Info.Settings,0);  %Find positive and negative spikes
           number_of_positive_spikes=size(evt_p,1);
           number_of_negative_spikes=size(evt_n,1);
           number_of_positive_wide_spikes=size(evt_psw,1);
           number_of_negative_wide_spikes=size(evt_nsw,1);
           number_of_positive_narrow_spikes=size(evt_psn,1);
           number_of_negative_narrow_spikes=size(evt_nsn,1);
           perc_narrow_spike_pos=number_of_positive_narrow_spikes/number_of_positive_spikes;
           perc_narrow_spike_neg=number_of_negative_narrow_spikes/number_of_negative_spikes;
           narrow_spike_density_pos=number_of_positive_narrow_spikes/duration;
           narrow_spike_density_neg=number_of_negative_narrow_spikes/duration;
           spike_density_pos=number_of_positive_spikes/duration;
           spike_density_neg=number_of_negative_spikes/duration;
           if isempty(evt_p)
               avg_spike_width_pos=0;
           else
               avg_spike_width_pos=mean(evt_p(:,3));
           end
           if isempty(evt_n)
               avg_spike_width_neg=0;
           else
               avg_spike_width_neg=mean(evt_n(:,3));
           end
           
           newstat=[y_start, m_start, d_start, h_start, mn_start, s_start, y_finish, m_finish, d_finish, h_finish, mn_finish, s_finish, duration,...
                    number_of_positive_spikes, number_of_negative_spikes,number_of_positive_wide_spikes,number_of_negative_wide_spikes,...
                    number_of_positive_narrow_spikes,number_of_negative_narrow_spikes,avg_spike_width_pos,avg_spike_width_neg,...
                    perc_narrow_spike_pos,perc_narrow_spike_neg,narrow_spike_density_pos,narrow_spike_density_neg,spike_density_pos,spike_density_neg,...
                    FileNumber,fs,now()];
           evt_stats=[evt_stats; newstat];
       end
       if ~isempty(evt_stats)
           TimeInEvent=sum(evt_stats(:,13));
       else
           TimeInEvent=0;
       end
       PercInEvent=TimeInEvent/filelength;
       Info.FileInfo(FileNumber).TimeInEvent=TimeInEvent;
       Info.FileInfo(FileNumber).PercInEvent=PercInEvent;
       Info.FileInfo(FileNumber).EventStats=evt_stats;
              
       %Determine the number of spikes NOT in noise 
       if isempty(noise) %ALL SPIKES COUNT
           totalposspikes=size(psn,1);
           totalnegspikes=size(nsn,1);
           
           totalnoise=0;
       else %ONLY SPIKES NOT IN NOISE COUNT
           inoise=1;totalnoise=0;
           iposspikes=1;totalposspikes=0;
           inegspikes=1;totalnegspikes=0;
           while inoise<=length(noise(:,1))    %Go through noise events
               
               %Count the number of positive spikes prior to this noise
               while ~isempty(psn) && iposspikes<=size(psn,1) && psn(iposspikes,1)<noise(inoise,1)
                   totalposspikes=totalposspikes+1;
                   iposspikes=iposspikes+1;
               end
               %Count the number of negative spikes prior to this noise
               while ~isempty(nsn) && inegspikes<=size(nsn,1) && nsn(inegspikes,1)<noise(inoise,1)
                   totalnegspikes=totalnegspikes+1;
                   inegspikes=inegspikes+1;
               end

               %There are now no more spikes left prior to noise. See if
               %there are any IN the noise
               while ~isempty(psn) && iposspikes<=size(psn,1) && psn(iposspikes,1)<noise(inoise,2)
                   iposspikes=iposspikes+1;
               end
               while ~isempty(nsn) && inegspikes<=size(nsn,1) && nsn(inegspikes,1)<noise(inoise,2)
                   inegspikes=inegspikes+1;
               end
               
               %There are now no more spikes left prior or during to noise
               %Go to next noise
               totalnoise=totalnoise+noise(inoise,2)-noise(inoise,1);
               inoise=inoise+1;
           end
           
           %There might be spikes left after the last noise event
           if ~isempty(psn)
               totalposspikes=totalposspikes+size(psn,1)-iposspikes+1;
           else
               totalposspikes=0;
           end
           if ~isempty(nsn)
               totalnegspikes=totalnegspikes+size(nsn,1)-inegspikes+1;
           else
               totalnegspikes=0;
           end
       end
       %Calculate the #spike-statistics
       time_in_noise=totalnoise/fs;
       time_without_noise=filelength-time_in_noise;     %time in s
       posspikerate=totalposspikes/time_without_noise;  %#/s
       negspikerate=totalnegspikes/time_without_noise;  %#/s
       
       Info.FileInfo(FileNumber).TotalTime=filelength;
       Info.FileInfo(FileNumber).TimeInNoise=time_in_noise;
       Info.FileInfo(FileNumber).totalposspikes=totalposspikes;
       Info.FileInfo(FileNumber).totalnegspikes=totalnegspikes;
       Info.FileInfo(FileNumber).posspikerate=posspikerate;
       Info.FileInfo(FileNumber).negspikerate=negspikerate;
       
       %Update the last analyzed time
       Info.FileInfo(FileNumber).LastAnalyzed=now();
       set(handles.MenuFiles,'UserData',Info);
    end
    
    
    
    DoSave=1;       %%%???
    if (DoSave)
        %Specify location of the analysis file
        outputfolder=Info.Settings.outputfolder;
        inputfolder=Info.CurrentDirectory;
        channel=Info.CurrentChannel;
        AnalysisFileName = char(fullfile(inputfolder, outputfolder, ['Analysis_ch ' num2str(channel) '.mat']));  
        AnalysisFileCSV = char(fullfile(inputfolder, outputfolder, ['Analysis_ch ' num2str(channel) '.csv']));  
        AnalysisFileXLS = char(fullfile(inputfolder, outputfolder, ['Analysis_ch ' num2str(channel) '.xls']));  
              
        %Try to load in any existing analysisfile
        try
            AnalysisFile=load(AnalysisFileName);
            FileExists=1;
        catch
            FileExists=0;
        end
        
        if FileExists
           %Load in the information stored so far
           LoadedSettings=AnalysisFile.Info.Settings;
            if ~isequal(LoadedSettings,Info.Settings)
                %h=warndlg('Analysis file will be overwritten');
                overwrite=1;
            else
                %append new information
                overwrite=0;
            end
        else
            overwrite=1;
        end
        %No file saved yet, so store current info       

        
        if ~overwrite
            LoadedFileInfo=AnalysisFile.Info.FileInfo;
            LoadedFileInfo(Info.ActiveFile)=Info.FileInfo(Info.ActiveFile);
            Info.FileInfo=LoadedFileInfo;           
        end
                
        
        %Save all matlab workspace variables to file
        save(AnalysisFileName,'Info'); 
        
        %Save the information on the events to a CSV file
        %This includes the stats on individual events
        try
            %If there is already a file, append
            olddata=csvread(AnalysisFileCSV);
            if ~Info.Settings.SettingsChangedSinceLastSave
                evt_stats=[olddata; evt_stats];
            end
        catch
            
        end
        %Save the CSV file
        csvwrite(AnalysisFileCSV,evt_stats);
        
        
        %Save the information on the events to a XLS file
        %This includes the stats on individual events
        try
            %If there is already a file, append
            [~,~,olddata]=xlsread(AnalysisFileXLS);
            if ~Info.Settings.SettingsChangedSinceLastSave
                evt_stats_with_header=[olddata; evt_stats];
            else
                evt_stats_with_header=[statsheader; num2cell(evt_stats)];
            end
        catch
            evt_stats_with_header=[statsheader; num2cell(evt_stats)];
        end
        %Save the XLS file
        xlswrite(AnalysisFileXLS,evt_stats_with_header);
        
        %Settings have not been changed since last save, since last save
        %just occurred.
        Info.Settings.SettingsChangedSinceLastSave=0; 
        set(handles.MenuFiles,'UserData',Info);
    end
    
    %Update display
    if DoGraph
        set(handles.ProgressUpdate,'String',[num2str(FileNumber) ': Update Graph']);
        
        %Update Graph
        DisplayFile(FileNumber,handles);
        
        %bookkeeping 
        EventSummaryName = char(fullfile(inputfolder, outputfolder, ['Channel ' num2str(Info.CurrentChannel) ' - ' 'File number ' num2str(FileNumber) ' - ' 'EventSummary.png']))  
        saveas(gcf,EventSummaryName);
        
    end
 
    set(handles.ProgressUpdate,'String',[num2str(FileNumber) ': Done!']);
  end
end


%%% -----------------------------------------------------------------
%%% Control and Refresh Display 
%%% -----------------------------------------------------------------


function ShowCutoffs(upperrange,handles)
    ShowIt=get(handles.MenuShowIt,'Value');
    if ShowIt
        axes(handles.MenuGraph);
        hold on
        chline=str2num(get(handles.BaselineUpper,'String'));
        clline=str2num(get(handles.BaselineLower,'String'));
        pchline=str2num(get(handles.SpikesUpper,'String'));
        nclline=str2num(get(handles.SpikesLower,'String'));
    
        thrbu=plot(handles.MenuGraph,[0,upperrange],[chline,chline],'m'); %high threshold (e.g. 90% data)
        thrbl=plot(handles.MenuGraph,[0,upperrange],[clline,clline],'m'); %low threshold (e.g. 90% data)
        thrsu=plot(handles.MenuGraph,[0,upperrange],[pchline,pchline],'y');%spike threshold (e.g. 2SD +data)
        thrsl=plot(handles.MenuGraph,[0,upperrange],[nclline,nclline],'y');%spike threshold (e.g. 2SD -data)
    
        hold off
    end

function ShowAnalysis(data,handles)
    %Obtain all the data
    Info=get(handles.MenuFiles,'UserData');                         
    FileNumber=Info.ActiveFile;
    fs=Info.FileInfo(FileNumber).Fs;
    domain=[0:1/fs:(length(data)-1)/fs]'; 
    posevents=Info.FileInfo(FileNumber).EventsPos;             
    negevents=Info.FileInfo(FileNumber).EventsNeg; 
    combinedevents=Info.FileInfo(FileNumber).EventsCombined; 
    noise=Info.FileInfo(FileNumber).Noise; 
    psn=Info.FileInfo(FileNumber).psn;
    nsn=Info.FileInfo(FileNumber).nsn;
    psw=Info.FileInfo(FileNumber).psw;
    nsw=Info.FileInfo(FileNumber).nsw;
    
    axes(handles.MenuGraph);
    hold on;
    
    %Add spikes to plot (if desired)
    %Black dots for regular spikes, pink dots for wide spikes
    ShowSpikes=get(handles.showspikes,'Value');
    if ShowSpikes
        if ~isempty(psn)
            plot(handles.MenuGraph,domain(psn(:,1)),psn(:,2),'k.');
        end
    
        if ~isempty(nsn)
            plot(handles.MenuGraph,domain(nsn(:,1)),nsn(:,2),'k.');
        end
    
        if ~isempty(psw)
            plot(handles.MenuGraph,domain(psw(:,1)),psw(:,2),'m.');
        end
        if ~isempty(nsw)
            plot(handles.MenuGraph,domain(nsw(:,1)),nsw(:,2),'m.');
        end
    end 
    
    %Add events to plot
    if ~isempty(posevents)
        loc=0.5*Info.Settings.BaselineUpper;
        for i=1:length(posevents(:,1))    
            plot([domain(posevents(i,1)) domain(posevents(i,2))], [loc loc],'g','LineWidth',2);
        end
        lastpos=plot([domain(posevents(i,1)) domain(posevents(i,2))], [loc loc],'g','LineWidth',2);
    end;
    
    if ~isempty(negevents)
        loc=0.5*Info.Settings.BaselineLower;
        for i=1:length(negevents(:,1))
            plot([domain(negevents(i,1)) domain(negevents(i,2))], [loc loc],'r','LineWidth',2);
        end
        lastneg=plot([domain(negevents(i,1)) domain(negevents(i,2))], [loc loc],'r','LineWidth',2);
    end
        
    if ~isempty(combinedevents)
        loc=0.2*Info.Settings.BaselineUpper;
        for i=1:length(combinedevents(:,1))
            plot([domain(combinedevents(i,1)) domain(combinedevents(i,2))], [loc loc],'m','LineWidth',2);
        end
        lastevent=plot([domain(combinedevents(i,1)) domain(combinedevents(i,2))], [loc loc],'m','LineWidth',2);
    end 
            
    %Add noise to plot
    if ~isempty(noise)
        for i=1:length(noise(:,1))
            plot([domain(noise(i,1)) domain(noise(i,2))],[0 0],'y','LineWidth',2);
        end
        lastnoise=plot([domain(noise(i,1)) domain(noise(i,2))],[0 0],'y','LineWidth',2);
    end
    
    hold off

    
function RefreshDisplay(handles)
%Updates the display for the current file. Will display data, analysis, and
%thresholds. Called whenever a display option is toggled or when analysis
%is run.
%
    %Get file information on current file
    Info=get(handles.MenuFiles,'UserData');          %get user data from handles
    ActiveFile=Info.ActiveFile;                      %load active file number
    data=LoadData(ActiveFile,1,handles);
    Info=get(handles.MenuFiles,'UserData');          %get user data from handles after update
    fs=Info.FileInfo(ActiveFile).Fs;
    domain=[0:1/fs:(length(data)-1)/fs]';            %specify domain
    upperrange=Info.FileInfo(ActiveFile).upperrange;

    %Plot raw data
    axes(handles.MenuGraph);
    plot(handles.MenuGraph,domain,data);

    %Plot threshold lines (if desired through check box)
    ShowIt=get(handles.MenuShowIt,'Value');
    if ShowIt
        ShowCutoffs(upperrange,handles);
    end

    %Plot event analysis (if desired and available) 
    try 
        %See if we opted to show analysis in graph
        DoShowAnalysis=Info.Settings.DoShowAnalysis;
    catch
        %If we didn't opt for it, we won't
    end
    if exist('DoShowAnalysis','var')
        if DoShowAnalysis
            ShowAnalysis(data,handles);
        end
    end

function DisplayFile(ActiveFile,handles)
%Loads in the active file. Then updates the information in the GUI about
%the file, and display's analysis if appropriate. 

%Extract info about active file
Info=get(handles.MenuFiles,'UserData');
NumberOfFiles=Info.NumberOfFiles;
CurrentDirectory=Info.CurrentDirectory;
AllFiles=Info.AllFiles;
CurrentChannel=Info.CurrentChannel;
set(handles.MenuFiles,'UserData',Info);


%If user selects previous file, while we were on first,
%or next file when we were on last, loop around if necessary. 
%If user chooses file number out of bounds, display warning
if ActiveFile>NumberOfFiles+1
    %Return error if out of bound
    set(handles.MenuWarning,'String',['There are only' num2str(NumberOfFiles) 'files']');
    return
elseif ActiveFile==0
   ActiveFile=NumberOfFiles; %Loop to end
elseif ActiveFile==NumberOfFiles+1
    ActiveFile=1; %Loop to beginning
end

%Update Filename and index on GUI
set(handles.MenuFileIndex,'String',num2str(ActiveFile));
set(handles.MenuFileTotal,'String',strcat('/',num2str(NumberOfFiles)));
newname=strcat('File ',num2str(ActiveFile),': ',AllFiles(ActiveFile));
set(handles.MenuFileName,'String',newname);
Info.ActiveFile=ActiveFile;
set(handles.MenuFiles,'UserData',Info);


%If file has not yet been loaded this session, or is not buffered, load it into buffer
data=LoadData(ActiveFile,1,handles);    %Load data, add to buffer if desired

%If analysis of the graph has already occurred, and settings have not been
%changed since, then load in those values for display as well.
try
    SettingsChangedSinceLastSave=Info.Settings.SettingsChangedSinceLastSave;
catch
end

if exist('SettingsChangedSinceLastSave','var')
   if SettingsChangedSinceLastSave
       %Analysis has to be redone. Cannot display it in current graph
       Info.Settings.DoShowAnalysis=false;
   else
       %Analysis is already present. No need to recompute, just display
       %However, maybe it wasn't desired, so only set it to true if it
       %didn't yet exist
       try 
           %see if it exists already
           dsa=Info.Settings.DoShowAnalysis;
       catch
           %if not, make it and set it to true
           Info.Settings.DoShowAnalysis=true;
       end
   end
   
else
   %Analysis has never been done, and cannot be displayed in current graph
   Info.Settings.DoShowAnalysis=false;
   
end
set(handles.MenuFiles,'UserData',Info);

%Display graph and optionally cutoffs
RefreshDisplay(handles);

function MenuShowIt_Callback(hObject, eventdata, handles)
% hObject    handle to MenuCutoffs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MenuCutoffs

%Call the show of the cutoff lines
Info=get(handles.MenuFiles,'UserData');
ShowIt=get(handles.MenuShowIt,'Value'); 
if ShowIt
    ShowCutoffs(Info.FileInfo(Info.ActiveFile).upperrange,handles);    
else
    try
        set(handles.thrbu,'Visible','off');
        set(handles.thrbl,'Visible','off');
        set(handles.thrsu,'Visible','off');
        set(handles.thrsl,'Visible','off');
    end
end
RefreshDisplay(handles);


% --- Executes on button press in MenuPrevious.
function MenuPrevious_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrevious (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Displays the previous file and perform calculations if desired
Info=get(handles.MenuFiles,'UserData');
DisplayFile(Info.ActiveFile-1,handles);

% --- Executes on button press in MenuNext.
function MenuNext_Callback(hObject, eventdata, handles)
% hObject    handle to MenuNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Display the next file, and perform calculations if desired
Info=get(handles.MenuFiles,'UserData');
DisplayFile(Info.ActiveFile+1,handles);



% --------------------------------------------------------------------
% SETTINGS
% --------------------------------------------------------------------


% --- Executes on button press in set_default.
function set_default_Callback(hObject, eventdata, handles)
% hObject    handle to set_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Restores the default settings
AllInfo=get(handles.MenuFiles,'UserData');
AllInfo.Settings=DefaultSettings();
AllInfo.Settings.SettingsChangedSinceLastSave=1; %assume settings changed
set(handles.MenuFiles,'UserData',AllInfo);
settingstool_ClickedCallback(hObject, eventdata, handles);

% --- Executes on button press in SettingsSave.
function SettingsSave_Callback(hObject, eventdata, handles)
% hObject    handle to SettingsSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Hide Settings Menu
set(handles.SettingsPanel,'Visible','off');

%Apply Changes
AllInfo=get(handles.MenuFiles,'UserData');  %obtain the structure AllInfo
Settings=AllInfo.Settings; 
Settings.SettingsChangedSinceLastSave=1; %assume settings changed

%update file settings
Settings.outputfolder=get(handles.set_output,'String');
Settings.buffer=str2num(get(handles.set_buffer,'String'));

%update automatic threshold settings
Settings.all_dist_pos=str2num(get(handles.set_posprespike,'String'));
Settings.all_dist_neg=str2num(get(handles.set_negprespike,'String'));
Settings.all_thrsh_pos=str2num(get(handles.set_pospreamp,'String'));
Settings.all_thrsh_pos=str2num(get(handles.set_negpreamp,'String'));
Settings.cl=str2num(get(handles.set_cutoff_all_low,'String'));
Settings.ch=str2num(get(handles.set_cutoff_all_high,'String'));
Settings.pcl=str2num(get(handles.set_cutoff_pos_low,'String'));
Settings.pch=str2num(get(handles.set_cutoff_pos_high,'String'));
Settings.ncl=str2num(get(handles.set_cutoff_neg_low,'String'));
Settings.nch=str2num(get(handles.set_cutoff_neg_high,'String'));
Settings.number_of_sd=str2num(get(handles.set_sd,'String'));

%update spike detection settings
Settings.dist_pos=str2num(get(handles.set_pos_dist,'String'));
Settings.dist_neg=str2num(get(handles.set_neg_dist,'String'));
Settings.ht_perc=str2num(get(handles.set_height,'String'));
Settings.req_width=str2num(get(handles.set_width,'String'));

%update event detection settings
Settings.min_szre_lngth=str2num(get(handles.set_szr_length,'String'));
Settings.min_nr_spikes=str2num(get(handles.set_szr_nrspikes,'String'));
Settings.min_szre_windw=str2num(get(handles.set_szr_nrseconds,'String'));
Settings.eventglue=str2num(get(handles.set_glue,'String'));
Settings.connect_events=get(handles.set_szr_glue,'Value');
if get(handles.set_events_union,'Value')
    dwm=0;
elseif get(handles.set_events_intersect,'Value')
    dwm=1;
elseif get(handles.set_events_pos,'Value')
    dwm=2;
else
    dwm=3;
end
Settings.deal_with_merge=dwm;



%update noise detection settings
Settings.glue_noise=str2num(get(handles.set_joinnoise,'String'));
Settings.noisecutoff=str2num(get(handles.set_ignorenoise,'String'));
if get(handles.set_noise_ignore,'Value')
    dwn=0;
elseif get(handles.set_noise_split,'Value')
    dwn=1;
else
    dwn=2;
end
Settings.deal_with_noise=dwn;
Settings.movingwin=[str2num(get(handles.set_noise_window,'String')),str2num(get(handles.set_noise_overlap,'String'))];

AllInfo.Settings=Settings;
set(handles.MenuFiles,'UserData',AllInfo);
set(handles.SettingsPanel,'Visible','off');


function settingstool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to settingstool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Open up settings menu
AllInfo=get(handles.MenuFiles,'UserData');  %obtain the structure AllInfo
Settings=AllInfo.Settings; 

%update file settings
set(handles.set_output,'String',Settings.outputfolder);
set(handles.set_buffer,'String',num2str(Settings.buffer));

%update automatic threshold settings
set(handles.set_posprespike,'String',num2str(Settings.all_dist_pos));
set(handles.set_negprespike,'String',num2str(Settings.all_dist_neg));
set(handles.set_pospreamp,'String',num2str(Settings.all_thrsh_pos));
set(handles.set_negpreamp,'String',num2str(Settings.all_thrsh_pos));
set(handles.set_cutoff_all_low,'String',num2str(Settings.cl));
set(handles.set_cutoff_all_high,'String',num2str(Settings.ch));
set(handles.set_cutoff_pos_low,'String',num2str(Settings.pcl));
set(handles.set_cutoff_pos_high,'String',num2str(Settings.pch));
set(handles.set_cutoff_neg_low,'String',num2str(Settings.ncl));
set(handles.set_cutoff_neg_high,'String',num2str(Settings.nch));
set(handles.set_sd,'String',num2str(Settings.number_of_sd));

%update spike detection settings
set(handles.set_pos_dist,'String',num2str(Settings.dist_pos));
set(handles.set_neg_dist,'String',num2str(Settings.dist_neg));
set(handles.set_height,'String',num2str(Settings.ht_perc));
set(handles.set_width,'String',num2str(Settings.req_width));

%update event detection settings
set(handles.set_szr_length,'String',num2str(Settings.min_szre_lngth));
set(handles.set_szr_nrspikes,'String',num2str(Settings.min_nr_spikes));
set(handles.set_szr_nrseconds,'String',num2str(Settings.min_szre_windw));
set(handles.set_szr_glue,'Value',Settings.connect_events);
set(handles.set_glue,'String',num2str(Settings.eventglue));

%Update toggle box about how to deal with noise 
dwn=Settings.deal_with_noise;
set(handles.set_noise_ignore,'Value',dwn==0);
set(handles.set_noise_split,'Value',dwn==1);
set(handles.set_noise_skip,'Value',dwn==2);
set(handles.set_noise_window,'String',num2str(Settings.movingwin(1)));
set(handles.set_noise_overlap,'String',num2str(Settings.movingwin(2)));


%Update the toggle box about how to deal with +/- event merges
dwm=Settings.deal_with_merge;
set(handles.set_events_union,'Value',dwm==0);
set(handles.set_events_intersect,'Value',dwm==1);
set(handles.set_events_pos,'Value',dwm==2);
set(handles.set_events_neg,'Value',dwm==3);  

%update noise detection settings
set(handles.set_joinnoise,'String',num2str(Settings.glue_noise));
set(handles.set_ignorenoise,'String',num2str(Settings.noisecutoff));
set(handles.SettingsPanel,'Visible','on');


% --------------------------------------------------------------------
% GUI FUNCTIONS
% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function MenuDir_Callback(hObject, eventdata, handles)
% hObject    handle to MenuDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function MenuDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MenuDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in MenuChannel.
function MenuChannel_Callback(hObject, eventdata, handles)
% hObject    handle to MenuChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function MenuChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MenuChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in MenuOverwriteMeta.
function MenuOverwriteMeta_Callback(hObject, eventdata, handles)
% hObject    handle to MenuOverwriteMeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in MenuLoadAnalysis.
function MenuLoadAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to MenuLoadAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function SpikesUpper_Callback(hObject, eventdata, handles)
% hObject    handle to SpikesUpper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function SpikesUpper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SpikesUpper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SpikesLower_Callback(hObject, eventdata, handles)
% hObject    handle to SpikesLower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function SpikesLower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SpikesLower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BaselineUpper_Callback(hObject, eventdata, handles)
% hObject    handle to BaselineUpper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function BaselineUpper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BaselineUpper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BaselineLower_Callback(hObject, eventdata, handles)
% hObject    handle to BaselineLower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function MenuFileIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BaselineLower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BaselineLower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BaselineLower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MenuHasPositives.
function MenuHasPositives_Callback(hObject, eventdata, handles)
% hObject    handle to MenuHasPositives (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on key press with focus on MenuFileTotal and none of its controls.
function MenuFileIndex_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to MenuFileTotal (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in MenuChangeCurrent.
function MenuChangeCurrent_Callback(hObject, eventdata, handles)
% hObject    handle to MenuChangeCurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ActiveFile=str2num(get(handles.MenuFileIndex,'String'));
DisplayFile(ActiveFile,handles);





    
function MenuFileIndex_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileTotal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MenuFileTotal as text
%        str2double(get(hObject,'String')) returns contents of MenuFileTotal as a double


% --- Executes during object creation, after setting all properties.
function MenuFileTotal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MenuFileTotal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


       





% --- Executes on button press in ChangeStep1.
function ChangeStep1_Callback(hObject, eventdata, handles)
% hObject    handle to ChangeStep1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.HidePanelLeft,'Visible','off');
set(handles.HidePanelMiddle,'Visible','on');
set(handles.ChangeStep2,'Visible','off');
set(handles.HidePanelRight,'Visible','on');



% --- Executes on button press in MenuFindNoise.
function MenuFindNoise_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFindNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MenuFindNoise


% --- Executes on button press in MenuFindSeizures.
function MenuFindSeizures_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFindSeizures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MenuFindSeizures




% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in MenuAllFiles.
function MenuAllFiles_Callback(hObject, eventdata, handles)
% hObject    handle to MenuAllFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MenuAllFiles


% --- Executes on button press in MenuCurrentFile.
function MenuCurrentFile_Callback(hObject, eventdata, handles)
% hObject    handle to MenuCurrentFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MenuCurrentFile


% --- Executes on button press in ChangeStep2.
function ChangeStep2_Callback(hObject, eventdata, handles)
% hObject    handle to ChangeStep2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
%Change visibility
set(handles.HidePanelMiddle,'Visible','off');
set(handles.HidePanelRight,'Visible','on');


%Also load in other values if settings didn't change and there are
%analysis-values already.


% --- Executes on button press in CloseupPrev.
function CloseupPrev_Callback(hObject, eventdata, handles)
% hObject    handle to CloseupPrev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CloseupNext.
function CloseupNext_Callback(hObject, eventdata, handles)
% hObject    handle to CloseupNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in showspikes.
function showspikes_Callback(hObject, eventdata, handles)
% hObject    handle to showspikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showspikes
RefreshDisplay(handles);


% --- Executes on button press in PopupThreshold.
function PopupThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to PopupThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PopupThreshold


% --- Executes on button press in PopupNoise.
function PopupNoise_Callback(hObject, eventdata, handles)
% hObject    handle to PopupNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PopupNoise


%


function set_sd_Callback(hObject, eventdata, handles)
% hObject    handle to set_sd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_sd as text
%        str2double(get(hObject,'String')) returns contents of set_sd as a double


% --- Executes during object creation, after setting all properties.
function set_sd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_sd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_cutoff_neg_low_Callback(hObject, eventdata, handles)
% hObject    handle to set_cutoff_neg_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_cutoff_neg_low as text
%        str2double(get(hObject,'String')) returns contents of set_cutoff_neg_low as a double


% --- Executes during object creation, after setting all properties.
function set_cutoff_neg_low_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_cutoff_neg_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_cutoff_all_low_Callback(hObject, eventdata, handles)
% hObject    handle to set_cutoff_all_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_cutoff_all_low as text
%        str2double(get(hObject,'String')) returns contents of set_cutoff_all_low as a double


% --- Executes during object creation, after setting all properties.
function set_cutoff_all_low_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_cutoff_all_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_cutoff_pos_low_Callback(hObject, eventdata, handles)
% hObject    handle to set_cutoff_pos_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_cutoff_pos_low as text
%        str2double(get(hObject,'String')) returns contents of set_cutoff_pos_low as a double


% --- Executes during object creation, after setting all properties.
function set_cutoff_pos_low_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_cutoff_pos_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_cutoff_all_high_Callback(hObject, eventdata, handles)
% hObject    handle to set_cutoff_all_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_cutoff_all_high as text
%        str2double(get(hObject,'String')) returns contents of set_cutoff_all_high as a double


% --- Executes during object creation, after setting all properties.
function set_cutoff_all_high_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_cutoff_all_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_cutoff_pos_high_Callback(hObject, eventdata, handles)
% hObject    handle to set_cutoff_pos_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_cutoff_pos_high as text
%        str2double(get(hObject,'String')) returns contents of set_cutoff_pos_high as a double


% --- Executes during object creation, after setting all properties.
function set_cutoff_pos_high_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_cutoff_pos_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_cutoff_neg_high_Callback(hObject, eventdata, handles)
% hObject    handle to set_cutoff_neg_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_cutoff_neg_high as text
%        str2double(get(hObject,'String')) returns contents of set_cutoff_neg_high as a double


% --- Executes during object creation, after setting all properties.
function set_cutoff_neg_high_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_cutoff_neg_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_negprespike_Callback(hObject, eventdata, handles)
% hObject    handle to set_negprespike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_negprespike as text
%        str2double(get(hObject,'String')) returns contents of set_negprespike as a double


% --- Executes during object creation, after setting all properties.
function set_negprespike_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_negprespike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_posprespike_Callback(hObject, eventdata, handles)
% hObject    handle to set_posprespike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_posprespike as text
%        str2double(get(hObject,'String')) returns contents of set_posprespike as a double


% --- Executes during object creation, after setting all properties.
function set_posprespike_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_posprespike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_pospreamp_Callback(hObject, eventdata, handles)
% hObject    handle to set_pospreamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_pospreamp as text
%        str2double(get(hObject,'String')) returns contents of set_pospreamp as a double


% --- Executes during object creation, after setting all properties.
function set_pospreamp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_pospreamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_negpreamp_Callback(hObject, eventdata, handles)
% hObject    handle to set_negpreamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_negpreamp as text
%        str2double(get(hObject,'String')) returns contents of set_negpreamp as a double


% --- Executes during object creation, after setting all properties.
function set_negpreamp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_negpreamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_output_Callback(hObject, eventdata, handles)
% hObject    handle to set_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_output as text
%        str2double(get(hObject,'String')) returns contents of set_output as a double


% --- Executes during object creation, after setting all properties.
function set_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_buffer_Callback(hObject, eventdata, handles)
% hObject    handle to set_buffer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_buffer as text
%        str2double(get(hObject,'String')) returns contents of set_buffer as a double


% --- Executes during object creation, after setting all properties.
function set_buffer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_buffer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.SettingsPanel,'Visible','off');


function set_joinnoise_Callback(hObject, eventdata, handles)
% hObject    handle to set_joinnoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_joinnoise as text
%        str2double(get(hObject,'String')) returns contents of set_joinnoise as a double


% --- Executes during object creation, after setting all properties.
function set_joinnoise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_joinnoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_ignorenoise_Callback(hObject, eventdata, handles)
% hObject    handle to set_ignorenoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_ignorenoise as text
%        str2double(get(hObject,'String')) returns contents of set_ignorenoise as a double


% --- Executes during object creation, after setting all properties.
function set_ignorenoise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_ignorenoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_szr_nrspikes_Callback(hObject, eventdata, handles)
% hObject    handle to set_szr_nrspikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_szr_nrspikes as text
%        str2double(get(hObject,'String')) returns contents of set_szr_nrspikes as a double


% --- Executes during object creation, after setting all properties.
function set_szr_nrspikes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_szr_nrspikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_szr_length_Callback(hObject, eventdata, handles)
% hObject    handle to set_szr_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_szr_length as text
%        str2double(get(hObject,'String')) returns contents of set_szr_length as a double


% --- Executes during object creation, after setting all properties.
function set_szr_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_szr_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_szr_nrseconds_Callback(hObject, eventdata, handles)
% hObject    handle to set_szr_nrseconds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_szr_nrseconds as text
%        str2double(get(hObject,'String')) returns contents of set_szr_nrseconds as a double


% --- Executes during object creation, after setting all properties.
function set_szr_nrseconds_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_szr_nrseconds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in set_szr_glue.
function set_szr_glue_Callback(hObject, eventdata, handles)
% hObject    handle to set_szr_glue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function set_width_Callback(hObject, eventdata, handles)
% hObject    handle to set_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function set_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_height_Callback(hObject, eventdata, handles)
% hObject    handle to set_height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function set_height_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_neg_dist_Callback(hObject, eventdata, handles)
% hObject    handle to set_neg_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_neg_dist as text
%        str2double(get(hObject,'String')) returns contents of set_neg_dist as a double


% --- Executes during object creation, after setting all properties.
function set_neg_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_neg_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_pos_dist_Callback(hObject, eventdata, handles)
% hObject    handle to set_pos_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_pos_dist as text
%        str2double(get(hObject,'String')) returns contents of set_pos_dist as a double


% --- Executes during object creation, after setting all properties.
function set_pos_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_pos_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_noise_overlap_Callback(hObject, eventdata, handles)
% hObject    handle to set_noise_overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_noise_overlap as text
%        str2double(get(hObject,'String')) returns contents of set_noise_overlap as a double


% --- Executes during object creation, after setting all properties.
function set_noise_overlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_noise_overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_noise_window_Callback(hObject, eventdata, handles)
% hObject    handle to set_noise_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_noise_window as text
%        str2double(get(hObject,'String')) returns contents of set_noise_window as a double


% --- Executes during object creation, after setting all properties.
function set_noise_window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_noise_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_glue_Callback(hObject, eventdata, handles)
% hObject    handle to set_glue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_glue as text
%        str2double(get(hObject,'String')) returns contents of set_glue as a double


% --- Executes during object creation, after setting all properties.
function set_glue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_glue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
