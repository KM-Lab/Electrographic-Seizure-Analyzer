function ConsiderFileList=AnalyzeDirectory(inputfolder,channel,DoCheck,outputfolder)
%Return a list of all *.mat files in the directory. It will ignore
%subdirectories and the files therein. It will check that all files for the
%given channel have the same parameters (ledon, ledoff, ledactive, fs)
%
%Input: 'inputfolder'         Location of the input files, with respect to
%                             current location
%       'channel'             The channel in the .mat files to consider
%       'DoCheck'             1 = Check for consistency between files
%                             0 = Do not check for consistency
%       'outputfolder'        Desired location of the output files, with
%                             respect to current location
%
%Output: 'ConsiderFileList'   A list of files {} or 0 if there are no files
%
%Calls:         none
%
%Called by:     'MenuPreanalysis_Callback' in 'findallofthem.m'
%

   %Obtain all files in the input directory
   drctry = [inputfolder '/*.mat'];
   matfiles = dir(drctry);
   matfiles([matfiles.isdir]) = [];   %Remove subdirectories from consideration
   numbertoconsider=length(matfiles); %Find Number of files

   %Screen the files for matching parameters
   FileList={matfiles.name};
   ConsiderFileList=FileList; %%%!Future change enables opt out of certain files

   %Make subdirectory for output files if it does not yet exist
   try
      mkdir(char(fullfile(inputfolder,outputfolder)));   
   end
   PreAnalysisFileName = char(fullfile(inputfolder, outputfolder, ['PreAnalysis_ch' num2str(channel) '.mat']));  
   try
      %Load in existing preanalysis file for the given channel if exists
      PreAnalysisFile=load(PreAnalysisFileName);
      FileExists=1;
   catch
      %File does not exist yet, so create it by analyzing files        
      FileExists=0;
   end
   
   %Preanalysis checks if files are internally consistent, and have the same stimulation
   %parameters and subjects associated with them. If we have not yet done
   %this for the given directory, or want to redo it, do so.
   if DoCheck    
      %Analyze the metadata, and check for consistency
      Metadata=struct;
      InconsistencyFlag=false;
      InconsistencyWarnings=[];
      ConsiderFileList={}; %Reset since we are analyzing from scratch
    
      %Loop through all files in the directory
      for i=1:numbertoconsider
         file_name=char(fullfile(inputfolder, FileList(i)));
         Metadata(i).name=FileList(i);
         Metadata(i).include=1;

         temp=load(file_name,'number_of_channels');Metadata(i).NrChannels=temp.number_of_channels;
         if channel<=Metadata(i).NrChannels
            cfl=length(ConsiderFileList); ConsiderFileList(cfl+1)=FileList(i);
            temp=load(file_name,'fs'); Metadata(i).fs=temp.fs;
            temp=load(file_name,'ledon'); Metadata(i).ledon=temp.ledon(channel);
            temp=load(file_name,'ledoff'); Metadata(i).ledoff=temp.ledoff(channel);
            temp=load(file_name,'ledactive'); Metadata(i).ledactive=temp.ledactive(channel);
    
            %Check for consistency of metadata
            if not(InconsistencyFlag)
               if Metadata(i).fs~=Metadata(i).fs
                  InconsistencyFlag=true;
                  InconsistencyWarnings='fs';
               end
               if Metadata(i).ledon~=Metadata(1).ledon
                  InconsistencyFlag=true;
                  InconsistencyWarnings='ledon';
               end
               if Metadata(i).ledoff~=Metadata(1).ledoff
                  InconsistencyFlag=true;
                  InconsistencyWarnings='ledoff';
               end
               if Metadata(i).ledactive~=Metadata(1).ledactive
                  InconsistencyFlag=true;
                  InconsistencyWarnings='ledactive';
               end 
            end
         end         
      end
      
      %save the metadata
      save(PreAnalysisFileName,'Metadata');
    
      %Report warnings, if any
      if and(InconsistencyFlag,DoCheck)
         Warn=warndlg(['The selected folder contains at least the following inconsistency for channel ' num2str(channel) ':' InconsistencyWarnings '\n You can see the PreAnalysis File in the created subdirectory for more information. You can also exclude files by changing this document.']);
      end 
   end
    
   %Check if there were any files for the given channel. If not, code this by
   %returning 0, in stead of {}
   if length(FileList)>0 && length(ConsiderFileList)==0
      ConsiderFileList=0; 
   end
   
end


