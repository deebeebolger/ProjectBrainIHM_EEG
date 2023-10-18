eeglab% Date: 01-03-2021   Programmed by: Deirdre BOLGER
% Script to segment the continuous data for the PEPs project data.
% The continuous datasets are segmented into "congruent" and, if present,
% "incongruent" and saved in a folder depending on the title of the video
% that corresponds to the dataset in question.
%**************************************************************************

sfiles = dir(strcat(dirpath,'*.set'));
findsets=find(~[sfiles.isdir]);        
Allsets = {sfiles(findsets).name}; 

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename',Allsets,'filepath',dirpath);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
eeglab redraw;

% Dialogue opens to query the time limits of the segmentation. 
deflts = {'0.8','0','-0.25'};
time_lim = inputdlg({'Enter trial upper limit (s)','Enter baseline upper limit (s)','Enter baseline lower limit (s)'},'Enter time limits',...
    [1 50;1 50;1 50],deflts);
lim_upper = str2double(time_lim{1,1});
limbl_low = str2double(time_lim{3,1});
limbl_upper = str2double(time_lim{2,1});

%%
for counter = 1:length(Allsets)

   
   [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',counter,'study',0);   % Retrieve the counter dataset. 
   EEG = eeg_checkset( EEG );
   eeglab redraw
   
   allconds = unique({EEG.event.type});
   conds2seg = allconds(ismember(allconds,Condsoi));   %Find the conditions for the current dataset.
   % Find the indices of the conds to segment.
   
   
   %% Check the name of the video corresponding to the current dataset and
   %check if it is necessary to create a new folder for this video name or
   %if such a folder already exists.
   
   dcurrs = dir(dirsave);
   alldirs = {dcurrs([dcurrs.isdir]).name};
   ishere = find(ismember(alldirs,{EEG.video_name}));
   
   if isempty(ishere)
       
      [status, msg, msgID] = mkdir(fullfile(dirsave,EEG.video_name));
      [status1, msg1, msgID1] = mkdir(fullfile(dirsave,EEG.video_name, 'BLCorrected'));
       
   else
       
       disp(['********The folder ,', EEG.video_name,',  already exists in current directory*******']);
       
   end
   
   savedir = fullfile(dirsave,EEG.video_name);
   savedirBL = fullfile(dirsave,EEG.video_name, 'BLCorrected');
   
   
   %% Carry out the segmentation
   
   for icond = 1:length(conds2seg)
       
       [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',counter,'study',0);   % Retrieve the counter dataset. 
       EEG = eeg_checkset( EEG );
       eeglab redraw
       
       % find the indices of the current condition to segment.
       indx1 = ismember({EEG.event.type},conds2seg{1,icond});
       if isfield(EEG.event,'feedback')
            FBinfo = {EEG.event(indx1).feedback};  % Isolate the exact feedbacks for the current condition.
            isfbfield = 1;
       else
           disp('*********No Feedback sub-field in current event field******************');
           isfbfield = 0;
       end
       TrialCnt = [EEG.event(indx1).trialnum];
       
       disp('--------------------Segment continuous data-----------------------------------');
       newepoch_title = [currsuj,'_',conds2seg{1,icond}];
       [EEG, epindx] = pop_epoch(EEG, {conds2seg{1,icond}}, [limbl_low lim_upper], 'newname', newepoch_title, 'epochinfo', 'yes','eventindices',find(indx1));
       [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',newepoch_title,'gui','off');
       EEG = eeg_checkset( EEG );
       eeglab redraw;
       
       % Assign the feedback and trialnum details to the events
       % structure of the current segmented EEG dataset.
       X ={EEG.epoch.event};
       condcurrind = zeros(length(X),1);
       for xcnt = 1:length(X)
           condcurrind(xcnt) = X{1,xcnt}(1);
           
       end
       
       for epcnt = 1:length(condcurrind) 
           
           if isfbfield ==1
            EEG.event(condcurrind(epcnt)).feedback = FBinfo{1,epcnt};
           end
           EEG.event(condcurrind(epcnt)).trialnum = TrialCnt(epcnt);
       end
       

       EEG = pop_saveset( EEG, 'filename',newepoch_title,'filepath',savedir);
       EEG = eeg_checkset( EEG );
       eeglab redraw;
       
       disp('--------------------Baseline correction-----------------------------------');
       baseline_low = limbl_low*1000;  % Change here to change baseline to use for correction.
       baseline_hi = limbl_upper*1000;
       
       Enom_bl=strcat(newepoch_title,'-bl');
       EEG = pop_rmbase( EEG, [baseline_low baseline_hi]);
       [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(Enom_bl),'gui','off');
       EEG = eeg_checkset( EEG );
       
       EEG = pop_saveset( EEG, 'filename',char(Enom_bl),'filepath',savedirBL);
       EEG = eeg_checkset( EEG );
       eeglab redraw
   
   end
   
end

disp('***************************All Done - Suckin Diesel****************************');
