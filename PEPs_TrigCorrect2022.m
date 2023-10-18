
% Load in the excel file with trigger information for the current
% participant.

curr_suj = 'S25';
curr_film = 'S25-Film4';
curr_vid  = 'Curare-Human-Incongruent';
compfile  = 'S25_CompareTriggers_May2022.xlsx';


Dirbase = fullfile(filesep,'Users','bolger','Documents','work','Projects','Project-BrainIHM','Data_for_trigger_correct',curr_suj,filesep);
Savebase = fullfile(filesep,'Users','bolger','Documents','work','Projects','Project-BrainIHM','Data_for_trigger_correct',curr_suj,curr_film,filesep);

trigcheck_file = fullfile(filesep, Dirbase,compfile);
TrigData = readtable(trigcheck_file,'FileType','spreadsheet','Sheet',curr_vid);  

datapath = fullfile(filesep,'Users','bolger','Documents','work','Projects','Project-BrainIHM','Data_for_trigger_correct',curr_suj, curr_film,filesep);

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);


sfiles = dir(strcat(datapath,['*',curr_film(end-4:end),'.set']));              % Load in the dataset with the patient feedback onsets first.

% Load the EEG data.
EEG = pop_loadset('filename',sfiles.name,'filepath',datapath);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
eeglab redraw;

if contains(curr_vid,'Congruent')

    Type     = TrigData.Type;
    FBtype   = TrigData.Feedbacks;
    Annots   = TrigData.Annotations;
    Lats     = TrigData.Latency;
    Trialnum = TrigData.Trialnum;
    Onsets   = TrigData.OriginalOnsets;

%     Type     = TrigData.Var1;
%     FBtype   = TrigData.Var2;
%     Annots   = TrigData.Var3;
%     Lats     = TrigData.Var4;
%     Trialnum = TrigData.Var5;
%     Onsets   = TrigData.Var7;

    EEG.event   = [];
    EEG.urevent = [];

    for ecnt = 1:length(Type)

        EEG.event(ecnt).type     = [FBtype{ecnt,1},'-',Type{ecnt,1}];
        EEG.event(ecnt).latency  = Lats(ecnt);
        EEG.event(ecnt).urevent  = ecnt;
        EEG.event(ecnt).trialnum = Trialnum(ecnt);
        EEG.event(ecnt).etimes   = Onsets(ecnt);
        EEG.event(ecnt).feedbacks = Annots(ecnt);

        EEG.urevent(ecnt).type = Type{ecnt,1};
        EEG.urevent(ecnt).latency = Lats(ecnt);

    end

elseif contains(curr_vid,'Incongruent')

    Type     = TrigData.Type;
    FBtype   = TrigData.Feedbacks;
    Annots   = TrigData.Annotations;
    Lats     = TrigData.Latency;
    Trialnum = TrigData.Trialnum;
    Onsets   = TrigData.OriginalOnsets;

%     Type     = TrigData.Var1;
%     FBtype   = TrigData.Var2;
%     Annots   = TrigData.Var3;
%     Lats     = TrigData.Var4;
%     Trialnum = TrigData.Var5;
%     Onsets   = TrigData.Var7;

    EEG.event   = [];
    EEG.urevent = [];

    for ecnt = 1:length(Type)
        

        EEG.event(ecnt).type     = [FBtype{ecnt,1},'-',Type{ecnt,1}];
        EEG.event(ecnt).latency  = Lats(ecnt);
        EEG.event(ecnt).urevent  = ecnt;
        EEG.event(ecnt).trialnum = Trialnum(ecnt);
        EEG.event(ecnt).etimes   = Onsets(ecnt);
        EEG.event(ecnt).feedbacks = Annots(ecnt);

        EEG.urevent(ecnt).type = Type{ecnt,1};
        EEG.urevent(ecnt).latency = Lats(ecnt);

    end

end

newfnom = [EEG.setname,'-Trigcorr1'];
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(newfnom),'gui','off');
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(newfnom),'filepath',Savebase);  % Saves a copy of the current resampled dataset to the current directory
eeglab redraw











