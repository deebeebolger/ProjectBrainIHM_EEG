

%% INTEGRATE THE NEW TRIGGERS TIMES INTO THE DATA
% Need to load in the following excel files:
% 1. PEPs_sessions_lists.xlsx : which lists correspond to which subject
% 2. Brain_IHM_lists.xlsx     : the names of the videos in each list

currsuj = 's01';

sujs_lists = 'PEPs_sessions_lists.xlsx';
lists_vids = 'Brain_IHM_Lists.xlsx';
currdir = fullfile(filesep,'Users','bolger','Documents','work','Projects','Project-BrainIHM');

X  = readtable(fullfile(currdir,sujs_lists),'FileType','spreadsheet','Sheet',2);                                         % Only interested in columns 1 and 2

% Search of name of list for the currsuj
sujnums   = X{:,1};
listnames = X{:,2};
ix        = [sujnums == str2double(currsuj(2:end))];
currlist  = listnames{ix,1};

% Find the names of the videos included in the list corresponding to the
% current suject. 

X1 = readtable(fullfile(currdir,lists_vids),'VariableNamingRule','preserve','FileType','spreadsheet','Sheet',currlist);  % Only interested in columns 1 and 2



%% Load in the annotations textfile for the first film

currfilms = X1{:,2};  % All the films for the current list. 

counter = 1;
currfilm = [currsuj,'-Film',num2str(counter)];  % Title of the subject sub-folder

currfile = currfilms{counter,1};
dataIn_file = fullfile(filesep,'Users','bolger','Desktop','Annotations-Philippe-excel',[currfile,'.txt']);
DataIn = readtable(dataIn_file);

% Sometimes the second column is NaN, so test for this just in case.
% If second column is NaN, then column 3 is start, column 4 is end and
% column 5 is the duration in milliseconds.

if isnan(DataIn{:,2})

    startms  = DataIn{:,3};
    endms    = DataIn{:,4};
    durms    = DataIn{:,5};
    feedback = DataIn{:,6};   % this should be a string.
else

    startms  = DataIn{:,2};
    endms    = DataIn{:,3};
    durms    = DataIn{:,4};
    feedback = DataIn{:,5};   % this should be a string.
end

% Isolate the feedbacks of the patient.

pat_fb      = cell2mat(cellfun(@isempty, feedback, 'UniformOutput',false));
startms_pat = startms(~pat_fb);

%% Load in the continuous data for the current subject and film.
%  The data should be resampled, rereferenced, filtered and the bad
%  channels should be removed.

datadir = fullfile(filesep,'Users','bolger','Documents','work','Projects','Project-BrainIHM','Data_for_trigger_correct');
datadir_curr = fullfile(filesep,datadir,currsuj,currfilm,filesep);

sfiles = dir(strcat(datadir_curr,'*.set'));
findsets=find(~[sfiles.isdir]);        
Allsets = {sfiles(findsets).name}; 

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename',Allsets,'filepath',datadir_curr);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
eeglab redraw;


%% We are only interested in the verbal cong (and incong).
%  We are not interested in the gestual as the new annotations take into
%  account both gestual and auditory dimensions of the feedbacks.

icong        = contains({EEG.event.type},'Verbal-cong');
fbtimes      = [EEG.event.etimes];
fbtimes_cong = fbtimes(icong);
lats         = [EEG.event.latency];
lats_cong    = lats(icong);


if contains(currfile,'_Incongruent')

    iincong        = contains({EEG.event.type},'Verbal-incong');
    fbtimes_incong = fbtimees(iincong);
    lats_incong    = lats(iicong);

else

    iincong        = [];
    fbtimes_incong = [];
end

%% Compare the new annotation onset times with the original onset times. 
%  Look at the congruent feedbacks to begin with.
%  Take into account the Inter-onset difference. 
%  The length of the two ITI variables should match. 
%  Write to file the original, new start times (ms) and the difference
%  between the two. 

if ~isempty(strfind(currfile,'Agent'))
    onsets_fname = 'Agent_Onsets_Original.xlsx';
elseif ~isempty(strfind(currfile,'Human'))
    onsets_fname = 'Human_Onsets_Original.xlsx';
end

XIn = readtable(fullfile(filesep,'Users','bolger','Desktop','Annotations-Philippe-excel',onsets_fname),'FileType','spreadsheet','Sheet',currfile);

diff_onsets = startms_pat - XIn.Onsets*1000;
onset_comp = [XIn.Onsets*1000,startms_pat, diff_onsets];   % Compare onsets of the original and newly annotated onsets. 

Tfb = table(XIn.Onsets*1000,startms_pat, diff_onsets);
curr_itifile = [currfile,'-onset_comp.csv'];

iti_fpath = fullfile(filesep,'Users','bolger','Desktop','Annotations-Philippe-excel',curr_itifile);
writetable(Tfb,iti_fpath,'Delimiter',';')  

% Find the ITI differences that exceed 0.5ms and calculate the new start
% feeback onset for these differences.

Y = [abs(new_orig_iti)>0.5];

fbtimes_cong(Y)





