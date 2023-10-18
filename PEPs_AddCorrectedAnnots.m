%% Date: 25/4/2022                   Programmed: Deirdre Bolger
% Script to take in the newly corrected feedback onsets for the
% Human-Machine videos and integrate them into the EEG datasets.
% As these new onsets do not correspond to the onsets extracted from the
% photodiode signal, these new onsets are integrated by corrected the
% already defined feedback onsets in the EEG datasets.
%**************************************************************************


%% INTEGRATE THE NEW TRIGGERS TIMES INTO THE DATA
% Need to load in the following excel files:
% 1. SubjectLists_Summary.xlsx : which lists correspond to which subject
% The first sheet of the excel file informs the list numbers for each
% subjects.
% The other 12 sheets, informs the Films that comprise each list number.
% There is a sheet per List number.


currsuj = 'S25';
fcnt    = 2;

sujs_lists = 'SubjectLists_Summary.xlsx';
currdir = fullfile(filesep,'Users','bolger','Documents','work','Projects','Project-BrainIHM');

X  = readtable(fullfile(currdir,sujs_lists),'FileType','spreadsheet','Sheet','subject-lists');

CurrList = X{ismember(X{:,1},currsuj),2};      % Find the list corresponding to the current participant.

Listdata = readtable(fullfile(currdir,sujs_lists),'FileType','spreadsheet','Sheet',CurrList{1,1});   % Import data from current list sheet.

Filmnums = Listdata.Film_Number;
Filmtype = Listdata.Film_Type;

%% DEFINE THE PATH TO THE TEXT FILE WITH THE NEWLY DEFINED FEEDBACK ONSETS AND ANNOTATIONS

newannot_path = '/Users/bolger/Documents/work/Projects/Project-BrainIHM/PEPs_AnnotationCompare1.xlsx';

%% LOAD IN THE DATA FOR THE CURRENT PARTICIPANT.
% Need to load in the dataset with the patient data

datapath = '/Users/bolger/Documents/work/Projects/Project-BrainIHM/Data_for_trigger_correct/';


% We will begin by correcting only the human videos
hindx = find(contains(Filmtype,'Human'));     % Find the index of human-as patient films.
hfilms = {Filmtype{hindx,:}};
nfilms = {Filmnums{hindx,:}};

% Start an eeglab session.
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% For the current film carry out onset correction. 

currfilm = strcat(currsuj,'-',nfilms{1,fcnt});
currpath = fullfile(datapath,currsuj,currfilm,filesep);
sfiles = dir(strcat(currpath,'*Trigcorr1.set'));              % Load in the dataset with the patient feedback onsets first.

% Load the EEG data.
EEG = pop_loadset('filename',sfiles.name,'filepath',currpath);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
eeglab redraw;

%% Load in the text file with the newly defined feedback onsets.
A        = readtable(newannot_path,'FileType','spreadsheet','Sheet',hfilms{1,fcnt});
notnan   = ~isnan(A.Version1_Version2_onset);
currdiff = A.Version1_Version2_onset;
currdiff = currdiff(notnan);                           % Difference between old and new feedback onsets (without the NaNs).


curronset_old = A.Version1;
curronset_old = curronset_old(~isnan(curronset_old));
curronset     = A.Version2_onset;                      % New annotation onsets (both strokes and preparation strokes)
currfbks_all  = A.Version2_annot;
currdur_all   = A.Version2_dur;
currcond_all  = A.Version2_cond;

notnan2       = cell2mat(cellfun(@isempty, currfbks_all, 'UniformOutput',false));
currfbks      = currfbks_all(~notnan2);
currstroke    = A.Version2_stroketype;
currstroke_v  = currstroke(~notnan2);
iprep         = contains(currstroke, 'preparation');
iprepv        = contains(currstroke_v, 'preparation');  % Index of the preparation strokes out of those strokes corresponding to feedback onsets.
curronset_fb  = curronset(~notnan2);                    % The onsets that corresponds to the feedback onsets.
currdur       = currdur_all(~notnan2);
currcond      = currcond_all(~notnan2);

eegtimes  = [EEG.event.etimes]';                        % Extract the onset times of the current EEG dataset (in seconds).
eegtype   = {EEG.event.type}';
eegfbks   = {EEG.event.feedbacks}';
iverb     = contains(eegtype, 'Verbal');
eegfbv    = eegfbks(iverb);
eegtimesv = eegtimes(iverb);


%% Find the difference between the preparation stroke and the following
% feedback onset.
prepcnt = 1;
currcnt = 1;
eegnew  = zeros(length(curronset_fb),1);
eeglat  = zeros(length(curronset_fb),1);
eegtype = cell(length(curronset_fb), 1);
eegtnum = zeros(length(curronset_fb),1);  % trial numbers
prepv_idx = find(iprepv);

for tcnt = 1:length(eegtimesv)


    if currcnt == prepv_idx(prepcnt)

        disp('preparation-stroke')

        eegnew(currcnt)    = eegtimesv(tcnt) + (currdiff(currcnt)*-1);
        eeglat(currcnt)    = eegnew(currcnt+1)*EEG.srate;
        eegtype{currcnt,1} = 'preparation';
        eegtnum(currcnt)   = tcnt;         % Assign the preparation stroke to the same trial as the following feedback.

        currcnt = currcnt +1;
        eegnew(currcnt) = eegtimesv(tcnt) + (currdiff(currcnt)*-1);
        eeglat(currcnt) = eegnew(currcnt)*EEG.srate;
        eegtype{currcnt,1} = 'Verbal-feedback';
        eegtnum(currcnt)   = tcnt;          % Assign the feedback to the same trial as the preceding preparation stroke.

        if prepcnt < length(prepv_idx)
            prepcnt = prepcnt +1;
        end

    else

        disp('not a preparation stroke')
        eegnew(currcnt) = eegtimesv(tcnt) + (currdiff(currcnt)*-1);
        eeglat(currcnt) = eegnew(currcnt)*EEG.srate;
        eegtype{currcnt,1} = 'Verbal-feedback';
        eegtnum(currcnt)   = tcnt;
    end

    currcnt = currcnt +1;

end

%% Need to construct the event field of the EEG structure by adding the following:
%  - latency
%  - type (verbal-feedback ou preparation)
%  - etimes (in seconds)
%  - feedback annotations.
%  - duration
%  - trial numbers
%  - urevent

EEG.event = [];
EEG.urevent = [];

for ecnt = 1:length(eegnew)

    EEG.event(ecnt).latency   = eeglat(ecnt);
    EEG.event(ecnt).type      = currcond{ecnt,1};
    EEG.event(ecnt).fbtype    = eegtype(ecnt);
    EEG.event(ecnt).duration  = currdur(ecnt);
    EEG.event(ecnt).urevent   = ecnt;
    EEG.event(ecnt).trialnum  = eegtnum(ecnt);
    EEG.event(ecnt).etimes    = eegnew(ecnt);
    EEG.event(ecnt).feedbacks = currfbks{ecnt,1};

    EEG.urevent(ecnt).type    = eegtype(ecnt);
    EEG.urevent(ecnt).latency = eeglat(ecnt);

end

%% Add the latency in samples for the added preparation strokes.

 i=find([EEG.event.latency]==0)';

 for icnt = 1:length(i)

     EEG.event(i(icnt)).latency = EEG.event(i(icnt)).etimes*EEG.srate;
     EEG.urevent(i(icnt)).latency = EEG.event(i(icnt)).etimes*EEG.srate;

 end

 %% Save to a new dataset.


newfnom = [EEG.setname,'-Trigcorr2'];
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(newfnom),'gui','off');
EEG = eeg_checkset( EEG );

EEG = pop_saveset( EEG, 'filename',newfnom,'filepath',currpath);  % Saves a copy of the current resampled dataset to the current directory
eeglab redraw





