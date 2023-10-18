
%% PEPS_EEGPreprocessing_MAIN         Programmed: D. Bolger
% Specific pre-processing script for PEPS project. 
% Script to carry out the different steps in the pre-processing of EEG
% data.
% It applies functions from EEGLAB, the PREP pipeline.
% Functions from the ADJUST toolbox are called in the "CREx_ICA_calc()".
% The following functions are called:
% Needs to be filled in...

%% OPEN THE PARAMETERS TEXTFILE. 
% The parameter textfile contains the necessary pre-processing parameters as well as the paths.
% 
close all
clear all

paramfile_nom = 'parameters-PEPS.txt';  % The title of the parameters file. 
% The path to the parameters path. This should be the only path that needs
% to be changed.

paramfile_path = fullfile(filesep,'Users','bolger','Brain-IHM','Data_Preproc',paramfile_nom); % Put the path to your parameters file here. 

fid2 = fopen(paramfile_path);  % Define a file identifier.
mydata = textscan(fid2,'%s %s');  % Scan textfile...

for i = 1:length(mydata{1,1})                     % generate a parameters structure from the parameters text file
    Params.(genvarname(mydata{1,1}{i})) = mydata{1,2}(i);
end

%% OPEN AN EEGLAB SESSION AND LOAD MANUALLY THE FILE TO BE PROCESSED
suj = 's40';
sujnom = 's40-Film4';
sujdir = fullfile(Params.Datadir{1,1},suj,sujnom,filesep);
currfile = dir(strcat(sujdir,'*doctrigs.set'));
currtitre = currfile.name;
currfnom = currfile.name(1:end-4);


[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename',currtitre,'filepath',sujdir);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
eeglab redraw;

%% ADD CHANNEL INFORMATION TO THE CURRENT DATASET.
% Channel coordinates and labels are added to the current dataset and
% the dataset is saved to the current subject-level directory.
% The Chaninfo.mat file is loaded as it contains the electrode labels.
% From EEGLAB plugins, the file, "standard-10-5-cap385.elp" is loaded
% as this contains the correct coordinates for the 10-20 system used here.
chanloc_path = Params.Chanlocdir{1,1};
DIRsave_curr = fullfile(Params.Savedir{1,1},suj,sujnom,filesep);

fnom_chans = strcat(currfnom,'-chan');
EEG = pop_chanedit(EEG, 'lookup',chanloc_path);                % Load channel path information
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(fnom_chans),'filepath',DIRsave_curr);
eeglab redraw

%% PREPARE INFORMATION TEXT-FILE FOR THE CURRENT SUBJECT.
% This textfile should be saved in each subject folder.
% Need to find which Film this dataset corresponds to...
ifilm = strfind(currtitre,'Film');
currfilm = currtitre(ifilm:ifilm+length('Film'));
fname = strcat(sujnom,'-info.txt');
fdir = strcat(DIRsave_curr,fname);
fid = fopen(fdir,'w');
fprintf(fid,['---------',sujnom,'-',currfilm,'----------\n']);

%% RESAMPLE DATA TO THE SAMPLING RATE DEFINED IN THE PARAMETER FILE (IF REQUIRED)
% Resamples using the EEG resampling function.
% If the user has the Matlab signal processing toolbox, it uses the
% Matlab resample() function.
% Write information regarding the resampling to the subject-level text
% file.

SR_orig = EEG.srate;
SR_new = str2double(Params.srate{1,1});

if SR_orig>SR_new

    fprintf(fid,'\nDownsampled from %dHz to %dHz\n',SR_orig,SR_new);
    disp('***********************************Resampling to 512Hz*******************************************')
    fnom_rs = strcat(fnom_chans,'-rs');

    EEG = pop_resample(EEG, SR_new);   %resample the data at sampling rate defined, sr.
    EEG =eeg_checkset(EEG);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom_rs),'gui','off'); % current set = xx;
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(fnom_rs),'filepath',DIRsave_curr);  % Saves a copy of the current resampled dataset to the current directory
    eeglab redraw
    
elseif SR_orig==SR_new
    fprintf(fid,'\nSampling rate already at %dHz. No need to resample.\n',SR_new);
    disp('***********************************Resampling can be skipped*******************************************');
    fnom_rs = fnom_chans;
end

%% IT IS NECESSARY TO KNOW THE INTER-TRIAL INTERVAL FOR THE VERBAL FEEDBACKS.
% For time-frequency decomposition, in particular, a large epoch duration
% is required to the resolve the low-frequency oscillations.
% The actual latencies of each event in seconds is calculated for all
% events (gestual and verbal). 
% The latencies in seconds of the verbal events is extracted.
% The inter-trial interval (in seconds) for the verbal feedbacks is
% calculated and saved as a field in the EEG structure of the current
% subject. 

event_times = [EEG.event.latency]./512;
event_times = event_times';

for ecnt = 1:length(EEG.event)
    EEG.event(ecnt).etimes = event_times(ecnt);
end

X = ismember({EEG.event.type},'Fb verbal');  % Find the latencies of the verbal events only.
verb_times = event_times(X);
verbtimes_diff = diff(verb_times);  % Time intervals between the verbal feedbacks.

EEG.vtimes_diff = verbtimes_diff;

[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
eeglab redraw;

%% LOAD IN ASSOCIATED AUDIO FILE AND RESAMPLE TO EEG SAMPLING RATE.
% 
% auddir = Params.Audiodir{1,1};           % The filepath to the audio files.
% audfiles = dir(strcat(auddir,'*.wav'));  % Find the wavfiles
% audfile_nom = {audfiles.name};
% fm = strfind(audfile_nom, EEG.video_name);
% iaud = [cell2mat(cellfun(@isempty,fm ,'UniformOutput',false))==0];
% audcurr = audfile_nom{1,iaud};
% 
% txtfiles = dir(strcat(auddir, '*.txt')); % Find the wavfiles.
% txtfiles_nom = {txtfiles.name};
% tfm = strfind(txtfiles_nom, EEG.video_name); 
% itxt = [cell2mat(cellfun(@isempty, tfm, 'UniformOutput',false))==0];
% txtcurr = txtfiles_nom{1, itxt};
% 
% audfulldir = fullfile(auddir,audcurr);  %Full path to auditory file
% txtfulldir = fullfile(auddir, txtcurr); %Full path to the txtfile with FB marker info.
% 
% [Y, fsamp] = audioread(audfulldir);      % Read in audio file (in *.wav format). 
% Y = mean(Y,2);
% time = 0:1/fsamp:(1/fsamp)*length(Y);    % Calculate the time vector of audio
% 
% % Need to align the EEG and the auditory stimulus.
% % EEG sampling rate = 512; audio sampling rate = 44100
% 
% Fs_new = 512;
% [Numer, Denom] = rat(Fs_new/fsamp);
% Ynew = resample(Y, Numer, Denom);
% tnew = 0:1/Fs_new:(1/Fs_new)*length(Ynew);   %New time vector
% tnew = tnew(1:end-1);

%% Visualise the auditory signal.
% fbinfo = readtable(txtfulldir,'ReadVariableNames',false);
% fbindx = zeros(length(fbinfo{:,1}),1);
% fbtime = nan(length(tnew),1);
% 
% for fcnt = 1:length(fbindx)
%     
%     fbindx(fcnt) = dsearchn(tnew',fbinfo{fcnt,1});
%     fbtime(fbindx(fcnt)) = 0.1;
%     
% end
% 
% %[Axes, H1, H2] = plotyy(tnew,Ynew,tnew,fbtime','plot','stem');
% plot(tnew, Ynew)
% hold on
% stem(tnew,fbtime,'r')
% scrollplot(gca,'WindowSizeX',5, 'MinX',0)

%% APPLY BAND-PASS FILTER BETWEEN THE LOWER AND UPPER LIMITS SPECIFIED IN PARAMETERS FILE.
% It applies a FIR windowed sinc filter using a blackman filter.
% The filter frequency response is plotted.
% The details of the filtering are written to subject information txt file.

f_low = str2double(Params.fc_low{1,1});
f_hi = str2double(Params.fc_hi{1,1});
fnom_filt = strcat(fnom_rs, '-filt');

disp('*********************************Bandpass filtering using a FIR windowed sinc filter***********************************')
[M, dev] = pop_firwsord('blackman',EEG.srate, 2);
[EEG,com,b] = pop_firws(EEG,'fcutoff',[f_low f_hi],'forder',M,'ftype','bandpass','wtype','blackman');
fvtool(b);                                      % Visualise the filter characteristics
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom_filt),'gui','off');   %save the resampled data as a newdata set.
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(fnom_filt),'filepath',DIRsave_curr);
eeglab redraw

fprintf(fid,'Band-pass filtered %f2 - %dHz with %d-order fir windowed sinc filter (blackman).\n',f_low,f_hi,M);

%% VISUALISE THE REFERENCE AND THEIR SPECTRA CALCULATED USING MULTI-TAPERS.
% Saves a figure of the spectra of the references to the current
% folder (in the CREx_SpectCalc_multitap() function.

if strcmp(Params.references{1,1},'average')
    Refs = Params.references{1,1}; 
    refs = [];
    sprintf('A mean global reference has been defined.\n Need to visual spectra of all electrodes')
else

    refs = [str2double(Params.references{1,1}(1:2)) str2double(Params.references{1,1}(3:4))];
    eegplot(EEG.data(refs,:),'srate',EEG.srate,'eloc_file',EEG.chanlocs([refs(1) refs(2)]),'events',EEG.event,'color',{'g' 'b'},'dispchans',2,...
        'winlength',20,'title','Visualise reference electrodes (10sec per window)');

    specnom_ref = fullfile(DIRsave_curr,strcat(fnom_filt,'-spectref'));
    CREx_SpectCalc_multitap(EEG,refs,[1 60],specnom_ref,.1);
end



%% VISUALISE THE SPECTRA OF ALL SCALP ELECTRODES BEFORE RE-REFERENCING.
% The spectrum is saved as a *.fig file. 

specnom_scalp = fullfile(DIRsave_curr,strcat(fnom_filt,'-spectscalp'));

CREx_SpectCalc_multitap(EEG,1:64,[1 60],specnom_scalp,.05);

%% RE-REFERENCE THE DATA TO THE ELECTRODES SPECIFIED IN THE PARAMETERS FILE.
% The channels used for referencing are generally EXG1 and EXG2,
% channels 65 and 66, respectively.
% The details of the re-referencing are written to the information text file.

disp('***********************Rereference to Defined Channel:  does zero potential exist?*****************************')
EEG = pop_reref(EEG, {'EXG1' 'EXG2'}, 'method','standard','keepref','on');
fnom_ref = strcat(fnom_filt,'-rref');
EEG = eeg_checkset( EEG );
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(fnom_ref),'gui','off');   %save the resampled data as a newdata set.
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(fnom_ref),'filepath',DIRsave_curr);
EEG = eeg_checkset( EEG );
eeglab redraw

here = CURRENTSET;   % Mark the current set.

if isempty(refs)
    fprintf(fid,'Rereferenced using mean global reference\n');
else  
    fprintf(fid,'Rereferenced using channels %s and %s.\n\n',EEG.chanlocs(refs(1)).labels,EEG.chanlocs(refs(2)).labels);

end

%%