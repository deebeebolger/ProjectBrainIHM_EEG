
audnom = 'Reg_Human_Incongruent_RS'; 
auddir = fullfile(filesep,'Volumes','deepassport','Projects','Project-PEPs','PEPS-protocol-phase2','Audio-Files',filesep);
audfiles = dir(strcat(auddir,'*.wav'));
findwav=find(~[audfiles.isdir]);        
Allwavs= {audfiles(findwav).name}; 
X = contains(string(Allwavs),audnom);  % Search auddir for the file with name matching that of audnom.

txtfiles = dir(strcat(auddir,'*.txt'));
findtxt=find(~[txtfiles.isdir]);        
Alltxt= {txtfiles(findtxt).name}; 
Xtxt = contains(string(Alltxt),audnom);

diraud_full = fullfile(auddir,Allwavs{X});
dirtxt_full = fullfile(auddir,Alltxt{Xtxt});
audtxtIn = readtable(dirtxt_full);

% Read in the *.wav file (audio file)
[Y, sF] = audioread(diraud_full);
info = audioinfo(diraud_full);
time = 0:1/sF:(1/sF)*length(Y);    % Create the time vector.
audtime = time; 




%% Mark the trigger time points on the audio time vector 

trigindx = zeros(length(audtrigs),1);
trigindx1 = zeros(length(audtrigs),1);
trigtime_onset = nan(length(audtime),1);
audtrigs_onset = audtxtIn{:,1};
audtrigs_offset = audtxtIn{:,2}; 
trigtime_offset = nan(length(audtime),1);

for fcnt = 1:length(trigindx)  
    trigindx(fcnt) = dsearchn(audtime',audtrigs_onset(fcnt));
    trigtime_onset(trigindx(fcnt)) = .01;
    
    trigindx1(fcnt) = dsearchn(audtime',audtrigs_offset(fcnt));
    trigtime_offset(trigindx1(fcnt)) = .01;
    
    
    
end


%%

bbandenv = PEPs_EnvelopeCalc(mean(Y,2), sF, audtime(1:size(Y,1)), audnom, trigtime_onset(1:size(Y,1)));

figure;
plot(audtime(1:length(bbandenv)), bbandenv)
hold on
stem(audtime,trigtime_onset,'r')
hold on
stem(audtime,trigtime_offset, 'b')

