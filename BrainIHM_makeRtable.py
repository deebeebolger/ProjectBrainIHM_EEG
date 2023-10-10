import numpy as np
import mne
import matplotlib.pyplot as plt
import glob
import os
import csv
import pandas as pd

def load_channels():
    currwd = os.getcwd()
    channel_path = os.path.join(currwd, 'Channels.csv')
    with open(channel_path, newline='') as f:
        readin = csv.reader(f)
        Channels = list(readin)
        return Channels

## This script creates a data table compatible with R (long format) with the following columns:
#  Subject/ Channel/Group/Condition/MV

dir_base = '/Users/bolger/Documents/work/Projects/project_Brain-IHM/'   # This is the base directory
Groups = ['Data_Human', 'Data_Agent']                                   # Define the groups (note will have to strip 'Data_')

# Get the subject names
g1_path = os.path.join(dir_base, Groups[0], 'EpochData')
g1_contents = os.listdir(g1_path)
sujfiles = [suj for suj in g1_contents if suj[0] == 'S']   # Find the files that correspond to subject files (begin with S).
sujfiles.sort(key= lambda x: float(x.strip('S')))          # Sort the files names according to subject number.

# Define the Group and Condition names. These will be used in the R table.
Group_names = [s.replace('Data_', '') for s in Groups]      # Just keep "Human" and "Agent"
Cond_names  = ['-Congru-Congru', '-InCongru-Congru', '-InCongru-InCongru']
AllChans    = load_channels()  # Call of function to load in channels.

# Initialize the data columns.
subjects  = []
electrodz = []
groupnom  = []
condnom   = []
data_mv   = []



for sujindx, sujcurr in enumerate(sujfiles):
    print(f'The the current participant is {sujcurr}')

    for chanidx, chancurr in enumerate(AllChans):
        print(f'The the current channel is {chancurr[0]}')

        for gpindx, gpcurr in enumerate(Groups):

            gpath_curr = os.path.join(dir_base, gpcurr, 'EpochData')
            spath_curr = os.path.join(gpath_curr, sujcurr)
            spath_contents = os.listdir(spath_curr)
            meanfiles = [sfile for sfile in spath_contents if sfile.endswith('ave.fif')]

            for condindx, condcurr in enumerate(Cond_names):

                condfile_curr = [c for c in meanfiles if condcurr in c]

                if condfile_curr == []:
                    print(f'The condition {condcurr} does not exist for the current subject, {sujcurr}, of group, {Group_names[gpindx]}.\n')
                else:
                    print(f'The condition {condcurr} exists for the current subject, {sujcurr}, of group, {Group_names[gpindx]}.\n')


                    avg2load_path = os.path.join(spath_curr, condfile_curr[0])
                    evokedIn = mne.read_evokeds(avg2load_path, baseline=None, kind='average', verbose=False)   # Load in the current evoked data, averaged over trials.
                    Evoked      = evokedIn[0]
                    time_vector = Evoked.times
                    evokedData  = Evoked.get_data()
                    chans_curr  = Evoked.info['ch_names']
                    cindx = chans_curr.index(chancurr[0])  # Get the index of the current channel
                    mv_chancurr = evokedData[cindx,:]

                    subjects.append(sujcurr)
                    electrodz.append(chancurr[0])
                    groupnom.append(Group_names[gpindx])
                    condnom.append(condcurr[1:])


                    if data_mv == []:
                        data_mv = mv_chancurr
                    else:
                        data_mv = np.dstack((data_mv, mv_chancurr))


Subject_col = np.transpose(subjects)
Channel_col = electrodz
Group_col   = np.transpose(groupnom)
Condition_col = np.transpose(condnom)
MV_col = data_mv[0].transpose()

# Intialize text file in which to write data.
Rfile_name = 'AllData_Rformat.csv'
Rfile_path = os.path.join(os.getcwd(), Rfile_name)

Big_array = [Subject_col, Channel_col, Group_col, Condition_col]
BigDF     = pd.DataFrame(Big_array).T
len(BigDF.columns)

s = np.shape(MV_col)
data2append = {}
for i in range(s[1]):
    data2append[i] = MV_col[:,i]
    data2append_s = pd.Series(data2append[i])
    BigDF = pd.concat([BigDF, data2append_s], axis=1)

BigDF.to_csv(Rfile_path)

