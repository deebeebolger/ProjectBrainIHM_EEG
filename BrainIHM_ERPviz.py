import numpy as np
import mne
import matplotlib.pyplot as plt
import glob
import os
import matplotlib
from matplotlib.patches import FancyArrowPatch, ArrowStyle, Rectangle
import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
from matplotlib.widgets import SpanSelector
import math
import pandas as pd


def cond_collapse():
    """
    Function to collapse the two types of congruent conditions (Cong-Cong, InCong-Cong) for those
    participants having both congruous conditions.
    :param
    :return:
    """


def sem_calc(evokedIn):
    """
    Function to calculate the standard error of the mean for each channel taking into account all subjects.
    :param evokedIn: the evoked data for each subject
    :return: the upper and lower limit of sem to plot
    """
    evshape = np.shape(evokedIn)
    sem_sujs = np.std(evokedIn, 2, ddof=1) / np.sqrt(evshape[2])
    sujmean = np.mean(evokedIn, 2)
    upperlim = sujmean + sem_sujs
    lowlim = sujmean - sem_sujs

    return upperlim, lowlim


def gfp_calc(evoked_data, t):
    dataev1 = evoked_data[:, :, 0]
    dataev2 = evoked_data[:, :, 1]
    dataev1_mean = np.mean(dataev1, axis=0) ** 2
    dataev2_mean = np.mean(dataev2, axis=0) ** 2
    dataev1_rm = []
    dataev2_rm = []
    for i in range(0, len(t)):
        dev1_rm = np.sum((dataev1[:, i] - dataev1_mean[i])) / 64
        dev2_rm = np.sum((dataev2[:, i] - dataev2_mean[i])) / 64
        dataev1_rm.append(dev1_rm)
        dataev2_rm.append(dev2_rm)

    gfp1 = np.sqrt(np.abs(dataev1_rm))
    gfp2 = np.sqrt(np.abs(dataev2_rm))

    return gfp1, gfp2

def stimLen_Calc(groupnom, condname, slenDir):
    if condname == 'CongruentAll':
        condname1 = 'Congru-Congru'
    elif condname == 'InCongru-Congru':
        condname1 = 'InCongru-InCongru'
    else:
        condname1 = condname

    fnomPart = groupnom + '-' + condname1 + '-markerchannel.csv'
    fnomPart_dir = os.path.join(slenDir, fnomPart)
    print(fnomPart_dir)
    contentF = pd.read_csv(fnomPart_dir, sep=";", header=None)
    condcol = contentF[2]
    stimDur = []

    if condname == 'InCongru-InCongru' and condname1 == 'InCongru-InCongru':
        incongIndx = [scnt for scnt, s in enumerate(condcol) if "Incong" in s]
        stimIncong = contentF.loc[contentF.index[incongIndx]]
        stimDur = stimIncong[1] - stimIncong[0]
    elif condname == 'InCongru-Congru' and condname1 == 'InCongru-InCongru':
        congIndx = [scnt1 for scnt1, s1 in enumerate(condcol) if "Cong" in s1]
        stimCong = contentF.loc[contentF.index[congIndx]]
        stimDur = stimCong[1] - stimCong[0]
    elif condname1 == 'Congru-Congru':
        stimDur = contentF[1] - contentF[0]

    stimDur_mean = np.mean(stimDur)
    stimDur_std = np.std(stimDur)
    stimDur_stdLow = stimDur_mean - stimDur_std
    stimDur_stdHigh = stimDur_mean + stimDur_std

    return stimDur_mean, stimDur_stdLow, stimDur_stdHigh


dir_base = '/Users/bolger/Documents/work/Projects/project_Brain-IHM'
dir_soundLen = os.path.join(dir_base,'SoundFile_Data')
Groups = ['Agent', 'Agent']

# Define the video-type and the feedback types.
Conds2plot = ['InCongru-Congru', 'Congru-Congru']  # Needs to be the same length as groups.

if len(Conds2plot) == 3:
    print(f'Need to collapse conditions: {Conds2plot[0]} and {Conds2plot[1]}')

# Initialize the data.
EvokedAll_cond = []
EvokedAll = []
EvokedAll_CSD = []
EvokedAllCSD_cond = []
condnames = []
Lims_upper = []
Lims_lower = []
channoms = []
all_meanSDur = []
all_stdSDur_low = []
all_stdSDur_high = []

for gcount, gcurr in enumerate(Groups):
    Group_dir = 'Data_' + gcurr
    Group_path = os.path.join(dir_base, Group_dir, 'EpochData')
    GPath_contents = os.listdir(Group_path)
    datadir_noms = [x for xind, x in enumerate(GPath_contents) if
                    x.startswith('S')]  # Extract names of the data folders.
    condcurr = Conds2plot[gcount]
    print(condcurr)
    [meanSDur, stdSDur_low, stdSDur_high] = stimLen_Calc(gcurr, condcurr,
                                                         dir_soundLen)  # Call of function to calculate the average and std stim lengths for current conditions.
    all_meanSDur.append(meanSDur)
    all_stdSDur_low.append(stdSDur_low)
    all_stdSDur_high.append(stdSDur_high)

    if condcurr == 'CongruentAll':
        condcurrf = 'CongruentAll'
    else:
        condcurrf = gcurr + '-' + condcurr
    condnames.append(condcurrf)
    EvokedLoad = {}
    EvokedCSD = {}
    evokeddata_cat = []

    for scount, scurr in enumerate(datadir_noms):
        Sujpath_curr = os.path.join(Group_path, scurr)
        sujcurr_contents = os.listdir(Sujpath_curr)
        evoked2load = [x1 for x1 in sujcurr_contents if condcurrf in x1]

        if len(evoked2load) > 0:
            Evked_load = [ev_id for ev_id in evoked2load if '-ave' in ev_id]
            print(f' Loading the evoked file: {Evked_load[0]}')
            evoked2load_path = os.path.join(Sujpath_curr, Evked_load[0])
            condsplit = condcurr.split('-')
            cond2load = condsplit[len(condsplit) - 1]
            E = mne.read_evokeds(evoked2load_path, condition=None, baseline=None, kind='average')
            EvokedLoad[scount] = E[0]
            EvokedLoad[scount].copy().set_eeg_reference(projection=True).apply_proj()
            EvokedCSD[scount] = mne.preprocessing.compute_current_source_density(EvokedLoad[scount])
            Evoked_data = EvokedLoad[scount].get_data()
            EvokedCSD_data = EvokedCSD[scount].get_data()

            if scount == 0:
                evokeddata_cat = Evoked_data
                evokedcsddata_cat = EvokedCSD_data
            elif scount > 0:
                evokeddata_cat = np.dstack((evokeddata_cat, Evoked_data))
                evokedcsddata_cat = np.dstack((evokedcsddata_cat, EvokedCSD_data))

        elif len(evoked2load) == 0:
            print(f'The condition file {condcurr} does not exist for {scurr}, {gcurr}.')

    Evokeddata_mean = np.average(evokeddata_cat, 2)
    EvokedCSDdata_mean = np.average(evokedcsddata_cat, 2)
    EvokedAll.append(EvokedLoad)
    EvokedAll_CSD.append(EvokedCSD)

    if gcount == 0:
        channoms = EvokedLoad[0].info['ch_names']
        times = EvokedLoad[0].times
        EvokedAll_cond = Evokeddata_mean
        EvokedAllCSD_cond = EvokedCSDdata_mean


    elif gcount > 0:
        EvokedAll_cond = np.dstack((EvokedAll_cond, Evokeddata_mean))
        EvokedAllCSD_cond = np.dstack((EvokedAllCSD_cond, EvokedCSDdata_mean))

    limup, limlow = sem_calc(evokeddata_cat)  # Should output standard deviation from mean for each channel.
    Lims_upper.append(limup)
    Lims_lower.append(limlow)

## Plot the ERP data for pre-defined electrodes
roi = ['F3', 'Fz', 'F4', 'FC3', 'FCz', 'FC4', 'C3', 'Cz', 'C4', 'CP3', 'CPz', 'CP4', 'P3', 'Pz', 'P4']  #
eindx = [channoms.index(elabels) for elabels in roi]
shape_ev = np.shape(EvokedAll_cond)

cols = 3
rows = int(len(roi) / cols)
fig, axes = plt.subplots(rows, cols, figsize=(30, 10))
counter = 0
colrs = ['r', 'g']
x0, y0, width, height = -0.5, 1.85, 0, 0   # bbox_to_anchor legend parameters.

for axs, ecurr in zip(axes.ravel(), eindx):

    upper1 = Lims_upper[0]
    upper2 = Lims_upper[1]

    lower1 = Lims_lower[0]
    lower2 = Lims_lower[1]

    y1 = EvokedAll_cond[ecurr, :, 0]
    y2 = EvokedAll_cond[ecurr, :, 1]
    x = times
    axs.plot(times, EvokedAll_cond[ecurr, :, 0], 'darkgray', label=condnames[0])
    axs.plot(times, upper1[ecurr, :], 'darkgray')
    axs.plot(times, lower1[ecurr, :], 'darkgray')
    axs.fill_between(times, lower1[ecurr, :], upper1[ecurr, :], facecolor='silver', alpha=0.75, linewidth=0.75)
    hdlA1 = axs.annotate('', xy=(all_meanSDur[0], -1 * (pow(10, -6))), xytext=(0, -1 * (pow(10, -6))),
                 arrowprops=dict(facecolor='red', width=2, headwidth=5), label='Mean FB duration(ms): '+condnames[0])
    styleA1 = ArrowStyle('<->', head_length=2, head_width=1.5)
    arrow1 = FancyArrowPatch((all_stdSDur_high[0], -1.5 * (pow(10, -6))), (all_stdSDur_low[0], -1.5 * (pow(10, -6))),
                             arrowstyle=styleA1, color='red', label='Mean and STD of FB length(ms): '+condnames[0])
    axs.add_patch(arrow1)

    axs.plot(times, EvokedAll_cond[ecurr, :, 1], 'steelblue', label=condnames[1])
    axs.plot(times, upper2[ecurr, :], 'steelblue')
    axs.plot(times, lower2[ecurr, :], 'steelblue')
    axs.fill_between(times, lower2[ecurr, :], upper2[ecurr, :], facecolor='steelblue', alpha=0.75, linewidth=0.5)
    hdlA2 = axs.annotate('', xy=(all_meanSDur[1], -2 * (pow(10, -6))), xytext=(0, -2 * (pow(10, -6))),
                 arrowprops=dict(facecolor='green', width=2, headwidth=5), label='Mean FB duration(ms): '+condnames[1])
    styleA2 = ArrowStyle('<->', head_length=2, head_width=1.5)
    arrow2 = FancyArrowPatch((all_stdSDur_high[1], -2.5 * (pow(10, -6))), (all_stdSDur_low[1], -2.5 * (pow(10, -6))), arrowstyle=styleA2,
                             color='green', label='Mean and STD of FB length(ms): '+condnames[1])
    axs.add_patch(arrow2)

    axs.axvline(x=0, c="black")
    axs.axhline(y=0, c="black")
    axs.set_title(channoms[ecurr], fontsize=10)
    axs.invert_yaxis()
    axs.set_frame_on(0)
    axs.set_ylim(bottom=6 * (pow(10, -6)), top=-6 * (pow(10, -6)))
    axs.set_xticks(np.arange(-0.25, 0.85, 0.2))
    axs.tick_params(axis='x', colors='white', labelsize=8)
    axs.set_ylabel(r"$\mu$V", fontsize=9)

    if counter == 0:
        axs.legend(loc='upper left', bbox_to_anchor=(x0, y0, width, height), fontsize='small', frameon=False)


    if counter >= len(roi) - 3:
        axs.tick_params(axis='x', colors='black', labelsize=8)
        axs.set_xlabel('time (seconds)', fontsize=9)

    counter += 1

curr_ax = []


def onclick(event):
    if not event.inaxes:
        return
    else:
        curr_ax[:] = [event.inaxes]


def onselect(xmin, xmax):
    xmin_round = np.rint(xmin * 1000)
    xmax_round = np.rint(xmax * 1000)
    indmin, indmax = np.searchsorted(times * 1000, (xmin_round, xmax_round))
    indmax = min(len(times) - 1, indmax)
    region_x = times[indmin:indmax]
    region_y1 = y1[indmin:indmax]
    region_y2 = y2[indmin:indmax]
    print(region_x)

    chnom = curr_ax[0].get_title()
    print(chnom)
    fig2 = plt.figure()
    Ax1 = plt.subplot2grid((2, 2), (0, 0))
    Ax2 = plt.subplot2grid((2, 2), (0, 1))
    Ax3 = plt.subplot2grid((2, 2), (1, 0), rowspan=1, colspan=2)
    fig2.tight_layout()

    # Call of function to calculate the GFP for each signal.
    GFP1, GFP2 = gfp_calc(EvokedAll_cond, times)
    data1 = np.mean(EvokedAllCSD_cond[:, indmin:indmax, 0], axis=1)
    data2 = np.mean(EvokedAllCSD_cond[:, indmin:indmax, 1], axis=1)
    im1, cn1 = mne.viz.plot_topomap(data1, EvokedCSD[0].info, vlim=(-800 * (pow(10, -6)), 800 * (pow(10, -6))),
                                    axes=Ax1)
    im2, cn2 = mne.viz.plot_topomap(data2, EvokedCSD[0].info, vlim=(-800 * (pow(10, -6)), 800 * (pow(10, -6))),
                                    axes=Ax2)
    tmin = round(times[indmin] * 1000, 1)
    tmax = round(times[indmax] * 1000, 1)
    Ax1.set_title(condnames[0] + ' (' + str(tmin) + '-' + str(tmax) + 'ms' + ')')
    Ax2.set_title(condnames[1] + ' (' + str(tmin) + '-' + str(tmax) + 'ms' + ')')
    cax1 = fig2.colorbar(im1, ax=Ax1)
    cax2 = fig2.colorbar(im2, ax=Ax2)
    cax1.set_label(r"Magnitude ($\mu$V)")
    cax2.set_label(r"Magnitude ($\mu$V)")

    gfpmax1 = np.max(GFP1)
    gfpmax2 = np.max(GFP2)
    bigmax = np.max([gfpmax1, gfpmax2])
    Ax3.plot(times, GFP1, 'b', label=condnames[0])
    Ax3.plot(times, GFP2, 'r', label=condnames[1])
    Ax3.set_ylabel(r"RMS")
    Ax3.axvline(x=0, c="black")
    Ax3.axhline(y=0, c="black")
    Ax3.set_title("Global Field Power")
    Ax3.fill_between(times[indmin:indmax], 0, bigmax, alpha=0.25, color='b')
    Ax3.set_frame_on(0)
    plt.legend()
    fig2.canvas.draw()


span_list = [SpanSelector(
    ax,
    onselect,
    "horizontal",
    useblit=False,
    props=dict(alpha=0.2, facecolor="tab:blue"),
    interactive=True,
    drag_from_anywhere=True
) for ax in axes.ravel()]

cid = fig.canvas.mpl_connect('button_press_event', onclick)