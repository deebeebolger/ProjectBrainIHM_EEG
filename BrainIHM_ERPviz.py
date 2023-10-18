import numpy as np
import mne
import matplotlib.pyplot as plt
import glob
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
from matplotlib.widgets import SpanSelector
import math

def sem_calc(evokedIn):
   """
   Function to calculate the standard error of the mean for each channel taking into account all subjects.
   :param evokedIn: the evoked data for each subject
   :return: the upper and lower limit of sem to plot
   """
   evshape = np.shape(evokedIn)
   sem_sujs = np.std(evokedIn, 2, ddof=1)/np.sqrt(evshape[2])
   sujmean = np.mean(evokedIn,2)
   upperlim = sujmean + sem_sujs
   lowlim   = sujmean - sem_sujs

   return upperlim, lowlim

def gfp_calc(evoked_data, t):

    dataev1 = evoked_data[:, :, 0]
    dataev2 = evoked_data[:, :, 1]
    dataev1_mean = np.mean(dataev1, axis=0)**2
    dataev2_mean = np.mean(dataev2, axis=0)**2
    dataev1_rm = []
    dataev2_rm = []
    for i in range(0,len(t)):
        dev1_rm = np.sum((dataev1[:,i] - dataev1_mean[i]))/64
        dev2_rm = np.sum((dataev2[:,i] - dataev2_mean[i]))/64
        dataev1_rm.append(dev1_rm)
        dataev2_rm.append(dev2_rm)

    gfp1 = np.sqrt(np.abs(dataev1_rm))
    gfp2 = np.sqrt(np.abs(dataev2_rm))

    return gfp1, gfp2


dir_base = '/Users/bolger/Documents/work/Projects/project_Brain-IHM'
Groups = ['Agent', 'Agent']
# Define the video-type and the feedback types.
Conds2plot = ['InCongru-Congru', 'InCongru-InCongru']  # Needs to be the same length as groups.

# Initialize the data.
EvokedAll_cond = []
EvokedAll = []
EvokedAll_CSD = []
EvokedAllCSD_cond = []
condnames = []
Lims_upper = []
Lims_lower = []
channoms = []

for gcount, gcurr in enumerate(Groups):
    Group_dir = 'Data_' + gcurr
    Group_path = os.path.join(dir_base, Group_dir, 'EpochData')
    GPath_contents = os.listdir(Group_path)
    datadir_noms = [x for xind, x in enumerate(GPath_contents) if x.startswith('S')] # Extract names of the data folders.
    condcurr = Conds2plot[gcount]
    condcurrf = gcurr + '-' + condcurr
    condnames.append(condcurrf)
    EvokedLoad = {}
    EvokedCSD = {}
    evokeddata_cat = []

    for scount, scurr in enumerate(datadir_noms):
        Sujpath_curr = os.path.join(Group_path, scurr)
        sujcurr_contents = os.listdir(Sujpath_curr)
        evoked2load = [x1 for x1 in sujcurr_contents if condcurrf in x1]

        if len(evoked2load)>0:
            print(f' Loading the evoked file: {evoked2load[0]}')
            evoked2load_path = os.path.join(Sujpath_curr, evoked2load[0])
            condsplit = condcurr.split('-')
            cond2load = condsplit[1]
            EvokedLoad[scount] = mne.read_evokeds(evoked2load_path,condition=cond2load, baseline=None, kind='average')
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

        elif len(evoked2load)==0:
            print(f'The condition file {condcurr} does not exist for {scurr}, {gcurr}.')



    Evokeddata_mean = np.average(evokeddata_cat,2)
    EvokedCSDdata_mean = np.average(evokedcsddata_cat, 2)
    EvokedAll.append(EvokedLoad)
    EvokedAll_CSD.append(EvokedCSD)

    if gcount == 0:
        channoms  = EvokedLoad[0].info['ch_names']
        times = EvokedLoad[0].times
        EvokedAll_cond = Evokeddata_mean
        EvokedAllCSD_cond = EvokedCSDdata_mean


    elif gcount>0:
        EvokedAll_cond = np.dstack((EvokedAll_cond, Evokeddata_mean))
        EvokedAllCSD_cond = np.dstack((EvokedAllCSD_cond, EvokedCSDdata_mean))

    limup, limlow = sem_calc(evokeddata_cat)  # Should output standard deviation from mean for each channel.
    Lims_upper.append(limup)
    Lims_lower.append(limlow)


## Plot the ERP data for pre-defined electrodes
roi = ['F3','Fz','F4','FC3','FCz','FC4','C3', 'Cz', 'C4', 'CP3','CPz', 'CP4', 'P3', 'Pz', 'P4']  #
eindx = [channoms.index(elabels) for elabels in roi]
shape_ev = np.shape(EvokedAll_cond)


cols = 3
rows = int(len(roi)/cols)
fig, axes = plt.subplots(rows, cols, figsize=(30,10))
counter = 0

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
    axs.fill_between(times, lower1[ecurr, :], upper1[ecurr,:], facecolor='silver', alpha=0.75, linewidth=0.75)

    axs.plot(times, EvokedAll_cond[ecurr, :, 1], 'steelblue', label=condnames[1])
    axs.plot(times, upper2[ecurr, :], 'steelblue')
    axs.plot(times, lower2[ecurr, :], 'steelblue')
    axs.fill_between(times, lower2[ecurr, :], upper2[ecurr, :], facecolor='steelblue', alpha=0.75, linewidth=0.5)
    axs.axvline(x=0, c="black")
    axs.axhline(y=0, c="black")
    axs.set_title(channoms[ecurr])
    axs.invert_yaxis()
    axs.set_frame_on(0)
    axs.set_ylim(bottom=6*(pow(10, -6)), top=-6*(pow(10, -6)))
    axs.set_xticks(np.arange(-0.25, 0.85, 0.2))
    axs.tick_params(axis='x', colors='white')
    axs.set_ylabel(r"$\mu$V")
    axs.legend(loc='upper right', fontsize='x-small', frameon=False)

    if counter >= len(roi)-3:
        axs.tick_params(axis='x', colors='black')
        axs.set_xlabel('time (seconds)')

    counter+=1



curr_ax = []
def onclick(event):
    if not event.inaxes:
        return
    else:
        curr_ax[:] = [event.inaxes]

def onselect(xmin, xmax):

    indmin, indmax = np.searchsorted(times, (xmin, xmax))
    indmax = min(len(times) - 1, indmax)
    region_x = times[indmin:indmax]
    region_y1 = y1[indmin:indmax]
    region_y2 = y2[indmin:indmax]
    print(region_x)

    chnom = curr_ax[0].get_title()
    print(chnom)
    fig2 = plt.figure()
    Ax1 = plt.subplot2grid((2,2), (0,0))
    Ax2 = plt.subplot2grid((2,2), (0,1))
    Ax3 = plt.subplot2grid((2,2), (1,0), rowspan=1, colspan=2)
    fig2.tight_layout()


    # Call of function to calculate the GFP for each signal.
    GFP1, GFP2 = gfp_calc(EvokedAllCSD_cond, times)
    data1 = np.mean(EvokedAllCSD_cond[:, indmin:indmax, 0], axis=1)
    data2 = np.mean(EvokedAllCSD_cond[:, indmin:indmax, 1], axis=1)
    im1, cn1 = mne.viz.plot_topomap(data1, EvokedCSD[0].info, vlim=(-800*(pow(10, -6)), 800*(pow(10, -6))), axes=Ax1)
    im2, cn2 = mne.viz.plot_topomap(data2, EvokedCSD[0].info, vlim=(-800* (pow(10, -6)), 800* (pow(10, -6))), axes=Ax2)
    tmin = round(times[indmin]*1000, 1)
    tmax = round(times[indmax] * 1000, 1)
    Ax1.set_title(condnames[0]+' ('+ str(tmin)+'-'+str(tmax)+'ms'+')')
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
    Ax3.fill_between(times[indmin:indmax], 0, bigmax, alpha=0.25, color= 'b')
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




