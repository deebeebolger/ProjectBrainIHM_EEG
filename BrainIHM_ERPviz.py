import numpy as np
import mne
import matplotlib.pyplot as plt
import glob
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
from matplotlib.widgets import SpanSelector




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


dir_base = '/Users/bolger/Documents/work/Projects/project_Brain-IHM'
Groups = ['Human', 'Agent']
# Define the video-type and the feedback types.
Conds2plot = ['Congru-Congru', 'Congru-Congru']  # Needs to be the same length as groups.

# Initialize the data.
EvokedAll_cond = []
EvokedAll = []
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
            Evoked_data = EvokedLoad[scount].get_data()

            if scount == 0:
                evokeddata_cat = Evoked_data
            elif scount > 0:
                evokeddata_cat = np.dstack((evokeddata_cat, Evoked_data))

        elif len(evoked2load)==0:
            print(f'The condition file {condcurr} does not exist for {scurr}, {gcurr}.')



    Evokeddata_mean = np.average(evokeddata_cat,2)
    EvokedAll.append(EvokedLoad)

    if gcount == 0:
        channoms  = EvokedLoad[0].info['ch_names']
        times = EvokedLoad[0].times
        EvokedAll_cond = Evokeddata_mean

    elif gcount>0:
        EvokedAll_cond = np.dstack((EvokedAll_cond, Evokeddata_mean))

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
    axs.plot(times, EvokedAll_cond[ecurr, :, 0], 'b-', label=condnames[0])
    axs.plot(times, upper1[ecurr, :], 'b-')
    axs.plot(times, lower1[ecurr, :], 'b-')
    axs.fill_between(times, lower1[ecurr, :], upper1[ecurr,:])

    axs.plot(times, EvokedAll_cond[ecurr, :, 1], 'r-', label=condnames[1])
    axs.plot(times, upper2[ecurr, :], 'r-')
    axs.plot(times, lower2[ecurr, :], 'r-')
    axs.fill_between(times, lower2[ecurr, :], upper2[ecurr, :])
    axs.axvline(x=0, c="black")
    axs.axhline(y=0, c="black")
    axs.set_title(channoms[ecurr])
    axs.invert_yaxis()
    axs.set_frame_on(0)
    axs.set_ylim(bottom=6*(pow(10, -6)), top=-3*(pow(10, -6)))
    axs.set_xticks(np.arange(-0.25, 0.85, 0.2))
    axs.tick_params(axis='x', colors='white')
    axs.set_ylabel(r"$\mu$V")

    if counter >= len(roi)-3:
        axs.tick_params(axis='x', colors='black')
        axs.set_xlabel('time (seconds)')

    counter+=1

plt.legend()

chnom = []
curr_ax = []
axes = []
def onclick(event):
    if not event.inaxes:
        return
    else:
        curr_ax[:] = [event.inaxes]

def onselect(xmin, xmax):

    indmin, indmax = np.searchsorted(x, (xmin, xmax))
    indmax = min(len(times) - 1, indmax)
    region_x = times[indmin:indmax]
    region_y1 = y1[indmin:indmax]
    region_y2 = y2[indmin:indmax]
    print(region_x)
    for ax, span in zip(axes, span_list):
        if ax != curr_ax[0]:
            span.set_visible(False)
    fig.canvas.draw_idle()

    chnom = curr_ax[0].get_title()
    print(chnom)
    fig2, axes = plt.subplots(1, 2)

    data1 = np.mean(EvokedAll_cond[:, indmin:indmax, 0], axis=1)
    data2 = np.mean(EvokedAll_cond[:, indmin:indmax, 1], axis=1)
    im1, cn1 = mne.viz.plot_topomap(data1, EvokedLoad[0].info, vlim=(-1*(pow(10, -6)), 7*(pow(10, -6))), axes=axes1)
    im2, cn2 = mne.viz.plot_topomap(data2, EvokedLoad[0].info, vlim=(-1 * (pow(10, -6)), 7 * (pow(10, -6))), axes=axes2)
    tmin = round(times[indmin]*1000, 1)
    tmax = round(times[indmax] * 1000, 1)
    axes[0].set_title(condnames[0]+' ('+ str(tmin)+'-'+str(tmax)+'ms'+')')
    axes[1].set_title(condnames[1] + ' (' + str(tmin) + '-' + str(tmax) + 'ms' + ')')
    cax1 = fig2.colorbar(im1, ax=axes[0])
    cax2 = fig2.colorbar(im2, ax=axes[1])
    cax1.set_label(r"Magnitude ($\mu$V)")
    cax2.set_label(r"Magnitude ($\mu$V)")
    fig2.canvas.draw()


span_list = [SpanSelector(
    ax,
    onselect,
    "horizontal",
    useblit=True,
    props=dict(alpha=0.2, facecolor="tab:blue"),
    interactive=True,
    drag_from_anywhere=True
) for ax in axes]

fig.canvas.mpl_connect('button_press_event', onclick)

plt.show()



