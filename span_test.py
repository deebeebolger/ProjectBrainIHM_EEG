import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
import numpy as np

#sample data
fig, axis = plt.subplots(3)
for i, ax in enumerate(axis):
    t=np.linspace(-i, i+1 , 100)
    ax.plot(t, np.sin(2 * np.pi * t))

#list to store the axis last used with a mouseclick
curr_ax = []

#detect the currently modified axis
def on_click(event):
    if event.inaxes:
        curr_ax[:] = [event.inaxes]

#modify the current axis objects
def onselect(xmin, xmax):
    #ignore if accidentally clicked into an axis object
    if xmin==xmax:
        return
    #set all span selectors invisible accept the current
    for ax, span in zip(axis, list_of_spans):
        if ax != curr_ax[0]:
            span.set_visible(False)
    #do something with xmin, xmax
    print(xmin, xmax)
    fig.canvas.draw_idle()


#collect span selectors in a list in the same order as their axes objects
list_of_spans = [SpanSelector(
        ax,
        onselect,
        "horizontal",
        useblit=True,
        props=dict(alpha=0.5, facecolor="tab:blue"),
        interactive=True,
        drag_from_anywhere=True
        )
        for ax in axis]

plt.connect('button_press_event', on_click)

plt.show()

