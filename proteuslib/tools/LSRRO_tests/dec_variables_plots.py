import os
import numpy as np
import matplotlib.pyplot as plt
from sim_LSRRO_Nstage import build_model, simulate, optimization


def autolabel(rects,ax,fmt='{}'):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate(fmt.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', rotation=0)

def collect_data(m,N):
    ### Collect output data ###
    pp = []
    ep = []
    ma = []
    sr = []
    for i in range(N):
        pump = getattr(m.fs,"P"+repr(i+1))
        pp.append(pump.outlet.pressure[0].value / 1e5)
    for i in range(N-1):
        pump = getattr(m.fs,"EqP"+repr(i+2))
        ep.append(pump.outlet.pressure[0].value / 1e5)
    ep.append(float('nan'))
    for i in range(N):
        stage = getattr(m.fs,"Stage"+repr(i+1))
        ma.append(stage.area.value)
        sr.append(stage.permeate.flow_mass_comp[0, 'NaCl'].value/stage.inlet.flow_mass_comp[0, 'NaCl'].value * 100) 
    return [pp,ep,ma,sr]

### Set simulation parameters ###
# n, salt_frac, recovery
# cases = [
#     [2, 0.035, 0.75],
#     [3, 0.035, 0.75],
#     [4, 0.035, 0.75]
# ]
cases = [
    [2, 0.035, 0.7],
    [3, 0.070, 0.5],
    [4, 0.105, 0.4],
]

### Initialize bar charts ###
pp_fig, pp_ax = plt.subplots() 
ep_fig, ep_ax = plt.subplots() 
ma_fig, ma_ax = plt.subplots() 
sr_fig, sr_ax = plt.subplots() 

### Initialize bar chart layout data###
bar_start = 0
bar_ticks = []
bar_labels = []
w = .85

### Iterate over n stage values
for i, [n, salt_frac, recovery] in enumerate(cases):

    ### Get opt data (fake) ###
    # pp = [n/2*(60-j) for j in range(n)]
    # ep = [n/2*(60-j) for j in range(n)]
    # ma = [n*60/2**j for j in range(n)]
    # sr = [5-j for j in range(n)]

    ### Get opt data (real) ###
    m = build_model(N=n)
    m = simulate(m, N=n)
    m = optimization(m, N=n, params=["flow_mass_comp","recovery"], values=[salt_frac,recovery])
    pp, ep, ma, sr = collect_data(m, n)

    ### Add data sets and labels to the bar charts
    rects = pp_ax.bar((np.arange(n)+bar_start), pp, label='Case '+repr(i+1),edgecolor='k', linewidth=2, width=w)
    autolabel(rects,pp_ax,fmt="{:.1f}")
    rects = ep_ax.bar((np.arange(n)+bar_start), ep, label='Case '+repr(i+1),edgecolor='k', linewidth=2, width=w)
    autolabel(rects,ep_ax,fmt="{:.1f}")
    rects = ma_ax.bar((np.arange(n)+bar_start), ma, label='Case '+repr(i+1),edgecolor='k', linewidth=2, width=w)
    autolabel(rects,ma_ax,fmt="{:.1f}")
    rects = sr_ax.bar((np.arange(n)+bar_start), sr, label='Case '+repr(i+1),edgecolor='k', linewidth=2, width=w)
    autolabel(rects,sr_ax,fmt="{:.1f}")

    ### Update values for bar chart formatting
    bar_ticks += [(bar_start+j) for j in range(n)]
    bar_labels += [j+1 for j in range(n)]
    bar_start += (n+.5)

### Create output folder ###
if not os.path.exists("plots/"): os.makedirs("plots/")

### Generate pump pressure chart ###
# pp_ax.set_title('Pump pressure for N-stage LSRRO')
pp_ax.set_ylabel('Pressure (bar)')
pp_ax.set_xlabel('Pump')
pp_ax.set_xticks(bar_ticks)
pp_ax.set_xticklabels(bar_labels)
pp_ax.margins(y=.15)
pp_ax.set_box_aspect(1)
pp_fig.tight_layout()
yt = pp_ax.get_yticks()
# pp_ax.set_ylim(yt[0],yt[-1])
pp_ax.set_ylim(yt[0],2*yt[-1]-yt[-2])
pp_ax.legend(loc=2,frameon=False)
pp_fig.savefig("plots/pressure.png",transparent=True,bbox_inches='tight')

### Generate eq pump pressure chart ###
# ep_ax.set_title('Eq Pump pressure for N-stage LSRRO')
ep_ax.set_ylabel('Pressure (bar)')
ep_ax.set_xlabel('Eq Pump')
ep_ax.set_xticks(bar_ticks)
ep_ax.set_xticklabels(np.array(bar_labels)+1)
ep_ax.margins(y=.15)
ep_ax.set_box_aspect(1)
ep_fig.tight_layout()
yt = ep_ax.get_yticks()
ep_ax.set_ylim(yt[0],yt[-1])
ep_ax.legend(loc=2,frameon=False)
ep_fig.savefig("plots/eq_pressure.png",transparent=True,bbox_inches='tight')

### Generate membrane area chart ###
# ma_ax.set_title('Membrane area for N-stage LSRRO')
ma_ax.set_ylabel('Membrane Area (m$^2$)')
ma_ax.set_xlabel('Stage')
ma_ax.set_xticks(bar_ticks)
ma_ax.set_xticklabels(bar_labels)
ma_ax.margins(y=.15)
ma_ax.set_box_aspect(1)
ma_fig.tight_layout()
yt = ma_ax.get_yticks()
ma_ax.set_ylim(yt[0],yt[-1])
ma_ax.legend(loc=2,frameon=False)
ma_fig.savefig("plots/area.png",transparent=True,bbox_inches='tight')

## Generate salt rejection chart ###
# sr_ax.set_title('Salt Passage for N-stage LSRRO')
sr_ax.set_ylabel('Salt Passage (%)')
sr_ax.set_xlabel('Stage')
sr_ax.set_xticks(bar_ticks)
sr_ax.set_xticklabels(bar_labels)
sr_ax.margins(y=.15)
sr_ax.set_box_aspect(1)
sr_fig.tight_layout()
yt = sr_ax.get_yticks()
sr_ax.set_ylim(yt[0],yt[-1])
sr_ax.legend(loc=2,frameon=False)
sr_fig.savefig("plots/salt_passage.png",transparent=True,bbox_inches='tight')

plt.show()