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

def sen_bar_chart(sinput,data,fmt="{:.1f}"):
    n_cases = len(data)
    n_sens = len(sinput)

    bar_width  = 1/(n_cases)
    padding = bar_width
    bar_ticks  = [(1+padding)*j for j in range(n_sens)]
    bar_labels = sinput

    fig, ax = plt.subplots() 
    for i, elem in enumerate(data):
        offset = i*bar_width-(0.5-bar_width/2)
        rects = ax.bar((np.array(bar_ticks)+offset), elem, label='Case '+repr(i+1),edgecolor='k', linewidth=2, width=bar_width)
        # autolabel(rects,ax,fmt=fmt)

    ax.set_xticks(bar_ticks)
    ax.set_xticklabels(bar_labels)
    ax.margins(y=.15)
    ax.set_box_aspect(1)
    yt = ax.get_yticks()
    ax.set_ylim(yt[0],6)
    # ax.set_ylim(yt[0],2*yt[-1]-yt[-2])
    ax.legend(frameon=False) #loc=2

    return [fig, ax]

### Set simulation parameters ###
# n, salt_frac, recovery
cases = [
    [2, 0.035, 0.7],
    [3, 0.070, 0.5],
    [4, 0.105, 0.4],
]

### Set the costs for sen analysis   0.07 $/kWh and 30 $/m2        ,
base_ele_cost   = 0.07    # $/kWh
base_mem_cost   = 30.0    # $/m2
base_water_perm = 4.2e-12 # m/s-Pa
ele_costs  = [0.03,0.05,0.07,0.09,0.11]       # 0.03-0.15 $/kWh
mem_costs  = [10.0,20.0,30.0,40.0,50.0]       # 20-50 $/m2
water_perms =[2e-12,4e-12,6e-12,8e-12,10e-12]       # 1 to 10 [1e-12 m/s-Pa]

### initialize blank data lists ###
sim_ele_list = []
sim_mem_list = []
sim_prm_list = []

### Iterate over n stage values
for i, [n, salt_frac, recovery] in enumerate(cases):

    ### Get opt data (fake)
    # sim_ele_list.append([])
    # for ele_cost in ele_costs:
    #     sim_ele_list[i].append(np.random.randint(1,100))

    # sim_mem_list.append([])
    # for mem_cost in mem_costs:
    #     sim_mem_list[i].append(np.random.randint(1,100))

    # sim_prm_list.append([])
    # for water_perm in water_perms:
    #     sim_prm_list[i].append(np.random.randint(1,100))

    ### Get opt data (real) ###
    m = build_model(N=n)
    m = simulate(m, N=n)

    sim_ele_list.append([])
    for ele_cost in ele_costs:
        m = optimization(m, N=n, 
            params=["flow_mass_comp","recovery","ele_cost","mem_cost","water_perm"], 
            values=[salt_frac,recovery,ele_cost,base_mem_cost,base_water_perm])
        sim_ele_list[i].append(float(m.fs.costing.LCOW.value))

    sim_mem_list.append([])
    for mem_cost in mem_costs:
        m = optimization(m, N=n, 
            params=["flow_mass_comp","recovery","ele_cost","mem_cost","water_perm"], 
            values=[salt_frac,recovery,base_ele_cost,mem_cost,base_water_perm])
        sim_mem_list[i].append(float(m.fs.costing.LCOW.value))

    sim_prm_list.append([])
    for water_perm in water_perms:
        m = optimization(m, N=n, 
            params=["flow_mass_comp","recovery","ele_cost","mem_cost","water_perm"], 
            values=[salt_frac,recovery,base_ele_cost,base_mem_cost,water_perm])
        sim_prm_list[i].append(float(m.fs.costing.LCOW.value))
        

### Create output folder ###
if not os.path.exists("plots/"): os.makedirs("plots/")

print()
print("    Ele List: "+repr(ele_costs))
print("    Mem List: "+repr(mem_costs))
print("    Per List: "+repr(water_perms))
print()
for i in range(len(cases)):
    print("Case "+repr(i+1))
    print("    Ele Costs: ")
    for v in sim_ele_list[i]:
        print("        "+repr(v))
    print("    Mem Costs: ")
    for v in sim_mem_list[i]:
        print("        "+repr(v))
    print("    Per Costs: ")
    for v in sim_prm_list[i]:
        print("        "+repr(v))
    print()

fig, ax = sen_bar_chart(ele_costs,sim_ele_list)
ax.set_ylabel(r'LCOW (\$/m$^3$)')
ax.set_xlabel(r'Electricity Cost (\$/kWh)')
fig.tight_layout()
fig.savefig("plots/ele_cost.png",transparent=True,bbox_inches='tight')

fig, ax = sen_bar_chart(mem_costs,sim_mem_list)
ax.set_ylabel(r'LCOW (\$/m$^3$)')
ax.set_xlabel(r'Membrane Cost (\$/m$^2$)')
fig.tight_layout()
fig.savefig("plots/mem_cost.png",transparent=True,bbox_inches='tight')

fig, ax = sen_bar_chart(np.array(water_perms)*1e12,sim_prm_list)
ax.set_ylabel(r'LCOW (\$/m$^3$)')
ax.set_xlabel(r'Water Permeability (1e-12 m/s-Pa)')
fig.tight_layout()
fig.savefig("plots/water_perm.png",transparent=True,bbox_inches='tight')

plt.show()


