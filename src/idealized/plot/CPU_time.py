import numpy as np
import matplotlib.pyplot as plt
# compute CPU time
def getdof(filename, nt):
    f = open(filename, 'r')
    print(filename)
    refine_time = 0.
    scheme_time = 0.
    dof = 0
    error = np.zeros([nt + 1, 3])
    mass  = np.zeros([nt + 1])
    line = f.readline()
    while line:
        string = line.split()
        if len(string) >= 2:
            if string[0] == 'nColumns' and string[1] != '0':
                refine_time += -float(string[-2])
                scheme_time += -float(string[-1])
                dof         += int(string[-3])
            # elif string[1] == 'time':
            #     error[int(string[2]), 0] = float(string[5])
            #     error[int(string[2]), 1] = float(string[8]) 
            #     error[int(string[2]), 2] = float(string[11])
            #     mass[int(string[2])]     = float(string[-1])
        line = f.readline()
    f.close()
    return refine_time/nt, scheme_time/nt, dof/nt
# 
nx        = np.array([72, 144, 288, 576])
dof       = np.array([2592, 10368, 41472, 165888])
snt       = [12*7, 12*7*2, 12*7*4, 12*7*8]

rf0_refine_time  = np.zeros(4)
rf0_scheme_time  = np.zeros(4)
rf0_dof = np.zeros(4)
for i, ix, nt in zip(range(4), nx, snt):
    rf0_refine_time[i], rf0_scheme_time[i], rf0_dof[i] = getdof('../solid0./rf0/plots_solid_'+str(ix)+'_rf0_'+str(nt)+'_0./nx'+str(ix)+'_rf0.log', nt)

rf1_refine_time  = np.zeros(4)
rf1_scheme_time  = np.zeros(4)
rf1_dof = np.zeros(4)
for i, ix, nt in zip(range(4), nx, snt):
    rf1_refine_time[i], rf1_scheme_time[i], rf1_dof[i] = getdof('../solid0./rf1/plots_solid_'+str(ix//2)+'_rf1_'+str(nt)+'_0./nx'+str(ix//2)+'_rf1.log', nt)

rf2_refine_time  = np.zeros(4)
rf2_scheme_time  = np.zeros(4)
rf2_dof = np.zeros(4)
for i, ix, nt in zip(range(4), nx, snt):
    rf2_refine_time[i], rf2_scheme_time[i], rf2_dof[i]  = getdof('../solid0./rf2/plots_solid_'+str(ix//4)+'_rf2_'+str(nt)+'_0./nx'+str(ix//4)+'_rf2.log', nt)


print(rf0_scheme_time/dof, rf1_scheme_time/dof, rf2_scheme_time/dof)
# dx = 360./nx

fig = plt.figure(1, figsize=(9, 6.75))
figlegend = plt.figure(figsize = (28, 6))
plt.clf()
ax = fig.add_subplot(111)
ax.loglog(dof[:], rf0_scheme_time, marker = '^', linestyle = 'solid',  markersize=15, color = 'k', label = 'no refinement')
ax.loglog(rf1_dof, rf1_scheme_time, marker = 'o', linestyle = 'dotted',  markersize=15, color = 'k', label = 'one level refinement')

ax.loglog(rf2_dof, rf2_scheme_time    , marker = '^', linestyle =  'dashed', fillstyle='none',  markersize=15,  color = 'k', label = 'two level refinement')

# ax.loglog([dof[0], dof[-1]], [rf0_scheme_time[0], rf0_scheme_time[0]*(dof[-1]/dof[0])]    , marker = '', linestyle =  'dashdot', fillstyle='none',  markersize=15,  color = 'k', label = 'linear')

lgnd = figlegend.legend(*ax.get_legend_handles_labels(), 'center', prop = {'size': 28}, ncol=2, handlelength=2.)
lgnd.legendHandles[0]._legmarker.set_markersize(20)
lgnd.legendHandles[1]._legmarker.set_markersize(20)
lgnd.legendHandles[2]._legmarker.set_markersize(20)
# lgnd.legendHandles[3]._legmarker.set_markersize(20)
ax.set_xticks(dof[:-1])
ax.set_xticklabels(dof[:-1])
# ax.set_yticks(np.linspace(0., 2., 11))
# ax.set_yticklabels(np.linspace(0., 2., 11))
ax.set_xlim(np.min(rf1_dof) - 300, dof[-1]+20000)
ax.set_ylim([0., 2.])
ax.set_xlabel("cell number",fontsize=20)
ax.set_ylabel("CPU time of numerical scheme",fontsize=20)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
# plt.show()
fig.savefig('CPU_scheme_time.pdf')
figlegend.savefig('CPU_legend.pdf')
plt.close()

fig = plt.figure(1, figsize=(9, 6.75))
# figlegend = plt.figure(figsize = (28, 6))
plt.clf()
ax = fig.add_subplot(111)
# ax.loglog(dof, rf0_refine_time, marker = '^', linestyle = 'solid',  markersize=15, color = 'k', label = 'rf0')
ax.plot(rf1_dof, rf1_refine_time/(rf1_refine_time + rf1_scheme_time)*100, marker = 'o', linestyle = 'dotted',  markersize=15, color = 'k', label = 'rf1')

ax.plot(rf2_dof, rf2_refine_time /(rf2_refine_time + rf2_scheme_time) *100  , marker = '^', linestyle = 'dashed', fillstyle='none',  markersize=15,  color = 'k', label = 'rf2')

# lgnd = figlegend.legend(*ax.get_legend_handles_labels(), 'center', prop = {'size': 28}, ncol=2, handlelength=2.)
# lgnd.legendHandles[0]._legmarker.set_markersize(20)
# lgnd.legendHandles[1]._legmarker.set_markersize(20)
# lgnd.legendHandles[2]._legmarker.set_markersize(20)
# lgnd.legendHandles[3]._legmarker.set_markersize(20)
ax.set_xticks(dof[:-1])
ax.set_xticklabels(dof[:-1])
# ax.set_yticks(np.linspace(0., 2., 11))
# ax.set_yticklabels(np.linspace(0., 2., 11))
ax.set_xlim(rf1_dof[0] - 1000, rf1_dof[-1] + 1000)
ax.set_ylim([0., 100.])
ax.set_xlabel("cell number",fontsize=20)
ax.set_ylabel("Percent of Refinement Time",fontsize=20)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
fig.savefig('CPU_percent_time.pdf')
# figlegend.savefig('CPU_legend.pdf')
plt.close()
