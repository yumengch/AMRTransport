import numpy as np
import matplotlib.pyplot as plt
# compute cell numbers (degrees of freedom)
def getdof(filename, nt):
    f = open(filename, 'r')
    print(filename)
    val = 0
    error = np.zeros([nt + 1, 3])
    line = f.readline()
    while line:
        line = f.readline()
        string = line.split()
        # print(string)
        if len(string) >= 2:
            if string[0] == 'nColumns':
                val += int(string[2])
            elif string[1] == 'time':
                error[int(string[2]), 0] = float(string[5])
                error[int(string[2]), 1] = float(string[8]) 
                error[int(string[2]), 2] = float(string[11]) 
    f.close()
    return val/(nt + 1), error[-1, :]

nx        = np.array([72, 144, 288, 576])
dof       = np.array([2592, 10368, 41472, 165888])
snt       = [55, 110, 220, 440]

inter_rf1_dof  = np.zeros(4)
inter_rf1_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, snt):
    inter_rf1_dof[i], inter_rf1_error[i, :]  = getdof('../solid0.15intermediate/rf1/plots_solid_'+str(ix//2)+'_rf1_'+str(nt)+'_0.15/nx'+str(ix//2)+'_rf1.log', nt)

inter_rf2_dof  = np.zeros(4)
inter_rf2_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, snt):
    inter_rf2_dof[i], inter_rf2_error[i, :]  = getdof('../solid0.15intermediate/rf2/plots_solid_'+str(ix//4)+'_rf2_'+str(nt)+'_0.15/nx'+str(ix//4)+'_rf2.log', nt)

s_rf1_dof  = np.zeros(4)
s_rf1_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, snt):
    s_rf1_dof[i], s_rf1_error[i, :]  = getdof('../solid0.15/rf1/plots_solid_'+str(ix//2)+'_rf1_'+str(nt)+'_0.15/nx'+str(ix//2)+'_rf1.log', nt)

s_rf2_dof  = np.zeros(4)
s_rf2_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, snt):
    s_rf2_dof[i], s_rf2_error[i, :]  = getdof('../solid0.15/rf2/plots_solid_'+str(ix//4)+'_rf2_'+str(nt)+'_0.15/nx'+str(ix//4)+'_rf2.log', nt)

print(inter_rf2_dof, s_rf2_dof)

snt       = [13, 26, 52, 104]

inter_rf1_dof_0  = np.zeros(4)
inter_rf1_error_0  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, snt):
    inter_rf1_dof_0[i], inter_rf1_error_0[i, :]  = getdof('../solid0.intermediate/rf1/plots_solid_'+str(ix//2)+'_rf1_'+str(nt)+'_0./nx'+str(ix//2)+'_rf1.log', nt)

inter_rf2_dof_0  = np.zeros(4)
inter_rf2_error_0  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, snt):
    inter_rf2_dof_0[i], inter_rf2_error_0[i, :]  = getdof('../solid0.intermediate/rf2/plots_solid_'+str(ix//4)+'_rf2_'+str(nt)+'_0./nx'+str(ix//4)+'_rf2.log', nt)

s_rf1_dof_0  = np.zeros(4)
s_rf1_error_0  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, snt):
    s_rf1_dof_0[i], s_rf1_error_0[i, :]  = getdof('../solid0./rf1/plots_solid_'+str(ix//2)+'_rf1_'+str(nt)+'_0./nx'+str(ix//2)+'_rf1.log', nt)

s_rf2_dof_0  = np.zeros(4)
s_rf2_error_0 = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, snt):
    s_rf2_dof_0[i], s_rf2_error_0[i, :]  = getdof('../solid0./rf2/plots_solid_'+str(ix//4)+'_rf2_'+str(nt)+'_0./nx'+str(ix//4)+'_rf2.log', nt)

dx = 360./nx

fig = plt.figure(1, figsize=(9, 6.75))
figlegend = plt.figure(figsize = (28, 6))
plt.clf()
ax = fig.add_subplot(111)
ax.plot(dx, (inter_rf1_dof - s_rf1_dof)/s_rf1_dof*100, marker = '^', linestyle = '-',  markersize=15, color = 'k', label = r'$\alpha = \frac{3}{20} \pi $')
ax.plot(dx, (inter_rf1_dof_0 - s_rf1_dof_0)/s_rf1_dof_0*100, marker = 'o', linestyle = '-',  markersize=15, color = 'k', label = r'$\alpha = 0 $')
# ax.plot(dx, (inter_rf2_dof - s_rf2_dof)/s_rf2_dof, marker = 'o', linestyle = '-',  markersize=15, color = 'k', label = 'two level refinement')

# ax.loglog(dx, s_rf1_dof    , marker = '^', linestyle = '-', fillstyle='none',  markersize=15,  color = 'k', label = 'one level refinement without intermediate step')
# ax.loglog(dx, s_rf2_dof    , marker = 'o', linestyle = '-', fillstyle='none',  markersize=15,  color = 'k', label = 'two level refinement without intermediate step')

lgnd = figlegend.legend(*ax.get_legend_handles_labels(), 'center', prop = {'size': 28}, ncol=2, handlelength=2.)
lgnd.legendHandles[0]._legmarker.set_markersize(20)
lgnd.legendHandles[1]._legmarker.set_markersize(20)
# lgnd.legendHandles[2]._legmarker.set_markersize(20)
# lgnd.legendHandles[3]._legmarker.set_markersize(20)
ax.set_xticks(dx)
ax.set_xticklabels(dx)
ax.set_xlim([0.5, 6])
ax.set_ylim([0, 20])
ax.set_xlabel("maximum resolution",fontsize=20)
ax.set_ylabel("cell number difference (%)",fontsize=20)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
fig.savefig('dof_diff.pdf')
figlegend.savefig('legend.pdf')
plt.close()
