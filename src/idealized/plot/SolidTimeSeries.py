import numpy as np
import matplotlib.pyplot as plt
# compute cell numbers (degrees of freedom)
def getdof(filename, nt):
    f = open(filename, 'r')
    print(filename)
    val = []
    error = np.zeros([nt + 1, 3])
    mass  = np.zeros([nt + 1])
    line = f.readline()
    while line:
        string = line.split()
        # print(string)
        if len(string) >= 2:
            if string[0] == 'nColumns':
                val.append(int(string[-1]))
            elif string[1] == 'time':
                error[int(string[2]), 0] = float(string[5])
                error[int(string[2]), 1] = float(string[8]) 
                error[int(string[2]), 2] = float(string[11])
                mass[int(string[2])]     = float(string[-1])
        line = f.readline()
    f.close()
    return np.array(val)

nx             = np.array([72, 144, 288, 576])
dof            = np.array([2592, 10368, 41472, 165888])

s05_c1_nt       = [1800, 1800*4, 1800*4*4, 1800*4*4*4]
s05_c5_nt       = [240, 240*4, 240*4*4, 240*4*4*4]
s05_c1_dof  = np.zeros(1800*4*4 + 1)
ix = 144
nt = 1800*4*4
s05_c1_dof  = getdof('../solid0.5/rf1/plots_solid_'+str( ix )+'_rf1_'+str(nt)+'_0.5/nx'+str( ix )+'_rf1.log', nt)
t1 = np.linspace(0, 12, nt + 1)


s05_c5_dof  = np.zeros(240*4*4 + 1)
ix = 144
nt = 240*4*4
s05_c5_dof  = getdof('../solid0.5/rf1/plots_solid_'+str( ix )+'_rf1_'+str(nt)+'_0.5/nx'+str( ix )+'_rf1.log', nt)
t5 = np.linspace(0, 12, nt + 1)

# print(len(t))
xt = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
fig = plt.figure(num=1, figsize=(16, 12))
figlegend = plt.figure(figsize = (20, 20))
plt.clf()
ax = fig.add_subplot(111)
ax.plot(t1, s05_c1_dof, 'k-', linewidth=4.0, label = 'Maximum Courant number around 1')
ax.plot(t5, s05_c5_dof, 'k:', linewidth=4.0, label = 'Maximum Courant number around 5')

lgnd = figlegend.legend(*ax.get_legend_handles_labels(), 'center', prop = {'size': 28}, ncol=2, handlelength=2)
lgnd.legendHandles[0]._legmarker.set_markersize(20)
lgnd.legendHandles[1]._legmarker.set_markersize(20)
ax.set_xticks(xt)
ax.set_xticklabels(xt)
ax.set_ylim([10000, 20000])
ax.set_xlabel("days",fontsize = 20)
ax.set_ylabel("cell number",fontsize = 20)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
fig.savefig('solid05_time.pdf')
figlegend.savefig('legend.pdf')
plt.close()

s0_c1_nt       = [84, 84*2, 84*4, 84*8]
s0_c5_nt       = [13, 26, 52, 104]
s0_c1_dof  = np.zeros(1800*4*4 + 1)
ix = 144
nt = 84*4
s0_c1_dof  = getdof('../solid0./rf1/plots_solid_'+str( ix )+'_rf1_'+str(nt)+'_0./nx'+str( ix )+'_rf1.log', nt)
t1 = np.linspace(0, 12, nt + 1)


s0_c5_dof  = np.zeros(240*4*4 + 1)
ix = 144
nt = 52
s0_c5_dof  = getdof('../solid0./rf1/plots_solid_'+str( ix )+'_rf1_'+str(nt)+'_0./nx'+str( ix )+'_rf1.log', nt)
t5 = np.linspace(0, 12, nt + 1)

# print(len(t))
xt = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
fig = plt.figure(1, figsize=(16, 12))
plt.clf()
ax = fig.add_subplot(111)
ax.plot(t1, s0_c1_dof, 'k-', linewidth=4.0, label = 'Maximum Courant number around 1')
ax.plot(t5, s0_c5_dof, 'k:', linewidth=4.0, label = 'Maximum Courant number around 1')

# plt.legend(loc = 'best')
ax.set_xticks(xt)
ax.set_xticklabels(xt)
ax.set_ylim([10000, 20000])
ax.set_xlabel("days",fontsize = 20)
ax.set_ylabel("cell number",fontsize = 20)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
fig.savefig('solid0_time.pdf')
plt.close()