import numpy as np
import matplotlib.pyplot as plt
# compute cell numbers (degrees of freedom)
def getdof(filename, nt):
    f = open(filename, 'r')
    print(filename)
    val = 0
    error = np.zeros([nt + 1, 3])
    mass  = np.zeros([nt + 1])
    line = f.readline()
    while line:
        string = line.split()
        # print(string)
        if len(string) >= 2:
            if string[0] == 'nColumns':
                val += int(string[-1])
            elif string[1] == 'time':
                error[int(string[2]), 0] = float(string[5])
                error[int(string[2]), 1] = float(string[8]) 
                error[int(string[2]), 2] = float(string[11])
                mass[int(string[2])]     = float(string[-1])
        line = f.readline()
    f.close()
    return val/(nt + 1), error[-1, :]

nx             = np.array([72, 144, 288, 576], dtype = 'int32')
dof            = np.array([2592, 10368, 41472, 165888])

d_c5_nt       = [24, 48, 96, 192]
d_c5_rf0_dof  = np.zeros(4)
d_c5_rf0_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, d_c5_nt):
    d_c5_rf0_dof[i], d_c5_rf0_error[i, :]  = getdof('../divergent0.5/rf0/plots_divergent_'+str( ix )+'_rf0_'+str(nt)+'_0.5/nx'+str( ix )+'_rf0.log', nt)

d_c5_rf1_dof  = np.zeros(4)
d_c5_rf1_error  = np.zeros([4, 3])
for i, il, nt in zip(range(4), range(2,6), d_c5_nt):
    d_c5_rf1_dof[i], d_c5_rf1_error[i, :]  = getdof('../divergent0.5/rf'+str(il)+'/plots_divergent_18_rf'+str(il)+'_' + str(nt) +'_0.5/nx18_rf' +str(il)+'.log', nt)
    
fig = plt.figure(num=1, figsize=(9, 6.75))
figlegend = plt.figure(figsize = (15, 15))
plt.clf()
ax = fig.add_subplot(111)
ax.loglog(dof[:], d_c5_rf0_error[:,  1]         , marker = 's', markersize=15, linestyle = '-',  color = 'k', label = '0 level refinement')
ax.loglog(d_c5_rf1_dof[:], d_c5_rf1_error[:,  1], marker = '^', markersize=15, linestyle = '-',  color = 'k', label = 'adaptive')

ax.loglog(np.array([dof[0],dof[-1]])//2, \
    1.6*np.array([d_c5_rf0_error[0, 1], d_c5_rf0_error[0, 1]*(float(nx[0])/nx[-1])**1]), 'k--', label='1st order')
ax.loglog(np.array([dof[0],dof[-1]])//2, \
    1.4*np.array([d_c5_rf0_error[0, 1], d_c5_rf0_error[0, 1]*(float(nx[0])/nx[-1])**2]), 'k-.', label='2nd order')
# ax.loglog(np.array([dof[0],dof[-1]])//2, \
#     1.2*np.array([d_c5_rf0_error[0, 1], d_c5_rf0_error[0, 1]*(float(nx[0])/nx[-1])**3]), 'k:', label='3rd order')
lgnd = figlegend.legend(*ax.get_legend_handles_labels(), 'center', prop = {'size': 28}, ncol=2, handlelength=2)
lgnd.legendHandles[0]._legmarker.set_markersize(20)
lgnd.legendHandles[1]._legmarker.set_markersize(20)
lgnd.legendHandles[2]._legmarker.set_markersize(20)
lgnd.legendHandles[3]._legmarker.set_markersize(20)
ax.set_xticks(dof)
ax.set_xticklabels(dof)
ax.set_xlim([500, 200000])
ax.set_ylim([5*10e-6, 1])
ax.set_xlabel("cell number",fontsize = 20)
ax.set_ylabel(r'$\ell_2$',fontsize = 20)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
fig.savefig('l2_div_c5_eff_multi_ref.pdf')
figlegend.savefig('legend.pdf')
plt.close()

fig = plt.figure(2, figsize=(9, 6.75))
plt.clf()
ax = fig.add_subplot(111)
plt.loglog(dof[:], d_c5_rf0_error[:,  2]         , marker = 's', markersize=15, linestyle = '-',  color = 'k', label = '0 level refinement')
plt.loglog(d_c5_rf1_dof[:], d_c5_rf1_error[:,  2], marker = '^', markersize=15, linestyle = '-',  color = 'k', label = 'adaptive')

plt.loglog(np.array([dof[0],dof[-1]])//2, \
    1.6*np.array([d_c5_rf0_error[0, 2], d_c5_rf0_error[0, 2]*(float(nx[0])/nx[-1])**1]), 'k--', label='1st order')
plt.loglog(np.array([dof[0],dof[-1]])//2, \
    1.4*np.array([d_c5_rf0_error[0, 2], d_c5_rf0_error[0, 2]*(float(nx[0])/nx[-1])**2]), 'k-.', label='2nd order')
# plt.loglog(np.array([dof[0],dof[-1]])//2, \
#     1.2*np.array([d_c5_rf0_error[0, 2], d_c5_rf0_error[0, 2]*(float(nx[0])/nx[-1])**3]), 'k:', label='3rd order')
# plt.legend(loc = 'best')
ax.set_xticks(dof)
ax.set_xticklabels(dof)
ax.set_xlim([500, 200000])
ax.set_ylim([5*10e-6, 1])
ax.set_xlabel("cell number",fontsize=20)
ax.set_ylabel(r'$\ell_\infty$',fontsize=20)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
fig.savefig('linf_div_c5_eff_multi_ref.pdf')
plt.close()