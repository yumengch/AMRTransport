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
                val += int(string[2])
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
s0_c1_nt       = np.array([84, 168, 336, 672])

s0_c1_rf0_dof  = np.zeros(4)
s0_c1_rf0_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, s0_c1_nt):
    s0_c1_rf0_dof[i], s0_c1_rf0_error[i, :]  = getdof('../solid0./rf0/plots_solid_'+str( ix )+'_rf0_'+str(nt)+'_0./nx'+str( ix )+'_rf0.log', nt)

s0_c1_rf1_dof  = np.zeros(4)
s0_c1_rf1_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, s0_c1_nt):
    s0_c1_rf1_dof[i], s0_c1_rf1_error[i, :]  = getdof('../solid0./rf1/plots_solid_'+str(ix//2)+'_rf1_'+str(nt)+'_0./nx'+str(ix//2)+'_rf1.log', nt)

s0_c1_rf2_dof  = np.zeros(4)
s0_c1_rf2_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, s0_c1_nt):
    s0_c1_rf2_dof[i], s0_c1_rf2_error[i, :]  = getdof('../solid0./rf2/plots_solid_'+str(ix//4)+'_rf2_'+str(nt)+'_0./nx'+str(ix//4)+'_rf2.log', nt)

s0_c5_nt       = [13, 26, 52, 104]
s0_c5_rf0_dof  = np.zeros(4)
s0_c5_rf0_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, s0_c5_nt):
    s0_c5_rf0_dof[i], s0_c5_rf0_error[i, :]  = getdof('../solid0./rf0/plots_solid_'+str( ix )+'_rf0_'+str(nt)+'_0./nx'+str( ix )+'_rf0.log', nt)

s0_c5_rf1_dof  = np.zeros(4)
s0_c5_rf1_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, s0_c5_nt):
    s0_c5_rf1_dof[i], s0_c5_rf1_error[i, :]  = getdof('../solid0./rf1/plots_solid_'+str(ix//2)+'_rf1_'+str(nt)+'_0./nx'+str(ix//2)+'_rf1.log', nt)

s0_c5_rf2_dof  = np.zeros(4)
s0_c5_rf2_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, s0_c5_nt):
    s0_c5_rf2_dof[i], s0_c5_rf2_error[i, :]  = getdof('../solid0./rf2/plots_solid_'+str(ix//4)+'_rf2_'+str(nt)+'_0./nx'+str(ix//4)+'_rf2.log', nt)


s05_c1_nt       = [1800, 1800*4, 1800*4*4, 1800*4*4*4]

s05_c1_rf0_dof  = np.zeros(4)
s05_c1_rf0_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, s05_c1_nt):
    s05_c1_rf0_dof[i], s05_c1_rf0_error[i, :]  = getdof('../solid0.5/rf0/plots_solid_'+str( ix )+'_rf0_'+str(nt)+'_0.5/nx'+str( ix )+'_rf0.log', nt)

s05_c1_rf1_dof  = np.zeros(4)
s05_c1_rf1_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, s05_c1_nt):
    s05_c1_rf1_dof[i], s05_c1_rf1_error[i, :]  = getdof('../solid0.5/rf1/plots_solid_'+str(ix//2)+'_rf1_'+str(nt)+'_0.5/nx'+str(ix//2)+'_rf1.log', nt)

s05_c1_rf2_dof  = np.zeros(4)
s05_c1_rf2_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, s05_c1_nt):
    s05_c1_rf2_dof[i], s05_c1_rf2_error[i, :]  = getdof('../solid0.5/rf2/plots_solid_'+str(ix//4)+'_rf2_'+str(nt)+'_0.5/nx'+str(ix//4)+'_rf2.log', nt)

s05_c5_nt       = [240, 240*4, 240*16, 240*64]
s05_c5_rf0_dof  = np.zeros(4)
s05_c5_rf0_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, s05_c5_nt):
    s05_c5_rf0_dof[i], s05_c5_rf0_error[i, :]  = getdof('../solid0.5/rf0/plots_solid_'+str( ix )+'_rf0_'+str(nt)+'_0.5/nx'+str( ix )+'_rf0.log', nt)

s05_c5_nt       = [240, 240*4, 240*16, 240*64]
s05_c5_rf1_dof  = np.zeros(4)
s05_c5_rf1_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, s05_c5_nt):
    s05_c5_rf1_dof[i], s05_c5_rf1_error[i, :]  = getdof('../solid0.5/rf1/plots_solid_'+str(ix//2)+'_rf1_'+str(nt)+'_0.5/nx'+str(ix//2)+'_rf1.log', nt)

s05_c5_nt       = [240, 240*4, 240*16, 240*64]
s05_c5_rf2_dof  = np.zeros(4)
s05_c5_rf2_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, s05_c5_nt):
    s05_c5_rf2_dof[i], s05_c5_rf2_error[i, :]  = getdof('../solid0.5/rf2/plots_solid_'+str(ix//4)+'_rf2_'+str(nt)+'_0.5/nx'+str(ix//4)+'_rf2.log', nt)

dx = 360./nx

fig = plt.figure(num=1, figsize=(9, 6.75))
figlegend = plt.figure(figsize = (20, 20))
plt.clf()
ax = fig.add_subplot(111)
ax.loglog(dof, s0_c1_rf0_error[:, 1]               , marker = 's', markersize=15, linestyle = '-',  color = 'k', label = '0 level refinement, '+r'$\alpha = 0$')
ax.loglog(s0_c1_rf1_dof, s0_c1_rf1_error[:, 1]     , marker = '^', markersize=15, linestyle = '-',  color = 'k', label = '1 level refinement, '+r'$\alpha = 0$')
ax.loglog(s0_c1_rf2_dof, s0_c1_rf2_error[:, 1]     , marker = 'o', markersize=15, linestyle = '-', color = 'k', label = '2 level refinement, '+r'$\alpha = 0$')

ax.loglog(dof[:], s05_c1_rf0_error[:, 1]           , marker = 's', fillstyle = 'none', markersize=15, linestyle = '-',  color = 'k', label = '0 level refinement, '+r'$\alpha = 0.5\pi$')
ax.loglog(s05_c1_rf1_dof[:], s05_c1_rf1_error[:, 1], marker = '^', fillstyle = 'none', markersize=15, linestyle = '-',  color = 'k', label = '1 level refinement, '+r'$\alpha = 0.5\pi$')
ax.loglog(s05_c1_rf2_dof[:], s05_c1_rf2_error[:, 1], marker = 'o', fillstyle = 'none', markersize=15, linestyle = '-', color = 'k', label = '2 level refinement, '+r'$\alpha = 0.5\pi$')

# ax.loglog(np.array([dof[0],dof[-1]]), \
#     1.6*np.array([s0_c5_rf0_error[0, 1], s0_c5_rf0_error[0, 1]*(float(nx[0])/nx[-1])**1]), 'k--', label='1st order')
ax.loglog(np.array([dof[0],dof[-1]]), \
   1.4*np.array([s0_c5_rf0_error[0, 1], s0_c5_rf0_error[0, 1]*(float(nx[0])/nx[-1])**2]), 'k-.', label='2nd order')
ax.loglog(np.array([dof[0],dof[-1]]), \
    1.2*np.array([s0_c5_rf0_error[0, 1], s0_c5_rf0_error[0, 1]*(float(nx[0])/nx[-1])**3]), 'k:', label='3rd order')
# plt.legend(loc = 'best')
lgnd = figlegend.legend(*ax.get_legend_handles_labels(), 'center', prop = {'size': 28}, ncol=3, handlelength=2)
lgnd.legendHandles[0]._legmarker.set_markersize(20)
lgnd.legendHandles[1]._legmarker.set_markersize(20)
lgnd.legendHandles[2]._legmarker.set_markersize(20)
lgnd.legendHandles[3]._legmarker.set_markersize(20)
lgnd.legendHandles[4]._legmarker.set_markersize(20)
lgnd.legendHandles[5]._legmarker.set_markersize(20)
lgnd.legendHandles[6]._legmarker.set_markersize(20)
lgnd.legendHandles[7]._legmarker.set_markersize(20)
# lgnd.legendHandles[8]._legmarker.set_markersize(20)
ax.set_xticks(dof)
ax.set_xticklabels(dof)
ax.set_xlim([500, 200000])
ax.set_ylim([5*10e-6, 1])
ax.set_xlabel("cell number",fontsize = 20)
ax.set_ylabel(r'$\ell_2$',fontsize = 20)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
fig.savefig('l2_solid_c1_eff.pdf')
figlegend.savefig('legend.pdf')
plt.close()

fig = plt.figure(2, figsize=(9, 6.75))
plt.clf()
ax = fig.add_subplot(111)
ax.loglog(dof, s0_c1_rf0_error[:, 2]               , marker = 's', markersize=15, linestyle = '-',  color = 'k', label = '0 level refinement, '+r'$\alpha = 0$')
ax.loglog(s0_c1_rf1_dof, s0_c1_rf1_error[:, 2]     , marker = '^', markersize=15, linestyle = '-',  color = 'k', label = '1 level refinement, '+r'$\alpha = 0$')
ax.loglog(s0_c1_rf2_dof, s0_c1_rf2_error[:, 2]     , marker = 'o', markersize=15, linestyle = '-', color = 'k', label = '2 level refinement, '+r'$\alpha = 0$')

ax.loglog(dof[:], s05_c1_rf0_error[:, 2]           , marker = 's', fillstyle = 'none', markersize=15, linestyle = '-',  color = 'k', label = '0 level refinement, '+r'$\alpha = 0.5\pi$')
ax.loglog(s05_c1_rf1_dof[:], s05_c1_rf1_error[:, 2], marker = '^', fillstyle = 'none', markersize=15, linestyle = '-',  color = 'k', label = '1 level refinement, '+r'$\alpha = 0.5\pi$')
ax.loglog(s05_c1_rf2_dof[:], s05_c1_rf2_error[:, 2], marker = 'o', fillstyle = 'none', markersize=15, linestyle = '-', color = 'k', label = '2 level refinement, '+r'$\alpha = 0.5\pi$')

# ax.loglog(np.array([dof[0],dof[-1]])//2, \
#     1.6*np.array([s0_c5_rf0_error[0, 1], s0_c5_rf0_error[0, 1]*(float(nx[0])/nx[-1])**1]), 'k--', label='1st order')
ax.loglog(np.array([dof[0],dof[-1]])//2, \
    1.4*np.array([s0_c5_rf0_error[0, 1], s0_c5_rf0_error[0, 1]*(float(nx[0])/nx[-1])**2]), 'k-.', label='2nd order')
ax.loglog(np.array([dof[0],dof[-1]])//2, \
    1.2*np.array([s0_c5_rf0_error[0, 1], s0_c5_rf0_error[0, 1]*(float(nx[0])/nx[-1])**3]), 'k:', label='3rd order')
# ax.legend(loc = 'best')
ax.set_xticks(dof)
ax.set_xticklabels(dof)
ax.set_xlim([500, 200000])
ax.set_ylim([5*10e-6, 1])
ax.set_xlabel("cell number",fontsize=20)
ax.set_ylabel(r'$\ell_\infty$',fontsize=20)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
fig.savefig('linf_solid_c1_eff.pdf')
plt.close()

fig = plt.figure(3, figsize=(9, 6.75))
plt.clf()
ax = fig.add_subplot(111)
ax.loglog(dof, s0_c5_rf0_error[:, 1]               , marker = 's', markersize=15, linestyle = '-',  color = 'k', label = '0 level refinement, '+r'$\alpha = 0$')
ax.loglog(s0_c5_rf1_dof, s0_c5_rf1_error[:, 1]     , marker = '^', markersize=15, linestyle = '-',  color = 'k', label = '1 level refinement, '+r'$\alpha = 0$')
ax.loglog(s0_c5_rf2_dof, s0_c5_rf2_error[:, 1]     , marker = 'o', markersize=15, linestyle = '-', color = 'k', label = '2 level refinement, '+r'$\alpha = 0$')

ax.loglog(dof[:], s05_c5_rf0_error[:, 1]           , marker = 's', fillstyle = 'none', markersize=15, linestyle = '-',  color = 'k', label = '0 level refinement, '+r'$\alpha = 0.5\pi$')
ax.loglog(s05_c5_rf1_dof[:], s05_c5_rf1_error[:, 1], marker = '^', fillstyle = 'none', markersize=15, linestyle = '-',  color = 'k', label = '1 level refinement, '+r'$\alpha = 0.5\pi$')
ax.loglog(s05_c5_rf2_dof[:], s05_c5_rf2_error[:, 1], marker = 'o', fillstyle = 'none', markersize=15, linestyle = '-', color = 'k', label = '2 level refinement, '+r'$\alpha = 0.5\pi$')

# ax.loglog(np.array([dof[0],dof[-1]])//2, \
#     1.6*np.array([s0_c5_rf0_error[0, 1], s0_c5_rf0_error[0, 1]*(float(nx[0])/nx[-1])**1]), 'k--', label='1st order')
ax.loglog(np.array([dof[0],dof[-1]])//2, \
    1.4*np.array([s0_c5_rf0_error[0, 1], s0_c5_rf0_error[0, 1]*(float(nx[0])/nx[-1])**2]), 'k-.', label='2nd order')
ax.loglog(np.array([dof[0],dof[-1]])//2, \
    1.2*np.array([s0_c5_rf0_error[0, 1], s0_c5_rf0_error[0, 1]*(float(nx[0])/nx[-1])**3]), 'k:', label='3rd order')
# ax.legend(loc = 'best')
ax.set_xticks(dof)
ax.set_xticklabels(dof)
ax.set_xlim([500, 200000])
ax.set_ylim([5*10e-6, 1])
ax.set_xlabel("cell number",fontsize = 20)
ax.set_ylabel(r'$\ell_2$',fontsize = 20)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
fig.savefig('l2_solid_c5_eff.pdf')
plt.close()


fig = plt.figure(4, figsize=(9, 6.75))
plt.clf()
ax = fig.add_subplot(111)
ax.loglog(dof, s0_c5_rf0_error[:, 2]               , marker = 's', markersize=15, linestyle = '-',  color = 'k', label = '0 level refinement, '+r'$\alpha = 0$')
ax.loglog(s0_c5_rf1_dof, s0_c5_rf1_error[:, 2]     , marker = '^', markersize=15, linestyle = '-',  color = 'k', label = '1 level refinement, '+r'$\alpha = 0$')
ax.loglog(s0_c5_rf2_dof, s0_c5_rf2_error[:, 2]     , marker = 'o', markersize=15, linestyle = '-', color = 'k', label = '2 level refinement, '+r'$\alpha = 0$')

ax.loglog(dof[:], s05_c5_rf0_error[:, 2]           , marker = 's', fillstyle = 'none', markersize=15, linestyle = '-',  color = 'k', label = '0 level refinement, '+r'$\alpha = 0.5\pi$')
ax.loglog(s05_c5_rf1_dof[:], s05_c5_rf1_error[:, 2], marker = '^', fillstyle = 'none', markersize=15, linestyle = '-',  color = 'k', label = '1 level refinement, '+r'$\alpha = 0.5\pi$')
ax.loglog(s05_c5_rf2_dof[:], s05_c5_rf2_error[:, 2], marker = 'o', fillstyle = 'none', markersize=15, linestyle = '-', color = 'k', label = '2 level refinement, '+r'$\alpha = 0.5\pi$')

# ax.loglog(np.array([dof[0],dof[-1]])//2, \
#     1.6*np.array([s0_c5_rf0_error[0, 2], s0_c5_rf0_error[0, 2]*(float(nx[0])/nx[-1])**1]), 'k--', label='1st order')
ax.loglog(np.array([dof[0],dof[-1]])//2, \
    1.4*np.array([s0_c5_rf0_error[0, 2], s0_c5_rf0_error[0, 2]*(float(nx[0])/nx[-1])**2]), 'k-.', label='2nd order')
ax.loglog(np.array([dof[0],dof[-1]])//2, \
    1.2*np.array([s0_c5_rf0_error[0, 2], s0_c5_rf0_error[0, 2]*(float(nx[0])/nx[-1])**3]), 'k:', label='3rd order')
# ax.legend(loc = 'best')
ax.set_xticks(dof)
ax.set_xticklabels(dof)
ax.set_xlim([500, 200000])
ax.set_ylim([5*10e-6, 1])
ax.set_xlabel("cell number",fontsize=20)
ax.set_ylabel(r'$\ell_\infty$',fontsize=20)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
fig.savefig('linf_solid_c5_eff.pdf')
plt.close()