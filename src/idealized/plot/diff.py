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

nx             = np.array([72, 144, 288, 576])
dof            = np.array([2592, 10368, 41472, 165888])

m_c1_nt       = [1320, 1320*4, 1320*16, 1320*64]
m_c1_rf0_dof  = np.zeros(4)
m_c1_rf0_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c1_nt):
    m_c1_rf0_dof[i], m_c1_rf0_error[i, :]  = getdof('../moving0.25/rf0/plots_moving_'+str( ix )+'_rf0_'+str(nt)+'_0.25/nx'+str( ix )+'_rf0.log', nt)

m_c1_rf1_dof  = np.zeros(4)
m_c1_rf1_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c1_nt):
    m_c1_rf1_dof[i], m_c1_rf1_error[i, :]  = getdof('../moving0.25interp/rf1/plots_moving_'+str(ix//2)+'_rf1_'+str(nt)+'_0.25/nx'+str(ix//2)+'_rf1.log', nt)

m_c1_rf2_dof  = np.zeros(4)
m_c1_rf2_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c1_nt):
    m_c1_rf2_dof[i], m_c1_rf2_error[i, :]  = getdof('../moving0.25interp/rf2/plots_moving_'+str(ix//4)+'_rf2_'+str(nt)+'_0.25/nx'+str(ix//4)+'_rf2.log', nt)

m_uni_c1_rf1_dof  = np.zeros(4)
m_uni_c1_rf1_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c1_nt):
    m_uni_c1_rf1_dof[i], m_uni_c1_rf1_error[i, :]  = getdof('../moving0.25interpuniform/rf1/plots_moving_'+str(ix//2)+'_rf1_'+str(nt)+'_0.25/nx'+str(ix//2)+'_rf1.log', nt)  
        
m_uni_c1_rf2_dof  = np.zeros(4)
m_uni_c1_rf2_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c1_nt):
    m_uni_c1_rf2_dof[i], m_uni_c1_rf2_error[i, :]  = getdof('../moving0.25interpuniform/rf2/plots_moving_'+str(ix//4)+'_rf2_'+str(nt)+'_0.25/nx'+str(ix//4)+'_rf2.log', nt)
#------------------------------------------------------------------------------------------------------------------------------------------------------------
m_noint_c1_rf1_dof  = np.zeros(4)
m_noint_c1_rf1_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c1_nt):
    m_noint_c1_rf1_dof[i], m_noint_c1_rf1_error[i, :]  = getdof('../moving0.25/rf1/plots_moving_'+str(ix//2)+'_rf1_'+str(nt)+'_0.25/nx'+str(ix//2)+'_rf1.log', nt)

m_noint_c1_rf2_dof  = np.zeros(4)
m_noint_c1_rf2_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c1_nt):
    m_noint_c1_rf2_dof[i], m_noint_c1_rf2_error[i, :]  = getdof('../moving0.25/rf2/plots_moving_'+str(ix//4)+'_rf2_'+str(nt)+'_0.25/nx'+str(ix//4)+'_rf2.log', nt)

m_noint_uni_c1_rf1_dof  = np.zeros(4)
m_noint_uni_c1_rf1_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c1_nt):
    m_noint_uni_c1_rf1_dof[i], m_noint_uni_c1_rf1_error[i, :]  = getdof('../moving0.25uniform/rf1/plots_moving_'+str(ix//2)+'_rf1_'+str(nt)+'_0.25/nx'+str(ix//2)+'_rf1.log', nt)  
        
m_noint_uni_c1_rf2_dof  = np.zeros(4)
m_noint_uni_c1_rf2_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c1_nt):
    m_noint_uni_c1_rf2_dof[i], m_noint_uni_c1_rf2_error[i, :]  = getdof('../moving0.25uniform/rf2/plots_moving_'+str(ix//4)+'_rf2_'+str(nt)+'_0.25/nx'+str(ix//4)+'_rf2.log', nt)



m_c5_nt       = [240, 240*4, 240*16, 240*64]
m_c5_rf0_dof  = np.zeros(4)
m_c5_rf0_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c5_nt):
    m_c5_rf0_dof[i], m_c5_rf0_error[i, :]  = getdof('../moving0.25/rf0/plots_moving_'+str( ix )+'_rf0_'+str(nt)+'_0.25/nx'+str( ix )+'_rf0.log', nt)

m_c5_rf1_dof  = np.zeros(4)
m_c5_rf1_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c5_nt):
    m_c5_rf1_dof[i], m_c5_rf1_error[i, :]  = getdof('../moving0.25interp/rf1/plots_moving_'+str(ix//2)+'_rf1_'+str(nt)+'_0.25/nx'+str(ix//2)+'_rf1.log', nt)

m_c5_rf2_dof  = np.zeros(4)
m_c5_rf2_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c5_nt):
    m_c5_rf2_dof[i], m_c5_rf2_error[i, :]  = getdof('../moving0.25interp/rf2/plots_moving_'+str(ix//4)+'_rf2_'+str(nt)+'_0.25/nx'+str(ix//4)+'_rf2.log', nt)

m_uni_c5_rf1_dof  = np.zeros(4)
m_uni_c5_rf1_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c5_nt):
    m_uni_c5_rf1_dof[i], m_uni_c5_rf1_error[i, :]  = getdof('../moving0.25interpuniform/rf1/plots_moving_'+str(ix//2)+'_rf1_'+str(nt)+'_0.25/nx'+str(ix//2)+'_rf1.log', nt)

m_uni_c5_rf2_dof  = np.zeros(4)
m_uni_c5_rf2_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c5_nt):
    m_uni_c5_rf2_dof[i], m_uni_c5_rf2_error[i, :]  = getdof('../moving0.25interpuniform/rf2/plots_moving_'+str(ix//4)+'_rf2_'+str(nt)+'_0.25/nx'+str(ix//4)+'_rf2.log', nt)
#------------------------------------------------------------------------------------------------------------------------------------------------------
m_noint_c5_rf1_dof  = np.zeros(4)
m_noint_c5_rf1_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c5_nt):
    m_noint_c5_rf1_dof[i], m_noint_c5_rf1_error[i, :]  = getdof('../moving0.25/rf1/plots_moving_'+str(ix//2)+'_rf1_'+str(nt)+'_0.25/nx'+str(ix//2)+'_rf1.log', nt)

m_noint_c5_rf2_dof  = np.zeros(4)
m_noint_c5_rf2_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c5_nt):
    m_noint_c5_rf2_dof[i], m_noint_c5_rf2_error[i, :]  = getdof('../moving0.25/rf2/plots_moving_'+str(ix//4)+'_rf2_'+str(nt)+'_0.25/nx'+str(ix//4)+'_rf2.log', nt)

m_noint_uni_c5_rf1_dof  = np.zeros(4)
m_noint_uni_c5_rf1_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c5_nt):
    m_noint_uni_c5_rf1_dof[i], m_noint_uni_c5_rf1_error[i, :]  = getdof('../moving0.25uniform/rf1/plots_moving_'+str(ix//2)+'_rf1_'+str(nt)+'_0.25/nx'+str(ix//2)+'_rf1.log', nt)

m_noint_uni_c5_rf2_dof  = np.zeros(4)
m_noint_uni_c5_rf2_error  = np.zeros([4, 3])
for i, ix, nt in zip(range(4), nx, m_c5_nt):
    m_noint_uni_c5_rf2_dof[i], m_noint_uni_c5_rf2_error[i, :]  = getdof('../moving0.25uniform/rf2/plots_moving_'+str(ix//4)+'_rf2_'+str(nt)+'_0.25/nx'+str(ix//4)+'_rf2.log', nt)


dx = 360./nx
fig = plt.figure(1, figsize=(9, 6.75))
figlegend = plt.figure(figsize = (20, 20))
plt.clf()
ax = fig.add_subplot(111)
ax.loglog(dx[:-1], abs(m_c1_rf0_error[:-1,  1] - m_noint_c1_rf1_error[:-1,  1] )           , marker = '^', fillstyle = 'none',markersize = 15,  linestyle = '-',  color = 'k', label = 'added error from 1 level adaptive refinement' )
ax.loglog(dx[:],   abs(m_c1_rf0_error[:,  1] - m_noint_c1_rf2_error[:,  1]     )           , marker = 'o', fillstyle = 'none',markersize = 15,  linestyle = '-', color = 'k',  label = 'added error from 2 level adaptive refinement' )

ax.loglog(dx[:-1], abs(m_c1_rf0_error[:-1,  1] - m_noint_uni_c1_rf1_error[:-1,  1])    , marker = '^',markersize = 15,  linestyle = '-',  color = 'k',                         label = 'added error from 1 level uniform refinement (coarse initial condition)' )
ax.loglog(dx[:],   abs(m_c1_rf0_error[:,  1] - m_noint_uni_c1_rf2_error[:,  1])        , marker = 'o',markersize = 15,  linestyle = '-', color = 'k',                          label = 'added error from 2 level uniform refinement (coarse initial condition)' )

lgnd = figlegend.legend(*ax.get_legend_handles_labels(), 'center', prop = {'size': 28}, ncol=2, handlelength=2)
lgnd.legendHandles[0]._legmarker.set_markersize(20)
lgnd.legendHandles[1]._legmarker.set_markersize(20)
lgnd.legendHandles[2]._legmarker.set_markersize(20)
lgnd.legendHandles[3]._legmarker.set_markersize(20)
# plt.legend(loc = 'best')
ax.set_xticks(dx)
ax.set_xticklabels(dx)
ax.set_xlim([0.5, 6])
ax.set_ylim([5*10e-7, 1])
ax.set_xlabel("maximum resolution",fontsize = 20)
ax.set_ylabel(r'$\ell_2$',fontsize = 20)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
fig.savefig('l2_moving_c1_diff.pdf')
figlegend.savefig('legend_diff.pdf')
plt.close()

fig = plt.figure(2, figsize=(9, 6.75))
plt.clf()
ax = fig.add_subplot(111)
ax.loglog(dx[:-1], abs(m_c1_rf0_error[:-1,  2] - m_noint_c1_rf1_error[:-1,  2]) , marker = '^', fillstyle = 'none',markersize = 15,  linestyle = '-',  color = 'k', label = 'added error from 1 level adaptive refinement' )
ax.loglog(dx[:],   abs(m_c1_rf0_error[:,  2] - m_noint_c1_rf2_error[:,  2]    ) , marker = 'o', fillstyle = 'none',markersize = 15,  linestyle = '-', color = 'k',  label = 'added error from 2 level adaptive refinement' )
                   
ax.loglog(dx[:-1], abs(m_c1_rf0_error[:-1,  2] - m_noint_uni_c1_rf1_error[:-1,  2]), marker = '^',markersize = 15,  linestyle = '-',  color = 'k',                  label = 'added error from 1 level uniform refinement (coarse initial condition)' )
ax.loglog(dx[:],   abs(m_c1_rf0_error[:,  2] - m_noint_uni_c1_rf2_error[:,  2])    , marker = 'o',markersize = 15,  linestyle = '-', color = 'k',                   label = 'added error from 2 level uniform refinement (coarse initial condition)' )

# plt.legend(loc = 'best')
ax.set_xticks(dx)
ax.set_xticklabels(dx)
ax.set_xlim([0.5, 6])
ax.set_ylim([5*10e-7, 1])
ax.set_xlabel("maximum resolution",fontsize = 20)
ax.set_ylabel(r'$\ell_\infty$',fontsize = 20)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
plt.savefig('linf_moving_c1_diff.pdf')
plt.close()

print(m_c1_rf0_error[:,  2])
print(m_c1_rf1_error[:-1,  2])
print(m_noint_c1_rf1_error[:-1,  2])
print(m_c1_rf1_error[:-1,  2] - m_noint_c1_rf1_error[:-1,  2])

fig = plt.figure(3, figsize=(9, 6.75))
plt.clf()
ax = fig.add_subplot(111)
ax.loglog(dx[:], abs(m_c5_rf0_error[:,  1]   - m_noint_c5_rf1_error[:,  1] ) , marker = '^', fillstyle = 'none',markersize = 15,  linestyle = '-',  color = 'k', label = 'added error from 1 level adaptive refinement' )
ax.loglog(dx[:], abs(m_c5_rf0_error[:,  1]   - m_noint_c5_rf2_error[:,  1] ) , marker = 'o', fillstyle = 'none',markersize = 15,  linestyle = '-', color = 'k',  label = 'added error from 2 level adaptive refinement' )

ax.loglog(dx[:], abs(m_c5_rf0_error[:,  1] - m_noint_uni_c5_rf1_error[:,  1]), marker = '^', markersize = 15, linestyle = '-',  color = 'k',                     label = 'added error from 1 level uniform refinement (coarse initial condition)' )
ax.loglog(dx[:], abs(m_c5_rf0_error[:,  1] - m_noint_uni_c5_rf2_error[:,  1]), marker = 'o',markersize = 15,  linestyle = '-', color = 'k',                      label = 'added error from 2 level uniform refinement (coarse initial condition)' )

# plt.legend(loc = 'best')
ax.set_xticks(dx)
ax.set_xticklabels(dx)
ax.set_xlim([0.5, 6])
ax.set_ylim([5*10e-7, 1])
ax.set_xlabel("maximum resolution",fontsize = 20)
ax.set_ylabel(r'$\ell_2$',fontsize = 20)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
plt.savefig('l2_moving_c5_diff.pdf')
plt.close()


fig = plt.figure(4, figsize=(9, 6.75))
plt.clf()
ax = fig.add_subplot(111)
ax.loglog(dx[:], abs(m_c5_rf0_error[:,  2]    - m_noint_c5_rf1_error[:,  2]) , marker = '^', fillstyle = 'none',markersize = 15,  linestyle = '-',  color = 'k', label = 'added error from 1 level adaptive refinement' )
ax.loglog(dx[:], abs(m_c5_rf0_error[:,  2]    - m_noint_c5_rf2_error[:,  2]) , marker = 'o', fillstyle = 'none',markersize = 15,  linestyle = '-', color = 'k',  label = 'added error from 2 level adaptive refinement' )

ax.loglog(dx[:], abs(m_c5_rf0_error[:,  2] - m_noint_uni_c5_rf1_error[:,  2]), marker = '^', markersize = 15, linestyle = '-',  color = 'k',                     label = 'added error from 1 level uniform refinement (coarse initial condition)' )
ax.loglog(dx[:], abs(m_c5_rf0_error[:,  2] - m_noint_uni_c5_rf2_error[:,  2]), marker = 'o',markersize = 15,  linestyle = '-', color = 'k',                      label = 'added error from 2 level uniform refinement (coarse initial condition)' )
# plt1legend(loc = 'best')
ax.set_xticks(dx)
ax.set_xticklabels(dx)
ax.set_xlim([0.5, 6])
ax.set_ylim([5*10e-7, 1])
ax.set_xlabel("maximum resolution",fontsize = 20)
ax.set_ylabel(r'$\ell_\infty$',fontsize = 20)
ax.tick_params(axis="x", labelsize=20)
ax.tick_params(axis="y", labelsize=20)
plt.savefig('linf_moving_c5_diff.pdf')
plt.close()