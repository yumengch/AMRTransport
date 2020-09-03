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
    i = 0
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
    return mass

nx             = np.array([72, 144, 288, 576])
dof            = np.array([2592, 10368, 41472, 165888])

d_c1_nt       = [120, 240, 480, 960]
d_c1_rf0_dof  = np.zeros(4)
d_c1_rf0_error  = np.zeros([4, 3])
# for i, ix, nt in zip(range(4), nx, d_c1_nt):
mass  = getdof('../divergent0.5/rf0/plots_divergent_576_rf0_960_0.5/nx576_rf0.log', 960)
meanMass = np.mean(mass)
mass = (mass - meanMass)/meanMass
print(np.max(mass), np.min(mass))
fig = plt.figure(1, figsize=(16, 12))
plt.clf()
ax = fig.add_subplot(111)
ax.plot(range(0, 961), mass, color = 'k', linewidth = 2.0)
ax.set_ylim([-1e-13, 1e-13])
ax.set_xlabel("steps",fontsize = 35)
ax.set_ylabel("mass",fontsize = 35)
ax.tick_params(axis="x", labelsize=35)
ax.tick_params(axis="y", labelsize=35)
text = ax.yaxis.get_offset_text() # Get the text object
text.set_size(35) # # Set the size.
fig.savefig("conserverf0.pdf")
plt.close()

mass  = getdof('../divergent0.5/rf1/plots_divergent_288_rf1_960_0.5/nx288_rf1.log', 960)
meanMass = np.mean(mass)
mass = (mass - meanMass)/meanMass
print(np.max(mass), np.min(mass))
fig = plt.figure(1, figsize=(16, 12))
plt.clf()
ax = fig.add_subplot(111)
ax.plot(range(0, 961), mass, color = 'k', linewidth = 2.0)
ax.set_ylim([-0.4*1e-11, 0.4*1e-11])
ax.set_xlabel("steps",fontsize = 35)
ax.set_ylabel("mass",fontsize = 35)
ax.tick_params(axis="x", labelsize=35)
ax.tick_params(axis="y", labelsize=35)
text = ax.yaxis.get_offset_text() # 35 the text object
text.set_size(35) # # Set the size.
fig.savefig("conserverf1.pdf")
plt.close()