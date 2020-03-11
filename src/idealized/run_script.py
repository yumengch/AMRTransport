import subprocess
import os
import shutil
import logging
import numpy as np
def results(name_dir, name_log, name_vtu):
    test = os.listdir(os.curdir)
    for item in test:
        if item.endswith(".mod") or item.endswith(".o"):
            os.remove(item)

    logging.basicConfig(filename=name_log, level=logging.INFO, format='%(message)s')
    process = subprocess.Popen(["./FFSL_AMR"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    for line in iter(process.stdout.readline, b''): # b'\n'-separated lines
        logging.info(line)
    log = logging.getLogger()
    for hdlr in log.handlers[:]:                    # remove all old handlers
        log.removeHandler(hdlr)
    if not os.path.isdir(name_dir):
        subprocess.call(["mkdir", name_dir])
    file = 'AMRDUST.nc'
    shutil.move(os.path.join(os.getcwd(),file), os.path.join(name_dir,file))
    os.rename(name_log, name_dir+"/"+name_log)

def run_conf_write(nx, ny, nt, maxreflvls, beta, initial, nof, theta_r, theta_c, exp = None):
    f = open("test.nml", "w+")
    f.write("&grid\n")
    s = " nx = "+str(nx)
    f.write(s+",\n")
    s = " ny = "+str(ny)
    f.write(s+",\n")
    s = " nt = "+str(nt)
    f.write(s+",\n")
    s = "MaxRef = "+str(maxreflvls)
    f.write(s+",\n")
    s = "nLimit = "+str(maxreflvls)
    f.write(s+",\n")
    s = "theta_r = " + str(theta_r)
    f.write(s+",\n")
    s = "theta_c = " + str(theta_c)
    f.write(s+",\n")        
    f.write("/\n")
    f.write("&output\n")
    s = "nof = " + str(nof)
    f.write(s +", \n")
    f.write("/\n")
    f.write("&test\n")    
    s = " beta = " + beta
    f.write(s +",\n")
    s = "TestName = "
    f.write(s + initial + ",\n")
    f.write("/\n")
    f.close()
    base_dir = initial[1:-1]+beta
    if exp is not None:
        base_dir += exp    
    if not os.path.isdir(base_dir):
        subprocess.call(["mkdir", base_dir])
    base_dir += "/rf"+str(maxreflvls)
    if not os.path.isdir(base_dir):
        subprocess.call(["mkdir", base_dir])
    name_dir = base_dir + "/plots_" + initial[1:-1] + "_" +str(nx) + "_rf" + str(maxreflvls) \
        + "_"+ str(nt) + "_" +beta
    name_log = "nx"+str(nx)+"_rf"+str(maxreflvls) + ".log"
    name_vtu = initial[1:-1]+"_"
    results(name_dir, name_log, name_vtu)


subprocess.call(["bash", "run.sh"])
# # solid body rotation beta = 0.

# # c1
# rf1
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*7,   maxreflvls = 1, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)
run_conf_write(nx = 72,    ny = 36,    nt = 12*7*2, maxreflvls = 1, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*7*4, maxreflvls = 1, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)
run_conf_write(nx = 72*4,  ny = 36*4,  nt = 12*7*8, maxreflvls = 1, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)
# rf2
run_conf_write(nx = 72//4, ny = 36//4, nt = 12*7,   maxreflvls = 2, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*7*2, maxreflvls = 2, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)
run_conf_write(nx = 72  ,  ny = 36,    nt = 12*7*4, maxreflvls = 2, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*7*8, maxreflvls = 2, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)

# c5.5
# rf1
run_conf_write(nx = 72//2, ny = 36//2, nt = 13,     maxreflvls = 1, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)
run_conf_write(nx = 72,    ny = 36,    nt = 13*2,   maxreflvls = 1, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 13*4,   maxreflvls = 1, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)
run_conf_write(nx = 72*4,  ny = 36*4,  nt = 13*8,   maxreflvls = 1, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)

# rf2
run_conf_write(nx = 72//4, ny = 36//4, nt = 13,     maxreflvls = 2, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)
run_conf_write(nx = 72//2, ny = 36//2, nt = 13*2,   maxreflvls = 2, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)
run_conf_write(nx = 72,    ny = 36,    nt = 13*4,   maxreflvls = 2, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 13*8,   maxreflvls = 2, beta = "0.", initial = "\"solid\"", nof = "12", theta_r = 1e-6, theta_c = 1e-5)

# beta = 0.5pi

# # c1
# rf1
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*150,    maxreflvls = 1, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72,    ny = 36,    nt = 12*150*4,  maxreflvls = 1, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*150*16, maxreflvls = 1, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72*4,  ny = 36*4,  nt = 12*150*64, maxreflvls = 1, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
# # rf2
run_conf_write(nx = 72//4, ny = 36//4, nt = 12*150,    maxreflvls = 2, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*150*4,  maxreflvls = 2, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72,    ny = 36,    nt = 12*150*16, maxreflvls = 2, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*150*64, maxreflvls = 2, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
# # c6
# rf1
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*20,     maxreflvls = 1, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72,    ny = 36,    nt = 12*20*4,   maxreflvls = 1, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*20*16,  maxreflvls = 1, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72*4,  ny = 36*4,  nt = 12*20*64,  maxreflvls = 1, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)

# # rf2
run_conf_write(nx = 72//4, ny = 36//4, nt = 12*20,     maxreflvls = 2, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*20*4,   maxreflvls = 2, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72,    ny = 36,    nt = 12*20*16,  maxreflvls = 2, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*20*64,  maxreflvls = 2, beta = "0.5", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)


# intermediate
# rf1
run_conf_write(nx = 72//2, ny = 36//2, nt = 55,   maxreflvls = 1, beta = "0.15", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72,    ny = 36,    nt = 55*2, maxreflvls = 1, beta = "0.15", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 55*4, maxreflvls = 1, beta = "0.15", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72*4,  ny = 36*4,  nt = 55*8, maxreflvls = 1, beta = "0.15", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)

# rf2
run_conf_write(nx = 72//4, ny = 36//4, nt = 55,   maxreflvls = 2, beta = "0.15", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72//2, ny = 36//2, nt = 55*2, maxreflvls = 2, beta = "0.15", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72,    ny = 36,    nt = 55*4, maxreflvls = 2, beta = "0.15", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 55*8, maxreflvls = 2, beta = "0.15", initial = "\"solid\"", nof = "12", theta_r = 1e-3, theta_c = 1e-2)

# divergent test
# c1
# rf0
run_conf_write(nx = 72*8, ny = 36*8, nt = 12*10*8, maxreflvls = 0, beta = "0.5", initial = "\"divergent\"", nof = "12", theta_r = 0.2, theta_c = 0.15)
# # rf1
run_conf_write(nx = 72*4,  ny = 36*4,  nt = 12*10*8, maxreflvls = 1, beta = "0.5", initial = "\"divergent\"", nof = "12", theta_r = 0.2, theta_c = 0.15)
# c5
# rf0
run_conf_write(nx = 72,   ny = 36,   nt = 12*2,   maxreflvls = 0, beta = "0.5", initial = "\"divergent\"", nof = "12", theta_r = 0.2, theta_c = 0.15)
run_conf_write(nx = 72*2, ny = 36*2, nt = 12*2*2, maxreflvls = 0, beta = "0.5", initial = "\"divergent\"", nof = "12", theta_r = 0.2, theta_c = 0.15)
run_conf_write(nx = 72*4, ny = 36*4, nt = 12*2*4, maxreflvls = 0, beta = "0.5", initial = "\"divergent\"", nof = "12", theta_r = 0.2, theta_c = 0.15)
run_conf_write(nx = 72*8, ny = 36*8, nt = 12*2*8, maxreflvls = 0, beta = "0.5", initial = "\"divergent\"", nof = "12", theta_r = 0.2, theta_c = 0.15)
# multiple levels of refinement
for i in range(2, 6):
    run_conf_write(nx = 18, ny = 9, nt = 12*2**(i-1),   maxreflvls = i, beta = "0.5", initial = "\"divergent\"", nof = "12", theta_r = 0.2, theta_c = 0.15)

# moving (low initial condition + exact wind)
# c1
# rf0
run_conf_write(nx = 72,   ny = 36  , nt = 12*110,    maxreflvls = 0, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72*2, ny = 36*2, nt = 12*110*4,  maxreflvls = 0, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72*4, ny = 36*4, nt = 12*110*16, maxreflvls = 0, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72*8, ny = 36*8, nt = 12*110*64, maxreflvls = 0, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
# rf1
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*110,    maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72,    ny = 36,    nt = 12*110*4,  maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*110*16, maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72*4,  ny = 36*4,  nt = 12*110*64, maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
# # rf2
run_conf_write(nx = 72//4, ny = 36//4, nt = 12*110,    maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*110*4,  maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72,    ny = 36,    nt = 12*110*16, maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*110*64, maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)

# c5.5
# rf0
run_conf_write(nx = 72,   ny = 36,   nt = 12*20,     maxreflvls = 0, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72*2, ny = 36*2, nt = 12*20*4,   maxreflvls = 0, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72*4, ny = 36*4, nt = 12*20*16,  maxreflvls = 0, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72*8, ny = 36*8, nt = 12*20*64,  maxreflvls = 0, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
# rf1
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*20,     maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72,    ny = 36,    nt = 12*20*4,   maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*20*16,  maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72*4,  ny = 36*4,  nt = 12*20*64,  maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
# rf2
run_conf_write(nx = 72//4, ny = 36//4, nt = 12*20,     maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*20*4,   maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72,    ny = 36,    nt = 12*20*16,  maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*20*64,  maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)

os.chdir("interp")
exec(open("run_script.py").read())
os.chdir("../uniform/")
exec(open("run_script.py").read())
os.chdir("../intermediate/")
exec(open("run_script.py").read())
os.chdir("..")