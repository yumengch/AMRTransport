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
    
    base_dir = "../"+initial[1:-1]+beta+"interp"
    if exp is not None:
        base_dir += exp
    if not os.path.isdir(base_dir):
        subprocess.call(["mkdir", base_dir])
    base_dir += "/rf"+str(maxreflvls)
    if not os.path.isdir(base_dir):
        subprocess.call(["mkdir", base_dir])    
    name_dir = base_dir+"plots_" + initial[1:-1] + "_" +str(nx) + "_rf" + str(maxreflvls) \
        + "_"+ str(nt) + "_" +beta
    name_log = "nx"+str(nx)+"_rf"+str(maxreflvls) + ".log"
    name_vtu = initial[1:-1]+"_"
    results(name_dir, name_log, name_vtu)

subprocess.call(["bash", "run.sh"])

# Moving Vortices (uniform refinement + interpolated wind)
# beta = 0.25pi
# c1
# rf1
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*110,    maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")
run_conf_write(nx = 72,    ny = 36,    nt = 12*110*4,  maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*110*16, maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")
#run_conf_write(nx = 72*4,  ny = 36*4,  nt = 12*110*64, maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")
# rf2
run_conf_write(nx = 72//4, ny = 36//4, nt = 12*110,    maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*110*4,  maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")
run_conf_write(nx = 72,    ny = 36,    nt = 12*110*16, maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")
#run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*110*64, maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")
# c5.5
# rf1
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*20,    maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")
run_conf_write(nx = 72,    ny = 36,    nt = 12*20*4,  maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*20*16, maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")
#run_conf_write(nx = 72*4,  ny = 36*4,  nt = 12*20*64, maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")
# rf2
run_conf_write(nx = 72//4, ny = 36//4, nt = 12*20,    maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*20*4,  maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")
run_conf_write(nx = 72,    ny = 36,    nt = 12*20*16, maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")
#run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*20*64, maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = -1, theta_c = -1, exp = "uniform")

# adaptive:(adaptive refinement + interpolated wind)
# c1
# rf1
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*110,    maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72,    ny = 36,    nt = 12*110*4,  maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*110*16, maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
#run_conf_write(nx = 72*4,  ny = 36*4,  nt = 12*110*64, maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
# rf2
run_conf_write(nx = 72//4, ny = 36//4, nt = 12*110,    maxreflvls = 2, beta = "0.25",  initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*110*4,  maxreflvls = 2, beta = "0.25",  initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72,    ny = 36,    nt = 12*110*16, maxreflvls = 2, beta = "0.25",  initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
#run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*110*64, maxreflvls = 2, beta = "0.25",  initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
# c5.5
# rf1
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*20,    maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72,    ny = 36,    nt = 12*20*4,  maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*20*16, maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
#run_conf_write(nx = 72*4,  ny = 36*4,  nt = 12*20*64, maxreflvls = 1, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
# rf2
run_conf_write(nx = 72//4, ny = 36//4, nt = 12*20,    maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72//2, ny = 36//2, nt = 12*20*4,  maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
run_conf_write(nx = 72,    ny = 36,    nt = 12*20*16, maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)
#run_conf_write(nx = 72*2,  ny = 36*2,  nt = 12*20*64, maxreflvls = 2, beta = "0.25", initial = "\"moving\"", nof = "12", theta_r = 0.8, theta_c = 0.4)