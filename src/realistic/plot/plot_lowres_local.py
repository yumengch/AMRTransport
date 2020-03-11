import Ngl, Nio
import numpy as np


p0mb = 101325.
p0mb1 = 1013.25
interp = 1
extrap = False
pnew = [800., 750.]
lev = 0
# open the file
f   = Nio.open_file("LowRes/LowRes_200610.01_tracer.nc")
# interpolation from sigma coordinate to pressure coordinate
# get parameters and variables for interpolation
hyam = f.variables["hyam"][:]/p0mb
hybm = f.variables["hybm"][:]
DU_CI = f.variables["DU_CI"][:, :, :, :]*1e6
DU_AI = f.variables["DU_AI"][:, :, :, :]*1e6
PS    = f.variables["aps"][:, :, :]
lon = f.variables["lon"]
lon_low   = f.variables["lon"][:]

lat_low   = f.variables["lat"][:]

lon_list = np.where((lon_low >= 330) & (lon_low <=360))[0].tolist()
lon_list += np.where((lon_low >= 0) & (lon_low <=90))[0].tolist()
lat_list = np.where((lat_low >= -10) & (lat_low <=50))[0].tolist()
lon_low  = np.where(lon_low > 180., lon_low - 360., lon_low)
# start the interpolation
NewTracer_CI_Low = Ngl.vinth2p(DU_CI[:,:,:,:],hyam,hybm,pnew,PS[:,:,:],interp,p0mb1, 1,extrap)
NewTracer_AI_Low = Ngl.vinth2p(DU_AI[:,:,:,:],hyam,hybm,pnew,PS[:,:,:],interp,p0mb1, 1,extrap)
NewTracer_AI_Low = np.ma.masked_where(NewTracer_AI_Low == 1.e30, NewTracer_AI_Low)

# day 10 to 20
# set up colormap
rlist    = Ngl.Resources()
rlist.wkColorMap = "WhiteYellowOrangeRed"
for it in range(30):
    print(it)
    # use colormap and type information to setup output background
    wks_type = "pdf"
    wks = Ngl.open_wks(wks_type,"LowRes"+str(it), rlist)
    # resources for the contour
    res = Ngl.Resources()
    res.lbLabelBarOn          = False
    # Filled contour/labels/labelinfos
    res.cnLinesOn             = False
    res.cnFillOn              = True
    res.cnLineLabelsOn        = False
    res.cnInfoLabelOn         = False
    res.mpGridAndLimbOn       = False
    res.mpLimitMode       = "LatLon"              #-- must be set using minLatF/maxLatF/minLonF/maxLonF
    res.mpMinLatF         =  -10.                #-- sub-region minimum latitude
    res.mpMaxLatF         =  50.               #-- sub-region maximum latitude
    res.mpMinLonF         = -30.                 #-- sub-region minimum longitude
    res.mpMaxLonF         =  90.                  #-- sub-region maximum longitude    
    #anotherway to define colormap (overlay the predefined color map)
    # cmap = Ngl.read_colormap_file("WhiteBlueGreenYellowRed")
    #Level selection
    # res.lbBoxLinesOn  = False
    res.cnLevelSelectionMode   = "ExplicitLevels"   # Set explicit contour levels
    res.cnLevels               = np.arange(1e-5,8.e-4, 4e-5)      # 0,5,10,...,70
    # maximize the plot
    res.nglMaximize = True
    res.nglFrame = False
    # coordinate settings
    sliceTracer = NewTracer_AI_Low[it, lev, lat_list, :]
    res.sfXArray = lon_low[lon_list]
    res.sfYArray = lat_low[lat_list]
    res.vpWidthF = 1
    res.vpHeightF = 0.5
    # Ngl.contour_map(wks,NewTracerPlot,res)
    Ngl.contour_map(wks,sliceTracer[:, lon_list],res)
    Ngl.frame(wks)
    # Ngl.Draw(wks)
    Ngl.destroy(wks)
Ngl.end()
