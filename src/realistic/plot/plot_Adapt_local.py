import Ngl, Nio
import numpy as np


p0mb = 101325.
p0mb1 = 1013.25
interp = 1
extrap = False
pnew = [800., 750.]
lev = 0
r2d = 57.2957795
f   = Nio.open_file("AMRDUST_Adapt.nc")
ncolsAdapt = f.variables['nCols'][:]
DU_AI   = f.variables['AI'][:, :, 1:3]
print(np.shape(DU_AI))
latAdapt = f.variables['lat'][:, :]*r2d
lonAdapt = f.variables['lon'][:, :]*r2d
vertx  = f.variables['vertx'][:, :, :]*r2d
verty  = f.variables['verty'][:, :, :]*r2d
NewTracer_AI_Adapt = np.ma.masked_where(DU_AI == 1.e30, DU_AI)

# day 10 to 20
# set up colormap

rlist    = Ngl.Resources()
rlist.wkColorMap = "WhiteYellowOrangeRed"
for it in range(15):
    print(it)
    # use colormap and type information to setup output background
    wks_type = "pdf"
    wks = Ngl.open_wks(wks_type,"Adapt"+str(it), rlist)
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
    #Level selection
    # res.lbBoxLinesOn  = False
    res.cnLevelSelectionMode   = "ExplicitLevels"   # Set explicit contour levels
    res.cnLevels               = np.arange(1e-5,8.e-4, 4e-5)
  
    # maximize the plot
    res.nglMaximize = True
    res.nglFrame = False
    # coordinate settings
    lontmp = lonAdapt[it, 0:ncolsAdapt[it]]
    lattmp = latAdapt[it, 0:ncolsAdapt[it]]
    lontmp  = np.where(lontmp > 180., lontmp - 360., lontmp)
    lon_list = np.where((lontmp >= -30) & (lontmp <=90))[0].tolist()
    # lon_list += np.where((lontmp >= 0) & (lontmp <=90))[0].tolist()
    lat_list = np.where((lattmp >= -10) & (lattmp <=50))[0].tolist()
    
    latlon_list = list(set(lat_list) & set(lon_list))
    res.sfXArray = lontmp[latlon_list]
    res.sfYArray = lattmp[latlon_list]
    # print(lon_list)
    # print(lat_list)
    contour = Ngl.contour_map(wks,NewTracer_AI_Adapt[it,latlon_list, lev], res)

    gsres                  = Ngl.Resources()
    gsres.gsLineColor      = "Gray25"
    gsres.gsLineThicknessF = 2.0
    for nc in latlon_list:
        Ngl.polyline(wks,contour,vertx[it, nc, :],verty[it, nc, :],gsres)
    Ngl.frame(wks)
    # Ngl.Draw(wks)
    Ngl.destroy(wks)
Ngl.end()
