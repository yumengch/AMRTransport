import Ngl, Nio
import numpy as np


p0mb = 101325.
p0mb1 = 1013.25
interp = 1
extrap = False
pnew = [800., 750.]
lev = 0
r2d = 57.2957795
f   = Nio.open_file("/home/yumeng/WriteMesh Real/AMRDUST.nc")
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
for it in range(25):
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
    #Level selection
    # res.lbBoxLinesOn  = False
    res.cnLevelSelectionMode   = "ExplicitLevels"   # Set explicit contour levels
    res.cnLevels               = np.arange(1e-5,8.e-4, 4e-5)
  
    # maximize the plot
    res.nglMaximize = True
    res.nglFrame = False
    # coordinate settings
    res.sfXArray = lonAdapt[it, 0:ncolsAdapt[it]]
    res.sfYArray = latAdapt[it, 0:ncolsAdapt[it]]
    contour = Ngl.contour_map(wks,NewTracer_AI_Adapt[it,0:ncolsAdapt[it], lev], res)

    gsres                  = Ngl.Resources()
    gsres.gsLineColor      = "Gray25"
    gsres.gsLineThicknessF = 2.0
    for nc in range(ncolsAdapt[it]):
        Ngl.polyline(wks,contour,vertx[it, nc, :],verty[it, nc, :],gsres)
    Ngl.frame(wks)
    # Ngl.Draw(wks)
    Ngl.destroy(wks)
Ngl.end()
