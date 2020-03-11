import Ngl, Nio
import numpy as np
p0mb = 101325.
p0mb1 = 1013.25
interp = 1
extrap = False
pnew = [800., 750.]
lev = 0
# open the file
f   = Nio.open_file("/home/yumeng/realthesis/HighRes/HighRes_200610.01_tracer.nc")
hyam = f.variables["hyam"][:]/p0mb
hybm = f.variables["hybm"][:]
DU_CI = f.variables["DU_CI"][:, :, :, :]*1e6
DU_AI = f.variables["DU_AI"][:, :, :, :]*1e6
PS    = f.variables["aps"][:, :, :]
lon_high   = f.variables["lon"][:]
lat_high   = f.variables["lat"][:]
# # interpolation from sigma coordinate to pressure coordinate
NewTracer_CI_High = Ngl.vinth2p(DU_CI[:,:,:,:],hyam,hybm,pnew,PS[:,:,:],interp,p0mb1, 1,extrap)
NewTracer_AI_High = Ngl.vinth2p(DU_AI[:,:,:,:],hyam,hybm,pnew,PS[:,:,:],interp,p0mb1, 1,extrap)
NewTracer_AI_High = np.ma.masked_where(NewTracer_AI_High == 1.e30, NewTracer_AI_High)
print(np.max(NewTracer_AI_High[:, 0, :, :]))
# day 10 to 20
# set up colormap
rlist    = Ngl.Resources()
rlist.wkColorMap = "WhiteYellowOrangeRed"
for it in range(30):
    print(it)
    # use colormap and type information to setup output background
    wks_type = "pdf"
    wks = Ngl.open_wks(wks_type,"HighRes"+str(it), rlist)
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
    res.cnLevelSelectionMode   = "ExplicitLevels"   # Set explicit contour levels
    res.cnLevels               = np.arange(1e-5,8.e-4, 4e-5)      # 0,5,10,...,70
    # maximize the plot
    contour = []
    res.nglMaximize = True
    res.nglFrame = False
    # coordinate settings
    NewTracerPlot, LonPlot = Ngl.add_cyclic(NewTracer_AI_High[it, lev, :, :], lon_high)
    res.sfXArray = LonPlot
    res.sfYArray = lat_high[:]
    res.vpWidthF = 1
    res.vpHeightF = 0.5
    Ngl.contour_map(wks,NewTracerPlot,res)

    Ngl.frame(wks)
    Ngl.destroy(wks)
Ngl.end()
