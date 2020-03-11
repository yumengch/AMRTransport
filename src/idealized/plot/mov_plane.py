import Ngl, Nio
import numpy as np
# open the file
f   = Nio.open_file("../moving0.25interp/rf1/plots_moving_72_rf1_960_0.25/AMRDUST.nc")
q  = f.variables['q'][:, 0, :]
r2d = 57.2957795
lat = f.variables['lat'][:, :]*r2d
lon = f.variables['lon'][:, :]*r2d
vertx = f.variables['vertx'][:, :, :]*r2d
verty = f.variables['verty'][:, :, :]*r2d
ncols = f.variables['nCols'][:]
ae  = 6.371229e6
#  Open a workstation.
wks_type = "pdf"
rlist    = Ngl.Resources()
rlist.wkColorMap = "WhBlReWh"
# print(l)
# l0 = np.ma.masked_values(l, np.min(l))
# l1 = np.ma.masked_values(l0, np.min(l) + 1)
# for i in range(len(lon)):
#     for j in range(len(lat)):
#         if not np.ma.is_masked(l1[0, j, i]):
#             l1[0, j, i] = j
        # print(, l1[0, j, i])
# print(np.min(l), np.max(l))
# # set up colormap
# rlist    = Ngl.Resources()
# rlist.wkColorMap = "wh-bl-gr-ye-re"

for it in range(13):
    # use colormap and type information to setup output background
    wks_type = "pdf"
    wks = Ngl.open_wks(wks_type,"mov_plane"+str(it), rlist)
    # resources for the contour
    res = Ngl.Resources()
    res.lbLabelBarOn = False
    # Filled contour/labels/labelinfos
    res.cnFillOn              = True
    res.cnLinesOn             = False
    res.cnLineLabelsOn        = False
    res.cnInfoLabelOn         = False
    #anotherway to define colormap (overlay the predefined color map)
    # cmap = Ngl.read_colormap_file("WhiteBlueGreenYellowRed")
    #Level selection
    res.cnLevelSelectionMode = "ManualLevels"
    res.cnMinLevelValF       = 0.4
    res.cnMaxLevelValF       = 1.6
    res.cnLevelSpacingF      = 0.05
    # maximize the plot
    res.nglMaximize = True
    res.nglFrame = False
    # res.nglDraw = False
    # map background
    # res.mpGridAndLimbOn        = False
    # res.mpGeophysicalLineColor = "Black"
    # coordinate settings
    res.sfXArray = lon[it, 0:ncols[it]]
    res.sfYArray = lat[it, 0:ncols[it]]
    res.vpWidthF = 1
    res.vpHeightF = 0.5
    # Set resources necessary to get map projection correct.
    # res.mpCenterLonF       = 180 

    # label settings
    lon_values  = Ngl.fspan(np.min(lon), np.max(lon), 5)
    lon_labels  = ["0","90","180","270","360"]
    lat_values  = Ngl.fspan(np.min(lat), np.max(lat), 5)
    lat_labels = ["-90", "-45", "0", "45", "90"]
    res.tmXBMode                 = "Explicit"    # Label X and Y axes
    res.tmYLMode                 = "Explicit"    # with explicit labels.
    res.tmXBValues               = lon_values
    res.tmYLValues               = lat_values
    res.tmXBLabels               = lon_labels
    res.tmYLLabels               = lat_labels

    contour = Ngl.contour(wks,q[it, 0:ncols[it]],res)

    gsres                  = Ngl.Resources()
    gsres.gsLineColor      = "Gray25"
    gsres.gsLineThicknessF = 3.0
    # print(vertx[0, :, 0],verty[0, :, 0])    
    for nc in range(ncols[it]):
        Ngl.polyline(wks,contour,vertx[it, :, nc],verty[it, :, nc],gsres)    
    Ngl.frame(wks)
Ngl.end()
