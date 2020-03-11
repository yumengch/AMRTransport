import numpy as np
import Ngl, Nio
def plot(filename, outputname, centerlon, centerlat):
    # set up parameters
    f   = Nio.open_file(filename)
    q  = f.variables['q'][:, 0, :]
    r2d = 57.2957795
    lat = f.variables['lat'][:, :]*r2d
    lon = f.variables['lon'][:, :]*r2d
    vertx = f.variables['vertx'][:, :, :]*r2d
    verty = f.variables['verty'][:, :, :]*r2d
    ncols = f.variables['nCols'][:]

    #  Open a workstation.
    wks_type = "pdf"
    wks = Ngl.open_wks(wks_type,outputname)
    # Create contour and map resource lists.
    for it in range(13):
        cnres = Ngl.Resources()
        cnres.nglFrame               = False
        cnres.nglDraw                = False
        cnres.mpGridAndLimbOn        = False
        cnres.mpOutlineOn           = False
        cnres.mpPerimOn              = False#             ; turn off box around plot
        cnres.nglMaximize            = True
        cnres.mpProjection           = "Orthographic" # Change the map projection.
        cnres.mpCenterLonF           = centerlon
        cnres.mpCenterLatF           = centerlat
        cnres.mpOutlineOn            = False
        cnres.sfXArray = lon[it, 0:ncols[it]]
        cnres.sfYArray = lat[it, 0:ncols[it]]
        cnres.cnLinesOn              = True
        cnres.cnFillOn               = False
        cnres.cnLineLabelsOn         = False
        cnres.cnInfoLabelOn          = False
        cnres.cnLevelSelectionMode   = "ExplicitLevels"   # Set explicit contour levels
        cnres.cnLevels               = np.arange(0.1, 1., 0.1)      # 0,5,10,...,70
        cnres.cnLineThicknessF       = 3.
        cnres.cnLineColor            = "red"
        cnres.pmTickMarkDisplayMode  = "Never"
        contour1 = Ngl.contour_map(wks,q[it, 0:ncols[it]],cnres)
        gsres                  = Ngl.Resources()
        gsres.gsLineColor      = "Gray25"
        gsres.gsLineThicknessF = 2.0    
        for nc in range(ncols[it]):
            Ngl.polyline(wks,contour1,vertx[it, :, nc],verty[it, :, nc],gsres)    
        Ngl.draw(contour1)
    Ngl.frame(wks)
    Ngl.end()
plot('../solid0./plots_solid_72_rf1_26_0./AMRDUST.nc', 'solid0illus270', 270., 0.)
plot('../solid0./plots_solid_72_rf1_26_0./AMRDUST.nc', 'solid0illus90', 90., 0.)
plot('../solid0.5/plots_solid_72_rf1_960_0.5/AMRDUST.nc', 'solid05illus_90', 0., -90.)
plot('../solid0.5/plots_solid_72_rf1_960_0.5/AMRDUST.nc', 'solid05illus90', 0., 90.)