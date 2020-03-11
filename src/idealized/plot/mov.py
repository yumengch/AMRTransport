import numpy as np
import Ngl, Nio
def qmoving(t, T, beta, ae):
    u0 = 2*np.pi*ae/T
    x0 = 1.5*np.pi; y0 = 0.

    x_r = np.arctan2(np.cos(y0)*np.sin(x0 - np.pi), np.cos(y0)*np.sin(0.5*np.pi - beta)*np.cos(x0 - np.pi) - np.cos(0.5*np.pi - beta)*np.sin(y0))
    if (x_r < 0.):
        x_r = x_r + 2.*np.pi
    y_r = np.arcsin(np.sin(y0)*np.sin(0.5*np.pi - beta) + np.cos(y0)*np.cos(0.5*np.pi - beta)*np.cos(x0 - np.pi))

    x_c = np.arctan2(np.cos(y_r)*np.sin(x_r + u0*t/ae), np.sin(y_r)*np.cos(0.5*np.pi - beta) + np.cos(y_r)*np.cos(x_r + u0*t/ae)*np.sin(0.5*np.pi - beta)) + np.pi
    y_c = np.arcsin(np.sin(y_r)*np.sin(0.5*np.pi - beta) - np.cos(y_r)*np.cos(0.5*np.pi - beta)*np.cos(x_r + u0*t/ae))

    return x_c, y_c

def plot(filename, outputend):
    # set up parameters
    f   = Nio.open_file(filename)
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

    # Create contour and map resource lists.
    for it in range(13):
        wks = Ngl.open_wks(wks_type,"Moving_"+str(it)+outputend, rlist)
        cnres = Ngl.Resources()
        cnres.nglFrame               = False
        # cnres.nglDraw                = False
        cnres.mpGridAndLimbOn        = False
        cnres.mpOutlineOn            = False
        cnres.mpPerimOn              = False#             ; turn off box around plot
        cnres.nglMaximize            = True
        cnres.mpProjection           = "Orthographic" # Change the map projection.
        # x_c, y_c = 0, 90.#qmoving(it*24.*3600., 12.*24.*3600., 0.25*np.pi, ae)
        cnres.mpCenterLonF           = x_c*r2d           # Rotate the projection.
        cnres.mpCenterLatF           = y_c*r2d           # Rotate the projection.
        cnres.mpOutlineOn            = False
        cnres.sfXArray = lon[it, 0:ncols[it]]
        cnres.sfYArray = lat[it, 0:ncols[it]]
        cnres.cnLinesOn              = False
        cnres.cnFillOn               = True
        cnres.cnLineLabelsOn         = False
        cnres.cnInfoLabelOn          = False
        cnres.lbLabelBarOn           = True
        cnres.cnLevelSelectionMode = "ManualLevels"
        cnres.cnMinLevelValF       = 0.4
        cnres.cnMaxLevelValF       = 1.6
        cnres.cnLevelSpacingF      = 0.05

        cnres.pmTickMarkDisplayMode  = "Never"
        contour1 = Ngl.contour_map(wks,q[it, 0:ncols[it]],cnres)
        gsres                  = Ngl.Resources()
        gsres.gsLineColor      = "Gray25"
        gsres.gsLineThicknessF = 3.0
        # print(vertx[0, :, 0],verty[0, :, 0])    
        for nc in range(ncols[it]):
            Ngl.polyline(wks,contour1,vertx[it, :, nc],verty[it, :, nc],gsres)    
        # Ngl.draw(contour1)
        Ngl.frame(wks)
        Ngl.destroy(wks)
    Ngl.end()
plot("../moving0.25interp/rf1/plots_moving_72_rf1_960_0.25/AMRDUST.nc", "interpAMR"):
plot("../moving0.25/rf0/plots_moving_72_rf0_240_0.25/AMRDUST.nc", "noAMR"):