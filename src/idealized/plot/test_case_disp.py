import numpy as np
import Ngl, Nio
# function for streamfunction
def fpsi(x, y, T, beta, ae):
    u0 = 2.0*np.pi*ae/T
    return -u0*(np.sin(y)*np.cos(beta) - np.cos(x)*np.cos(y)*np.sin(beta))

def solid(x, y, dx, dy, T, beta, ae):
    ypsi = np.zeros([ny + 1, nx])
    ypsi[:-1, :] = y[:, :] - dy/2.
    ypsi[-1, :]  = y[-1, :] + dy/2.
    xpsi = np.zeros([ny+1, nx])
    xpsi[:-1, :] = x[:, :]
    xpsi[-1, :]  = x[0, :]
    psi = fpsi(xpsi, ypsi, T, beta, ae)
    u = -(psi[1:, :] - psi[:-1, :])/dy
    ypsi = np.zeros([ny, nx + 1])
    ypsi[:, :-1] = y[:, :]
    ypsi[:, -1]  = y[:, -1]
    xpsi = np.zeros([ny, nx + 1])
    xpsi[:, 1:] = x[:, :] + dx/2.
    xpsi[:, 0]  = 0.
    cosp = np.cos(y)
    psi = fpsi(xpsi, ypsi, T, beta, ae)
    v = (psi[:, 1:] - psi[:, :-1])/dx/cosp
    return u, v

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

def plot(filename, outputname):
    # set up parameters
    f   = Nio.open_file(filename)
    q  = f.variables['q'][:, 0, :]
    qY  = f.variables['qY'][:, 0, :]
    ae  = 6.371229e6
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
    for it in range(1,2):
        cnres = Ngl.Resources()
        cnres.nglFrame               = False
        cnres.nglDraw                = False
        cnres.mpGridAndLimbOn        = False
        cnres.mpOutlineOn           = False
        cnres.mpPerimOn              = False#             ; turn off box around plot
        cnres.nglMaximize            = True
        cnres.mpProjection           = "Orthographic" # Change the map projection.
        x_c, y_c = qmoving(it*24.*3600., 12.*24.*3600., 0.*np.pi, ae)
        cnres.mpCenterLonF           = x_c*r2d
        cnres.mpCenterLatF           = y_c*r2d
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
        cnres.pmTickMarkDisplayMode  = "Never"
        contour1 = Ngl.contour_map(wks,q[it, 0:ncols[it]],cnres)
        cnres.cnLineColor            = "red"
        contour2 = Ngl.contour_map(wks,qY[it, 0:ncols[it]],cnres)
        Ngl.draw(contour1)
        Ngl.draw(contour2)
        # Ngl.draw(vc)
        Ngl.frame(wks)

        Ngl.end()

plot("../solid0.15intermediate/plots_solid_144_rf1_220_0.15/AMRDUST.nc", "solid0.15inter")
plot("../solid0.intermediate/plots_solid_144_rf1_52_0./AMRDUST.nc", "solid0.inter")