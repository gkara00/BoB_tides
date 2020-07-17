import xarray as xr
from scipy.io import netcdf
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import cmocean

import datetime 
import os

from IPython.display import set_matplotlib_formats
set_matplotlib_formats('png', 'pdf')

plt.rcParams["font.size"] = 12


def read_netcdf_xr(fid):
    ds = xr.open_mfdataset(fid)
    return ds


if __name__ == "__main__": 
    root_dir = '/home/gkara/Documents/for-tide_paper/'

    nc = netcdf.netcdf_file(root_dir + 'areatabmask.nc','a')
    lon = nc.variables['nav_lon'][:]
    lat = nc.variables['nav_lat'][:]
    area = nc.variables['area'][:]
    nc.close()
    nc = netcdf.netcdf_file(root_dir + 'bathy_meter.nc','r')
    bathy = nc.variables['Bathymetry'][:]
    nc.close()

    curl_ref = read_netcdf_xr(root_dir + 'CURL_ref/BISCAY-T60_1h_CURL*')
    curl_ref = curl_ref.set_coords(['nav_lon', 'nav_lat'])
    #
    curl_exp1 = read_netcdf_xr(root_dir + 'CURL_exp1/BISCAY-T60_1h_CURL*')
    curl_exp1 = curl_exp1.set_coords(['nav_lon', 'nav_lat'])

    vmin = -.5e-04
    vmax = .5e-04
    levels = np.linspace(vmin,vmax,41, endpoint=True)

    for n in range(380,720):
        print('***' + str(n) + '***')
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(13.5,6))
        fig.subplots_adjust(wspace=0.08)

        curl_ref.socurl[n,::].plot.contourf(ax=axes[0], x='nav_lon', y='nav_lat', levels=levels, vmin=vmin, vmax=vmax, cmap=cmocean.cm.balance, add_colorbar=False)
        curl_exp1.socurl[n,::].plot.contourf(ax=axes[1], x='nav_lon', y='nav_lat', levels=levels, vmin=vmin, vmax=vmax, cmap=cmocean.cm.balance, add_colorbar=False)
        mask = np.isnan(curl_ref.socurl[n,:,:])
        mask.plot.contour(ax=axes[0], x='nav_lon', y='nav_lat', linewidths=0.2,colors='black')
        mask.plot.contour(ax=axes[1], x='nav_lon', y='nav_lat', linewidths=0.2,colors='black')

        axes[1].yaxis.tick_right()
        axes[0].yaxis.set_ticks_position('both')
        axes[1].yaxis.set_ticks_position('both')
        axes[1].set_ylabel(None)

        p0 = axes[0].get_position().get_points().flatten()
        p1 = axes[1].get_position().get_points().flatten()
        ax_cbar = fig.add_axes([p0[0]+.04, 0, p1[2]-p0[0]-.08, 0.025])
        cb = plt.colorbar(im0, cax=ax_cbar, orientation='horizontal', extend='both', label='Relative_Vorticity $\zeta$')
        cb.set_ticks(np.arange(-.5e-04,1e-03,.2e-04))
        cb.formatter.set_powerlimits((0,0))
        cb.update_ticks()

        #save figure
        fig.savefig(root_dir + '_plots/curl%04d.png'%n,dpi=300,bbox_inches='tight')
        plt.close()