from netCDF4 import Dataset
from scipy.io import netcdf #### <--- This is the library to import.
import xarray as xr
import xesmf as xe
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import cmocean.cm as cmo
import cmaps
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import cmocean
import multiprocessing

import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature

import datetime 
import os

%matplotlib inline
 
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('png', 'pdf')

plt.rcParams["font.size"] = 12


def read_netcdf(fid, var):
    from netCDF4 import Dataset
    nc = Dataset(fid,'r')
    out = nc[var][:]
    return out

def read_netcdf_xr(fid):
    ds = xr.open_mfdataset(fid, combine='by_coords')
#     out = ds[var]
    return ds

def plot_map(ax,lon,lat,c,vmin,vmax,color, ind, title=''):
    x1 = -10
    x2 = 0
    y1 = 42
    y2 = 51
    
#     mm = ax.pcolormesh(lon,lat,c, vmin=vmin,vmax=vmax,transform=ccrs.PlateCarree(),cmap=color)
    v = np.linspace(vmin,vmax,40, endpoint=True)
    mm = ax.contourf(lon,lat,c, v , vmin=vmin,vmax=vmax, extend='both', transform=ccrs.PlateCarree(),cmap=color)

    ax.contour(lon,lat,bathy,levels=[200], colors='k', linestyles='-', linewidth=1)
    
    ax.set_xticks([0,-2.5,-5,-7.5,-10], crs=ccrs.PlateCarree())
    ax.set_yticks([42.5, 45, 47.5, 50], crs=ccrs.PlateCarree())
    ax.set_title(title)
#
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
#
#     ax.xaxis.set_major_formatter(lon_formatter)
#     ax.yaxis.set_major_formatter(lat_formatter)

    land_10m = cartopy.feature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])
    ax.add_feature(land_10m)
    ax.coastlines(resolution='10m',)
    ax.set_extent([x1,x2,y1,y2])
    
    ax.text(-1.5, 42.4, ind,weight = 'bold', fontsize=12,
     horizontalalignment='left',
     transform=ccrs.PlateCarree())
#     cb = plt.colorbar(im,ax=ax,extend='both', pad=0.02)

#     plt.colorbar(mm,ax=ax,orientation='horizontal',extend='max',pad=0.1)
    
    return mm


def map_generator(ax, lon, lat, c, vmin, vmax, color, ind, title=''):
    x1 = -10
    x2 = 0
    y1 = 42
    y2 = 51

#     mm = ax.pcolormesh(lon,lat,c, vmin=vmin,vmax=vmax,transform=ccrs.PlateCarree(),cmap=color)
    v = np.linspace(vmin,vmax,40, endpoint=True)
    mm = ax.contourf(lon,lat,c, v , vmin=vmin,vmax=vmax, extend='both', transform=ccrs.PlateCarree(),cmap=color)

    ax.contour(lon,lat,bathy,levels=[200], colors='k', linestyles='-', linewidth=1)

    ax.set_xticks([0,-2.5,-5,-7.5,-10], crs=ccrs.PlateCarree())
    ax.set_yticks([42.5, 45, 47.5, 50], crs=ccrs.PlateCarree())
    ax.set_title(title)
#
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
#
#     ax.xaxis.set_major_formatter(lon_formatter)
#     ax.yaxis.set_major_formatter(lat_formatter)

    land_10m = cartopy.feature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])
    ax.add_feature(land_10m)
    ax.coastlines(resolution='10m',)
    ax.set_extent([x1,x2,y1,y2])

    ax.text(-6.5, 42.4, ind,weight = 'bold', fontsize=12,
     horizontalalignment='left',
     transform=ccrs.PlateCarree())
#     cb = plt.colorbar(im,ax=ax,extend='both', pad=0.02)

#     plt.colorbar(mm,ax=ax,orientation='horizontal',extend='max',pad=0.1)

    return mm

def plot_generator(k, ds1, ds2):
    '''fname: filename
       ds1: first xarray dataset
       ds2: second xarray dataset'''
    
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(13.5,6),subplot_kw={'projection': ccrs.PlateCarree()})
    fig.subplots_adjust(wspace=0.08)

    for n in range(24):

        im0 = map_generator(axes[0],ds1.nav_lon.data,ds1.nav_lat.data,ds1.socurl.data[n,:,:],vmin=-.6e-04,vmax=.6e-04, color=cmocean.cm.balance, ind=str(ds1.time_counter[n].data)[:19], title='surf. Relative Vorticity $\zeta$ -- With Tides')
        im1 = map_generator(axes[1],ds2.nav_lon.data, ds2.nav_lat.data,ds2.socurl.data[n,:,:],vmin=-.6e-04,vmax=.6e-04, color=cmocean.cm.balance, ind=str(ds2.time_counter[n].data)[:19], title='surf. Relative Vorticity $\zeta$ -- Without Tides')

        axes[1].yaxis.tick_right()
        axes[0].yaxis.set_ticks_position('both')
        axes[1].yaxis.set_ticks_position('both')

        p0 = axes[0].get_position().get_points().flatten()
        p1 = axes[1].get_position().get_points().flatten()
        ax_cbar = fig.add_axes([p0[0]+.04, .03, p1[2]-p0[0]-.08, 0.025])
        cb = plt.colorbar(im0, cax=ax_cbar, orientation='horizontal', extend='both')
        cb.set_ticks(np.arange(-.6e-04,1e-03,.2e-04))
        cb.formatter.set_powerlimits((0,0))
        cb.update_ticks()
        fig.savefig(root_dir + 'curl%04d.png'%(24*k+n),dpi=300,bbox_inches='tight')
        plt.clf()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
root_dir = '/home/gkara/Documents/for-tide_paper/'
_f = ['BISCAY-T60_1h_CURL_20110928-20110928.nc', 'BISCAY-T60_1h_CURL_20120224-20120224.nc']
#
curl_ref = read_netcdf_xr([root_dir + 'CURL_ref/' + f for f in _f])
curl_ref = curl_ref.set_coords(['nav_lon', 'nav_lat'])
#
curl_exp1 = read_netcdf_xr([root_dir + 'CURL_exp1/' + f for f in _f])
curl_exp1 = curl_exp1.set_coords(['nav_lon', 'nav_lat'])
#
bathy = read_netcdf_xr(root_dir + 'bathy_meter.nc')
bathy = bathy.set_coords(['nav_lon', 'nav_lat'])




#
def plot_200_mask(axes, mask, bathy=bathy):
    for ax in axes:
        mask.plot.contour(ax=ax, x='nav_lon', y='nav_lat', linewidths=0.1, colors='k')
        bathy.Bathymetry.plot.contour(ax=ax, x='nav_lon', y='nav_lat', levels=[200], linewidths=1.1, colors='k')
    return
    
#
vmin = -.6e-04
vmax = .6e-04
levels = np.linspace(vmin,vmax,49, endpoint=True)

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12,10))
fig.subplots_adjust(wspace=0.08, hspace=.14)
#
mask = np.isnan(curl_ref.socurl[0,:,:])
im0 = curl_ref.sel(time_counter='2011-09-28T15:30:00.000000000').socurl.plot.contourf(ax=axes[0,0], x='nav_lon', y='nav_lat', levels=levels, vmin=vmin, vmax=vmax, cmap=cmocean.cm.balance, add_colorbar=False)
im1 = curl_exp1.sel(time_counter='2011-09-28T15:30:00.000000000').socurl.plot.contourf(ax=axes[0,1], x='nav_lon', y='nav_lat', levels=levels, vmin=vmin, vmax=vmax, cmap=cmocean.cm.balance, add_colorbar=False)
plot_200_mask(axes=[axes[0,0], axes[0,1]], mask=mask)
#
p0 = axes[0,0].get_position().get_points().flatten()
p1 = axes[0,1].get_position().get_points().flatten()
ax_cbar = fig.add_axes([p0[0]+.04, 0.035, p1[2]-p0[0]-.08, 0.025])
cb = plt.colorbar(im0, cax=ax_cbar, orientation='horizontal', extend='both', label='Relative_Vorticity $\zeta$')
cb.set_ticks(np.arange(-.6e-04,1e-03,.1e-04))
cb.formatter.set_powerlimits((0,0))
cb.update_ticks()
#
im2 = curl_ref.sel(time_counter='2012-02-24T15:30:00.000000000').socurl.plot.contourf(ax=axes[1,0], x='nav_lon', y='nav_lat', levels=levels, vmin=vmin, vmax=vmax, cmap=cmocean.cm.balance, add_colorbar=False)
im3 = curl_exp1.sel(time_counter='2012-02-24T15:30:00.000000000').socurl.plot.contourf(ax=axes[1,1], x='nav_lon', y='nav_lat', levels=levels, vmin=vmin, vmax=vmax, cmap=cmocean.cm.balance, add_colorbar=False)
plot_200_mask(axes=[axes[1,0], axes[1,1]], mask=mask)
#
axes[0,0].set_xlabel(None); axes[0,1].set_xlabel(None)
axes[0,0].set_xticklabels([]); axes[0,1].set_xticklabels([])
axes[0,1].set_ylabel(None); axes[1,1].set_ylabel(None)
axes[0,1].set_yticklabels([]); axes[1,1].set_yticklabels([])
#
axes[0,0].set_title("(a) TON: 2011-09-28"); axes[0,1].set_title("(b) TOFF: 2011-09-28")
axes[1,0].set_title("(c) TON: 2012-02-24"); axes[1,1].set_title("(d) TOFF: 2012-02-24")

#save figure
fig.savefig(root_dir + '_plots/fig_vorticity.png',dpi=300,bbox_inches='tight')


#ANIMAITONS
curl_ref = read_netcdf_xr(root_dir + 'CURL_ref/BISCAY-T60_1h_CURL*')
curl_ref = curl_ref.set_coords(['nav_lon', 'nav_lat'])
#
curl_exp1 = read_netcdf_xr(root_dir + 'CURL_exp1/BISCAY-T60_1h_CURL*')
curl_exp1 = curl_exp1.set_coords(['nav_lon', 'nav_lat'])

vmin = -.5e-04
vmax = .5e-04
levels = np.linspace(vmin,vmax,41, endpoint=True)

for n in range(1):
    print('***' + str(n) + '***')
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(13.5,6))
    fig.subplots_adjust(wspace=0.08)

    im0 = curl_ref.socurl[n,::].plot.contourf(ax=axes[0], x='nav_lon', y='nav_lat', levels=levels, vmin=vmin, vmax=vmax, cmap=cmocean.cm.balance, add_colorbar=False)
    im1 = curl_exp1.socurl[n,::].plot.contourf(ax=axes[1], x='nav_lon', y='nav_lat', levels=levels, vmin=vmin, vmax=vmax, cmap=cmocean.cm.balance, add_colorbar=False)
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
    cb.set_ticks(np.arange(-.5e-04,1e-03,.1e-04))
    cb.formatter.set_powerlimits((0,0))
    cb.update_ticks()
    
    #save figure
#     fig.savefig(root_dir + '_plots/curl%04d.png'%n,dpi=300,bbox_inches='tight')
#     plt.close()