import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import pandas as pd
from pandas import Timestamp
from wrf import to_np, getvar, smooth2d, get_basemap, latlon_coords, extract_times, ALL_TIMES, interplevel
import sys, os, glob

mpl.font_manager._rebuild()

mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.style'] = 'normal'

def landmarks():
    landmark_dict = {'C':(-64.2123,-31.3154),
                    'M':(-68.7987,-32.8278),
                    '3':(-64.1131,-32.1767),
                    '4':(-64.34992,-33.13067), 
                    'Y':(-64.7545,-32.1062),
                    'S':(-63.978435,-31.715689),
                    'SR':(-68.285522,-34.587997),
                    'SL':(-66.33694,-33.30278)}
    return landmark_dict

landmark_dict=landmarks()

#run='2015-11-07_00:00:00'
#tindex=np.all([])
#tindex=np.arange(0,24,1) 

run=sys.argv[1]
tindex=int(sys.argv[2])


############################tCONFIGURATION PARAMETERS###############################

titlestr='Servicio Meteorológico Nacional GFS-WRF '
modname='SMN_WRF'
filenames=sorted(glob.glob('/data/meso/a/jmulhol2/RELAMPAGO_Dry_Run_Data/7_8_Nov_2015/wrf_smn/wrfout_d01_*'))
outpath='/data/meso/a/snesbitt/wrf3911/fcsts/'+modname+'/'+run+'/images'

############################tCONFIGURATION PARAMETERS###############################

#print(filenames)
print(run,tindex)

os.system('mkdir -p '+outpath)

basetime=pd.to_datetime(run, format='%Y-%m-%d_%H:%M:%S')

files=[]
times=[]
#times=pd.date_range(start=basetime,periods=43,freq='H')

for file in filenames:
        #files.append(Dataset('/data/meso/a/jmulhol2/RELAMPAGO_Dry_Run_Data/12_Nov_2015/UIUC_WRF_runs/wrfout_d01_'+params['times'][params['time_index']].strftime('%Y-%m-%d %H:%M')))
        files.append(Dataset(file))
        print(os.path.basename(file)[11:])
        times.append(pd.to_datetime(os.path.basename(file)[11:],format='%Y-%m-%d_%H:%M:%S'))
# In[544]:

def make_plot(cffield,lfield,lfield2,ufld,vfld,params):
# Get the latitude and longitude points

    print(params['time_index'])
    
    lats, lons = latlon_coords(cffield)

    # Get the basemap object

    #Original large domain ---    
    bm = Basemap(projection='lcc',width=3000*550,height=3000*375,
             resolution='i',lat_1=-32.8,lat_2=-32.8,lat_0=-32.8,lon_0=-67.0)

    #Cordoba inset ---- 
    #bm = Basemap(projection='lcc',width=1500*550,height=1500*375,
    #          resolution='i',lat_1=-32.2,lat_2=-32.2,lat_0=-32.2,lon_0=-65.0) 
    
    #Mendoza inset ---- 
    #bm = Basemap(projection='lcc',width=1500*550,height=1500*375,
    #          resolution='i',lat_1=-33.2,lat_2=-33.2,lat_0=-33.2,lon_0=-69.0) 
    
    # Create a figure
    fig = plt.figure(figsize=(12,9))

    # Add geographic outlines
    bm.drawcoastlines(linewidth=0.75)
    bm.drawstates(linewidth=1.)
    bm.drawcountries(linewidth=1.)

    # Convert the lats and lons to x and y.  Make sure you convert the lats and lons to
    # numpy arrays via to_np, or basemap crashes with an undefined RuntimeError.
    x, y = bm(to_np(lons), to_np(lats))


    if lfield is not None:
        CS=bm.contour(x, y, to_np(lfield), 10, colors="black", levels=params['llevels'],linewidths=1.0)
        plt.clabel(CS, inline=1, fontsize=12, fmt='%d')

    if lfield2 is not None:
        CS=bm.contour(x, y, to_np(lfield2), 10, colors="dimgrey", levels=params['llevels2'],linewidths=2.25)
        #plt.clabel(CS, inline=1, fontsize=12, fmt='%d')
    
    if ufld is not None:
        bm.barbs(x[::params['skip'],::params['skip']], 
                 y[::params['skip'],::params['skip']], 
                 to_np(ufld[::params['skip'],::params['skip']]),
                 to_np(vfld[::params['skip'],::params['skip']]), length=5, linewidth=0.75, zorder=10)

    if not('lalpha' in params):
        params['lalpha']=None
        

    # Draw the contours and filled contours
    bm.contourf(x, y, to_np(cffield), 10, cmap=get_cmap(params['ccmap']), levels=params['clevels'], extend='both',
               alpha=params['lalpha'])


    parallels = np.arange(-50.,-10.,2.)
    # labels = [left,right,top,bottom]
    bm.drawparallels(parallels,labels=[False,True,False,False],linewidth=0.5,dashes=[2,2])
    meridians = np.arange(-90.,-50.,2.)
    bm.drawmeridians(meridians,labels=[False,False,False,True],linewidth=0.5,dashes=[2,2])

    # Add a color bar
    plt.colorbar(shrink=.62, extend='both')

    timediff=params['times'][params['time_index']]-params['times'][0]
    timediff_secs=int(timediff.total_seconds()//3600)

    plt.title(titlestr+' '+cffield.description+' ('+cffield.units+')\n'+
             "Initialized: "+params['times'][0].strftime('%Y-%m-%d %H:%M')+"Z Forecast hour: "+'{:03d}'.format(timediff_secs)+" Valid: "+params['times'][params['time_index']].strftime('%Y-%m-%d %H:%M')+'Z')

    for key in landmark_dict.keys():
        kx,ky=bm(landmark_dict[key][0],landmark_dict[key][1])
        plt.text(kx,ky,key,fontsize=12,
                        ha='center',va='center',color='b')
    #fig.figimage(im, fig.bbox.xmax-290, fig.bbox.ymin,zorder=10)


#----
#    def draw_screen_poly( lats, lons, bm):
#             x, y = bm( lons, lats )
#             xy = zip(x,y)
#             poly = Polygon( xy, facecolor='black', alpha=0.45 )
#             plt.gca().add_patch(poly)
#            
#    lats = [ -33, -30, -30, -33 ]
#    lons = [ -65, -65, -62, -62 ]
#    draw_screen_poly( lats, lons, bm )

#----

#Drawing a 200-km range ring around S-PolKa site----
    def shoot(lon, lat, azimuth, maxdist=None):
        """Shooter Function
        Original javascript on http://williams.best.vwh.net/gccalc.htm
        Translated to python by Thomas Lecocq
        """
        glat1 = lat * np.pi / 180.
        glon1 = lon * np.pi / 180.
        s = maxdist / 1.852
        faz = azimuth * np.pi / 180.
 
        EPS= 0.00000000005
        if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
            alert("Only N-S courses are meaningful, starting at a pole!")
 
        a=6378.13/1.852
        f=1/298.257223563
        r = 1 - f
        tu = r * np.tan(glat1)
        sf = np.sin(faz)
        cf = np.cos(faz)
        if (cf==0):
            b=0.
        else:
            b=2. * np.arctan2 (tu, cf)
 
        cu = 1. / np.sqrt(1 + tu * tu)
        su = tu * cu
        sa = cu * sf
        c2a = 1 - sa * sa
        x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
        x = (x - 2.) / x
        c = 1. - x
        c = (x * x / 4. + 1.) / c
        d = (0.375 * x * x - 1.) * x
        tu = s / (r * a * c)
        y = tu
        c = y + 1
        while (np.abs (y - c) > EPS):
 
            sy = np.sin(y)
            cy = np.cos(y)
            cz = np.cos(b + y)
            e = 2. * cz * cz - 1.
            c = y
            x = e * cy
            y = e + e - 1.
            y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
                  d / 4. - cz) * sy * d + tu
 
        b = cu * cy * cf - su * sy
        c = r * np.sqrt(sa * sa + b * b)
        d = su * cy + cu * sy * cf
        glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
        c = cu * cy - su * sy * cf
        x = np.arctan2(sy * sf, c)
        c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
        d = ((e * cy * c + cz) * sy * c + y) * sa
        glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    
 
        baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
 
        glon2 *= 180./np.pi
        glat2 *= 180./np.pi
        baz *= 180./np.pi
 
        return (glon2, glat2, baz)


    def equi(bm, centerlon, centerlat, radius, *args, **kwargs):
        glon1 = centerlon
        glat1 = centerlat
        X = []
        Y = []
        for azimuth in range(0, 360):
            glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
            X.append(glon2)
            Y.append(glat2)
        X.append(X[0])
        Y.append(Y[0])
 
        #m.plot(X,Y,**kwargs) #Should work, but doesn't...
        X,Y = bm(X,Y)
        plt.plot(X,Y,**kwargs)

    radii = [200]
 
    # S-PolKa:
    centerlon = -63.978435
    centerlat = -31.715689
    for radius in radii:
        equi(bm, centerlon, centerlat, radius,lw=1.5,color='k')

    # San Rafael:
    centerlon = -68.285522
    centerlat = -34.587997
    for radius in radii:
        equi(bm, centerlon, centerlat, radius,lw=1.5,color='k')

    # Mendoza:
    centerlon = -68.7987
    centerlat = -32.8278
    for radius in radii:
        equi(bm, centerlon, centerlat, radius,lw=1.5,color='k')

#----

    os.system('mkdir -p '+outpath+'/'+params['modfld'])
    plt.savefig(outpath+'/'+params['modfld']+'/model.'+params['modname']+'.'+params['times'][0].strftime('%Y%m%d%H%M')+'.'+'{:03d}'.format(timediff_secs)+'_'+params['modfld']+'.png',dpi=225,bbox_inches='tight')


# In[512]:


#=====================PRECIPITABLE WATER=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'pwat',
        'cfield':'pw',
        'clevels':np.arange(15,75,5),
        'ccmap':"viridis_r",
        'llevels':None,
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}

cffield = getvar(files, params['cfield'], timeidx=params['time_index'])
lfield = None
lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
ufld = None
vfld = None

make_plot(cffield,lfield,lfield2,ufld,vfld,params)


# In[497]:


#=========================MIN UH=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'upheli',
        'cfield':'UP_HELI_MIN',
        'clevels':np.arange(0,220,20),
        'ccmap':"gist_stern_r",
        'llevels':None,
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}

try:
	cffield = getvar(files, params['cfield'], timeidx=params['time_index'])
	cffield.values=cffield.values*-1.
	cffield.attrs['description']='updraft helicity'
	cffield.attrs['units']='m2 s-2'
	lfield = None
	lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
	ufld = None
	vfld = None

	make_plot(cffield,lfield,lfield2,ufld,vfld,params)
except: print(params['cfield']," is not in the file, skipping")

# In[496]:

#=========================0-1 km SRH=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'0_1km_srh',
        'cfield':'SRH_MIN1',
        'clevels':np.arange(0,750,25),
        'ccmap':"cubehelix_r",
        'llevels':None,
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}

try:
	cffield = getvar(files, 'helicity', timeidx=params['time_index'],top=1000.0)
	cffield.values=cffield.values*-1.
	cffield.attrs['description']='0-1 km AGL storm relative helicity'
	cffield.attrs['units']='m2 s-2'
	lfield = None
	lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
	uvmet = getvar(files, 'uvmet10', timeidx=params['time_index'], units='kt')
	ufld = uvmet.isel(u_v=0)
	vfld = uvmet.isel(u_v=1)

	make_plot(cffield,lfield,lfield2,ufld,vfld,params)
except: print(params['cfield']," is not in the file, skipping")

# In[496]:

#=========================0-3 km SRH=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'0_3km_srh',
        'cfield':'SRH_MIN3',
        'clevels':np.arange(0,950,25),
        'ccmap':"cubehelix_r",
        'llevels':None,
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}

try:
	cffield = getvar(files, 'helicity', timeidx=params['time_index'],top=3000.0)
	cffield.values=cffield.values*-1.
	cffield.attrs['description']='0-3 km AGL storm relative helicity'
	cffield.attrs['units']='m2 s-2'
	lfield = None
	lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
	uvmet = getvar(files, 'uvmet10', timeidx=params['time_index'], units='kt')
	ufld = uvmet.isel(u_v=0)
	vfld = uvmet.isel(u_v=1)

	make_plot(cffield,lfield,lfield2,ufld,vfld,params)
except: print(params['cfield']," is not in the file, skipping")

#====================TOTAL PRECIPITAITON=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'precip_acc',
        'cfield':'RAINNC',
        'clevels':np.arange(0,165,10),
        'ccmap':"ocean_r",
        'llevels':None,
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}

cffield=None
cffield = getvar(files, params['cfield'], timeidx=params['time_index'])
cffield.attrs['description']='total precipitation'
cffield.attrs['units']='mm'
lfield = None
lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
ufld = None
vfld = None

make_plot(cffield,lfield,lfield2,ufld,vfld,params)


# In[495]:


#====================1 HOUR PRECIPITAITON=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'precip_1hr',
        'cfield':'PREC_ACC_NC',
        'clevels':np.arange(0,105,5),
        'ccmap':"ocean_r",
        'llevels':None,
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}
try:
	cffield=None
	cffield = getvar(files, params['cfield'], timeidx=params['time_index'])
	cffield.attrs['description']='hourly precipitation'
	cffield.attrs['units']='mm'
	lfield = None
	lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
	ufld = None
	vfld = None

	make_plot(cffield,lfield,lfield2,ufld,vfld,params)
except: print(params['cfield']," is not in the file, skipping")

# In[503]:


#====================MAX SFC WIND=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'max_sfc_wind',
        'cfield':'WSPD10MAX',
        'clevels':np.arange(0,100,2),
        'ccmap':"gist_stern_r",
        'llevels':None,
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}

try:
	cffield = getvar(files, params['cfield'], timeidx=params['time_index'], units='kt')
	cffield.values=cffield.values*1.94384
	#cffield.attrs['description']='maximum surface wind'
	cffield.attrs['units']='kt'
	lfield = None
	lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
	ufld = None
	vfld = None

	make_plot(cffield,lfield,ufld,vfld,params)
except: print(params['cfield']," is not in the file, skipping")


# In[504]:


#====================MAX UPDRAFT=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'max_w',
        'cfield':'W_UP_MAX',
        'clevels':np.arange(0,100,2),
        'ccmap':"gist_stern_r",
        'llevels':None,
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}

try:
	cffield = getvar(files, params['cfield'], timeidx=params['time_index'], units='kt')
	cffield.values=cffield.values*1.94384
	cffield.attrs['description']='maximum updraft velocity'
	cffield.attrs['units']='kt'
	lfield = None
	lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
	ufld = None
	vfld = None

	make_plot(cffield,lfield,lfield2,ufld,vfld,params)
except: print(params['cfield']," is not in the file, skipping")


# In[553]:


#=====================2 M TEMP=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'T2',
        'cfield':'T2',
        'clevels':np.arange(-15,40,2.5),
        'ccmap':"plasma",
        'llevels':np.arange(970,1040,4),
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}

cffield = getvar(files, params['cfield'], timeidx=params['time_index'],units='degC')
cffield.attrs['description']='2 m temperature'
cffield.attrs['temperature']='degC'
cffield.values=cffield.values-273.15
lfield = getvar(files, 'slp', timeidx=params['time_index'], units='hPa')
lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
uvmet = getvar(files, 'uvmet10', timeidx=params['time_index'], units='kt')
ufld = uvmet.isel(u_v=0)
vfld = uvmet.isel(u_v=1)

make_plot(cffield,lfield,lfield2,ufld,vfld,params)


# In[378]:


#=====================2 M DWPT=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'TD2',
        'cfield':'td2',
        'clevels':np.arange(-5,30,2.5),
        'ccmap':"YlGn",
        'llevels':np.arange(970,1040,4),
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}

cffield = getvar(files, params['cfield'], timeidx=params['time_index'])
cffield.attrs['description']='2 m dewpoint temperature'
cffield.attrs['temperature']='degC'
cffield.values=cffield.values
lfield = getvar(files, 'slp', timeidx=params['time_index'], units='hPa')
lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
uvmet = getvar(files, 'uvmet10', timeidx=params['time_index'], units='kt')
ufld = uvmet.isel(u_v=0)
vfld = uvmet.isel(u_v=1)

make_plot(cffield,lfield,lfield2,ufld,vfld,params)


# In[555]:


#=====================SFC THETAE=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'theta_e_sfc',
        'cfield':'eth',
        'clevels':np.arange(290,390,10),
        'ccmap':"viridis_r",
        'llevels':np.arange(970,1040,4),
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}

cffield = getvar(files, params['cfield'], timeidx=params['time_index'], units='K')
cffield=cffield.isel(bottom_top=0)
lfield = getvar(files, 'slp', timeidx=params['time_index'], units='hPa')
lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
uvmet = getvar(files, 'uvmet10', timeidx=params['time_index'], units='kt')
ufld = uvmet.isel(u_v=0)
vfld = uvmet.isel(u_v=1)

make_plot(cffield,lfield,lfield2,ufld,vfld,params)


# In[556]:


#=====================10 M WIND=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'WSP',
        'cfield':'wsp',
        'clevels':np.arange(8,40,4),
        'ccmap':"YlOrBr",
        'llevels':np.arange(970,1040,4),
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}


uvmet10 = getvar(files, 'wspd_wdir10', timeidx=params['time_index'], units="kt")

cffield = uvmet10.isel(wspd_wdir=0)
cffield.attrs['description']='10 m wind'
lfield = getvar(files, 'slp', timeidx=params['time_index'], units='hPa')
lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
uvmet = getvar(files, 'uvmet10', timeidx=params['time_index'], units='kt')
ufld = uvmet.isel(u_v=0)
vfld = uvmet.isel(u_v=1)

make_plot(cffield,lfield,lfield2,ufld,vfld,params)


# In[415]:


#=====================0-1 KM WIND=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'lcl_0-1kmwind',
        'cfield':'lcl_0-1kmwind',
        'clevels':np.arange(500,4000,250),
        'ccmap':"ocean",
        'llevels':np.arange(10,50,5),
        'llevels2':[1000],
        'lalpha':0.7,
        'time_index':tindex,
        'times':times,
        'skip':17}


uvmet10 = getvar(files, 'wspd_wdir10', timeidx=params['time_index'], units="kt")

ter = getvar(files, 'ter', timeidx=params['time_index'], units="m")
z = getvar(files, 'z', timeidx=params['time_index'], units="m")
p = getvar(files, 'p', timeidx=params['time_index'], units="hPa")
u = getvar(files, 'ua', timeidx=params['time_index'], units="kt")
v = getvar(files, 'va', timeidx=params['time_index'], units="kt")

ter_3d=np.tile(ter.values,[41,1,1])
z.values=z.values-ter.values
u2 = interplevel(u, z, 1000)
v2 = interplevel(v, z, 1000)
uvmet = getvar(files, 'uvmet10', timeidx=params['time_index'], units='kt')
ufld = uvmet.isel(u_v=0)
vfld = uvmet.isel(u_v=1)

u2.values=u2.values-ufld.values
v2.values=v2.values-vfld.values

lfield = uvmet10.isel(wspd_wdir=0)
lfield.values=smooth2d(np.sqrt(u2.values**2.+v2.values**2.),3)

lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')

sounding_parameters=getvar(files, 'cape_2d', timeidx=params['time_index'])
cffield = sounding_parameters.isel(mcape_mcin_lcl_lfc=3)

cffield.attrs['description']='lcl, 0-1 km bulk wind difference'
cffield.attrs['units']='m; kt'
uvmet = getvar(files, 'uvmet10', timeidx=params['time_index'], units='kt')
ufld = u2
vfld = v2

make_plot(cffield,lfield,lfield2,ufld,vfld,params)


# In[414]:


#=====================MLCAPE; 0-6 KM WIND=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'mlcape_0-6kmwind',
        'cfield':'mlcape_0-6kmwind',
        'clevels':np.arange(0,4000,200),
        'ccmap':"viridis_r",
        'llevels':np.arange(10,50,5),
        'llevels2':[1000],
        'lalpha':0.7,
        'time_index':tindex,
        'times':times,
        'skip':17}


uvmet10 = getvar(files, 'wspd_wdir10', timeidx=params['time_index'], units="kt")

ter = getvar(files, 'ter', timeidx=params['time_index'], units="m")
z = getvar(files, 'z', timeidx=params['time_index'], units="m")
p = getvar(files, 'p', timeidx=params['time_index'], units="hPa")
u = getvar(files, 'ua', timeidx=params['time_index'], units="kt")
v = getvar(files, 'va', timeidx=params['time_index'], units="kt")

ter_3d=np.tile(ter.values,[41,1,1])
z.values=z.values-ter.values
u2 = interplevel(u, z, 6000)
v2 = interplevel(v, z, 6000)
uvmet = getvar(files, 'uvmet10', timeidx=params['time_index'], units='kt')
ufld = uvmet.isel(u_v=0)
vfld = uvmet.isel(u_v=1)

u2.values=u2.values-ufld.values
v2.values=v2.values-vfld.values

lfield = uvmet10.isel(wspd_wdir=0)
lfield.values=smooth2d(np.sqrt(u2.values**2.+v2.values**2.),3)

lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')

sounding_parameters=getvar(files, 'cape_2d', timeidx=params['time_index'])
cffield = sounding_parameters.isel(mcape_mcin_lcl_lfc=0)

cffield.attrs['description']='MLCAPE, 0-6 km bulk wind difference'
cffield.attrs['units']='J kg-1; kt'
uvmet = getvar(files, 'uvmet10', timeidx=params['time_index'], units='kt')
ufld = u2
vfld = v2

make_plot(cffield,lfield,lfield2,ufld,vfld,params)


# In[338]:


#=========================500 M=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'wsp_500m',
        'cfield':'wsp_500m',
        'clevels':np.arange(8,40,4),
        'ccmap':"YlOrBr",
        'llevels':np.arange(510,606,6),
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}

uvmet = getvar(files, 'wspd_wdir', timeidx=params['time_index'], units="kt")
wspd = uvmet.isel(wspd_wdir=0)
ter = getvar(files, 'ter', timeidx=params['time_index'], units="m")
z = getvar(files, 'z', timeidx=params['time_index'], units="m")
p = getvar(files, 'p', timeidx=params['time_index'], units="hPa")
u = getvar(files, 'ua', timeidx=params['time_index'], units="kt")
v = getvar(files, 'va', timeidx=params['time_index'], units="kt")

ter_3d=np.tile(ter.values,[41,1,1])
z.values=z.values-ter.values
cffield = interplevel(wspd, z, 500)
cffield.attrs['description']='500 m AGL winds'
lfield = None
lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
ufld = interplevel(u, z, 500)
vfld = interplevel(v, z, 500)

make_plot(cffield,lfield,lfield2,ufld,vfld,params)


# In[481]:


#=========================1 km DBZ=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'radar_1km',
        'cfield':'radar_1km',
        'clevels':np.arange(5,75,5.),
        'ccmap':"gist_ncar",
        'llevels':np.arange(510,606,6),
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}

dbz = getvar(files, 'dbz', timeidx=params['time_index'])
ter = getvar(files, 'ter', timeidx=params['time_index'], units="m")
z = getvar(files, 'z', timeidx=params['time_index'], units="m")

ter_3d=np.tile(ter.values,[41,1,1])
z.values=z.values-ter.values
cffield = interplevel(dbz, z, 1000)
cffield.values=np.ma.masked_less(cffield.values,5.)
cffield.attrs['description']='1 km AGL radar reflectivity'
lfield = None
lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
ufld=None
vfld=None

make_plot(cffield,lfield,lfield2,ufld,vfld,params)


# In[418]:


#=========================500 HPA=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'500avo',
        'cfield':'avo',
        'clevels':[-12,-9,-6,-3,0,12,15,18,21,24],
        'ccmap':"RdBu",
        'lalpha':0.7,
        'llevels':np.arange(510,606,6),
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}

avo = getvar(files, 'avo', timeidx=params['time_index'])
z = getvar(files, 'z', timeidx=params['time_index'], units="dm")
p = getvar(files, 'p', timeidx=params['time_index'], units="hPa")
u = getvar(files, 'ua', timeidx=params['time_index'], units="kt")
v = getvar(files, 'va', timeidx=params['time_index'], units="kt")

cffield = interplevel(avo, p, 500)
cffield.values=cffield.values*-1
cffield.attrs['description']='500 hPa absolute vorticity'
lfield = interplevel(z, p, 500)
lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
ufld = interplevel(u, p, 500)
vfld = interplevel(v, p, 500)

make_plot(cffield,lfield,lfield2,ufld,vfld,params)


# In[510]:


#=========================200 HPA=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'200_V',
        'cfield':'200_V',
        'clevels':np.arange(40,160,20),
        'ccmap':"plasma_r",
        'lalpha':0.7,
        'llevels':np.arange(1100,1300,12),
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}

uvmet = getvar(files, 'wspd_wdir', timeidx=params['time_index'], units="kt")
wspd = uvmet.isel(wspd_wdir=0)
z = getvar(files, 'z', timeidx=params['time_index'], units="dm")
p = getvar(files, 'p', timeidx=params['time_index'], units="hPa")
u = getvar(files, 'ua', timeidx=params['time_index'], units="kt")
v = getvar(files, 'va', timeidx=params['time_index'], units="kt")

cffield = interplevel(wspd, p, 200)
cffield.values=cffield.values
cffield.attrs['description']='200 hPa isotachs'
cffield.values=np.ma.masked_less(cffield.values,40.)
lfield = interplevel(z, p, 200)
lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
ufld = interplevel(u, p, 200)
vfld = interplevel(v, p, 200)

make_plot(cffield,lfield,lfield2,ufld,vfld,params)


# In[546]:


#=========================850 HPA=========================
params={'outpath':outpath,
        'modname':modname,
        'modfld':'850_theta',
        'cfield':'850_theta',
        'clevels':np.arange(-300,340,40),
        'ccmap':"BrBG",
        'lalpha':0.7,
        'llevels':np.arange(250,400,5),
        'llevels2':[1000],
        'time_index':tindex,
        'times':times,
        'skip':17}

uvmet = getvar(files, 'wspd_wdir', timeidx=params['time_index'], units="kt")
wspd = uvmet.isel(wspd_wdir=0)
eth = getvar(files, 'eth', timeidx=params['time_index'])
q = getvar(files, 'QVAPOR', timeidx=params['time_index'])

dx=uvmet.XLAT[1,0]-uvmet.XLAT[0,0]
dx=dx.values*111100
#MFC_advect = -( u*(dq/dx)+v*(dq/dy) )    ; advection term
#    MFC_conv    = -q*( (du/dx)+  (dv/dy) )      ; con(div)-vergence
#    MFC = MFC_advect + MFC_conv

p = getvar(files, 'p', timeidx=params['time_index'], units="hPa")
u = getvar(files, 'ua', timeidx=params['time_index'])
v = getvar(files, 'va', timeidx=params['time_index'])
ufld = interplevel(u, p, 850)
vfld = interplevel(v, p, 850)
qfld = interplevel(q, p, 850)

grad_q_x,grad_q_y = np.gradient(qfld.values)
grad_u_x,grad_u_y = np.gradient(ufld.values)
grad_v_x,grad_v_y = np.gradient(vfld.values)

MFC_advect=-1.*(ufld.values*grad_q_x/dx)+(vfld.values*grad_q_y/dx)
MFC_conv=-1.*qfld.values*((grad_u_x/dx)+(grad_v_y/dx))

cffield = qfld
cffield.values=-1.*smooth2d(86400.*1000.*(MFC_advect + MFC_conv),5)
cffield.attrs['description']='850 hPa moisture convergence and theta-e'
cffield.attrs['units']='g kg-1 dy-1; K'
lfield = smooth2d(interplevel(eth, p, 850),3)
lfield2 = getvar(files, 'ter', timeidx=params['time_index'], units='m')
ufld = ufld*1.94
vfld = vfld*1.94

make_plot(cffield,lfield,lfield2,ufld,vfld,params)


# In[547]:



# In[ ]:




