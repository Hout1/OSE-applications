#!/usr/bin/env python
# coding: utf-8

# # Extreme Rainfall Detection and Alerting System (ERDAS) v 0.1

# In[1]:


__author__ = 'OSE GIS'
__project__ = 'Extreme Rainfall Detection and Alerting System (ERDAS)'
__contact__ = 'michael.manalili@wfp.org', 'wfp.hq.gis@wfp.org'


# In[2]:


#!/usr/bin/env python3
import requests, os, time, shutil
from time import strptime
import matplotlib.pyplot as plt
from datetime import datetime, timedelta, date
import pandas as pd
import folium
import geopandas as gpd
import descartes
import schedule
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart #python3
from email.mime.base import MIMEBase #python3
from email import encoders #python3
from sentinelsat.sentinel import SentinelAPI, read_geojson, geojson_to_wkt
import ee
import eeconvert
import rasterio
from rasterio.plot import show
from rasterstats import zonal_stats
import cdsapi
from geopandas import GeoSeries, GeoDataFrame
from sqlalchemy import create_engine
import psycopg2 
import io
from lxml import html
import wget
import cdsapi
from configparser import ConfigParser
import sentinelsat
from sqlalchemy import create_engine
import psycopg2 
import pandas.io.sql as psql
from io import BytesIO as StringIO
from osgeo import ogr
from ftplib import FTP
import pycountry
import uuid
#import pysftp
#from io import BytesIO as StringIO

try:
    ee.Initialize()
    print ('Earth Engine package initialized successfully..')
except ee.EEException as e:
    print ('Earth Engine package failed to initialize!')
except:
    print ('Unexpected error:', sys.exc_info()[0])
    raise

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

pd.set_option('mode.chained_assignment', None)
#pd.set_option('mode.chained_assignment', 'raise')

execTime = time.time()


# In[3]:


#import geopandas as gpd
#dir(gpd.GeoDataFrame)


# In[4]:


now = datetime.today()
d = str(now)
date = d[0:10]
date_cut = date.replace('-', '')
whitelist_data = 'extreme_precip_alert_' + date_cut
alert_data = 'red_alert_' + date_cut

script = os.path.dirname(os.path.realpath("__file__"))
root = os.path.dirname(script)

# Assumed this folders are pre-created
process_root = root + '/processing'
process_amsr2 = process_root + '/process_amsr2/'
process_cdsapi = process_root + '/process_cdsapi/'
process_ecmwf = process_root + '/process_ecmwf/'
process_pxr = process_root + '/process_pxr/'

alert_csv_path = root + '/alert/'
s1_txt_path = root + '/alert/'
map_view_path = root + '/view/'

alert_data_name = 'Extreme_rainfall_alert_beta_' + date_cut + '.csv'
map_view_name = 'Map_view_' + date_cut + '.html'
s1_name = 'sentinel_sched_' + date_cut + '.txt'


# In[5]:


config = ConfigParser()
config.read(os.path.join(root, "config.txt"))

esa_un = config.get('sentinelsat', 'ss_un')
esa_pw = config.get('sentinelsat', 'ss_pw')
esa_link = config.get('sentinelsat', 'ss_link')
api = SentinelAPI(esa_un, esa_pw, esa_link)
api

config.read(os.path.join(root, "config.txt"))
smtp_un = config.get('dbtrack','dbt_uname')
smtp_pw = config.get('dbtrack','dbt_pword')
print('config loaded.')


# In[6]:


#ECMWF_SDI_raster_cat = config.get("ECMWF","ECMWF_SDI_raster_cat")
#ECMWF_SDI_raster_cat = os.path.join(connfiles_folder,ECMWF_SDI_raster_cat)

#ECMWF_SDI_mosaic_dataset= config.get("ECMWF","ECMWF_SDI_mosaic_dataset")
#ECMWF_SDI_mosaic_dataset = os.path.join(connfiles_folder,ECMWF_SDI_mosaic_dataset)


# In[7]:


def spatial_red(ans, reducer, scale):
    data = ans.reduceRegions(
            collection = wfp_aoi,
            reducer = reducer,
            crs='EPSG:3857',
            scale=scale,
            tileScale=1)
    return data   

def write2file(filename, header, gdf):
    data = gdf.to_csv(filename, columns=header, encoding='utf-8')
    return data

def export2drive_raster(raster, bbox ,prefix,scale):
    run = ee.batch.Export.image.toDrive(
         image=raster, 
         scale=scale,
         fileNamePrefix = prefix,
         region = ee.Feature(bbox).geometry().bounds().getInfo()['coordinates'], #can also be single country
         maxPixels = 1e13,
         fileFormat='GeoTiff',
         cloudOptimized=True).start()
    return run

#Format Options: GeoJSON, SHP, KML
def export2drive_table(ee_featCol, description):
    table = ee.batch.Export.table.toDrive(
         collection=ee_featCol, 
         description=description,
         fileFormat='CSV')#.start() 
    return table

def month(n):
    start = ee.Date('2014-01-01').advance(n, 'month')
    end = start.advance(1,'month')
    return ee.ImageCollection('JAXA/GPM_L3/GSMaP/v6/operational').select('hourlyPrecipRate').filterDate(start,end).mean().set('system:time_start', start.millis())

def get_wc():
    #x = ee.Image(collectionList.get(band_month))
    wc = ee.Image('users/wfphqgis/CLIM/wc20_30s_prec_'+s_bm)
    wc_vector_mean = spatial_red(wc,ee.Reducer.mode(),900)
    df_wc = eeconvert.fcToGdf(wc_vector_mean)
    dff_wc = df_wc.drop(columns = ['adm0_name', 'adm1_name','adm2_name', 'adm0_code', 'adm1_code','disp_area','geometry','status'])
    dff_wc_r = dff_wc.rename(columns={'mode':'WC_Mean_Prec'})
    return dff_wc_r

def get_sm2map():
    #x = ee.Image(collectionList.get(band_month))
    sm2 = ee.Image('users/wfphqgis/CLIM/clm_precipitation_sm2rain_m_1km_s00cm_20072018_v02_'+s_bm)
    sm2_vector_mean = spatial_red(sm2,ee.Reducer.mode(),900)
    df_sm2 = eeconvert.fcToGdf(sm2_vector_mean)
    dff_sm2 = df_sm2.drop(columns = ['adm0_name', 'adm1_name','adm2_name', 'adm0_code', 'adm1_code','disp_area','geometry','status'])
    dff_sm2_r = dff_sm2.rename(columns={'mode':'Mean_Prec'})
    return dff_sm2_r

def get_ecmwf():
    data_ecmwf = ee.Image('users/wfphqgis/CLIM/ecmwf_reanalysis_ymonmeanR_ymonmean')
    crs = data_ecmwf.projection().crs()
    scale = 9000
    tp = data_ecmwf.select('b' + s_bm)
    ecmwf = tp.resample('bilinear').reproject(crs = crs, scale = scale)
    wc_vector_mean = spatial_red(tp,ee.Reducer.mode(),9000)
    df_wc = eeconvert.fcToGdf(wc_vector_mean)
    dff_wc = df_wc.drop(columns = ['adm0_name', 'adm1_name','adm2_name', 'adm0_code', 'adm1_code','disp_area','geometry','status'])
    dff_wc_r = dff_wc.rename(columns={'mode':'ECMWF_Mean_Prec'})
    return dff_wc_r

def get_amsr2():
    #data = 'http://www.gdacs.org/flooddetection/DATA/AMSR2/AvgSignalTiffs/2019' 400K max
    data = 'http://www.gdacs.org/flooddetection/DATA/AMSR2/MagTiffs/2019' #20K max
    now = datetime.today()
    d = str(now)
    date = d[0:10]
    date_frmt = date.replace('-', '')

    page = requests.get(data)
    webpage = html.fromstring(page.content)

    amsr2 = webpage.xpath('//a/@href')
    amsr2_last = amsr2[-1]
    cur_year = d[0:4]
    base = 'http://www.gdacs.org'
    #prod = 'signal_4days_avg_4days_'
    prod = 'mag_signal_'
    hd5 = '_HD5'
    ext = '.tif'
    lnk =  base + amsr2_last + prod + date_frmt +hd5 +ext
    dl = wget.download(lnk,process_amsr2)
    return dl

def get_s1_info(i):
    s1 = api.query(final_df.envelope[i],
                   date=sensing_date,
                   platformname='Sentinel-1',
                   producttype='GRD',
                   area_relation='Intersects') #Intersects, IsWithin
                   #orbitdirection='ASCENDING')
    s1_json = api.to_geojson(s1)
    s1_gdf = api.to_geodataframe(s1)
    s1_gdf

    s1_gdf = s1_gdf[['ingestiondate']]#,'beginposition','ingestiondate','platformname','orbitnumber', 'producttype', 'endposition']]
    s1_gdf.sort_values('ingestiondate', ascending=True).head()
    s1_gdf

    dates = s1_gdf["ingestiondate"].dt.date.values
    map_date = dates + timedelta(days=11)
    ld = len(dates)
    
    print_path = open(root + '/alert/' + 'sentinel_sched_' + date_cut + '.txt','a')
    
    p = print("Sentinel Acquisition Schedule: " + "{0}: {1} \n{2} alert! \nThere are {3} Sentinel-1 overpass \nMappable on:\n{4}\n"
             .format(final_df['adm0_name'][i],final_df['adm2_name'][i],final_df['status'][i],ld,map_date),
             file=print_path)
    return p

def job():
    msg = MIMEMultipart()
    msg['From'] = 'wfp.hq.dbtrack@gmail.com'
    msg['Subject'] = "Extreme Precipitation Alert"
    
    body = """<div style="color:rgb(28, 130, 196)"><font size="+2"><b><img id="myimage" src="https://mw1.google.com/crisisresponse/icons/un-ocha/disaster_heavy_rain_100px_icon.png"></img>
    <br /> Extreme Precipitation Alert</b> </font></div>
    <br/><b>Please See attached list and location of WFP country concerned.</b></b>
    <br />
    <br/> Number of <font color="YELLOW">YELLOW</font> Alerts: %s
    <br/> Number of <font color="ORANGE">ORANGE</font> Alerts: %s
    <br/> Number of <font color="RED">RED</font> Alerts: %s
    <br />
    <br/><b>Note on color categories:</b></b>
    <br />
    <br/>Yellow: The 3-days accumulated rainfall is greater than the average 3 days maximum total precipitation (moderate nowcast)
    <br/>Orange: The 3-days accumulated rainfall is greater than the long term mean rainfall for that month (high nowcast)
    <br/>Red*: The 3-days accumulated rainfall + 3-days forecast is greater than the long term mean rainfall (high forecast)
    <br/><i>*The red alert is valid until 3-5 days of forecast</i></b>
    <br />
    <br/><b>Emergency Mapping/Charter Activation Links:</b></b>
    <br />
    <br/><b>--------------------------------------------------------------------------------</b></b>
    <br/>This is an automatically generated email, please do not reply.<br/>
    <br/>Service provided by 
    <br/><i>UNITED NATIONS WORLD FOOD PROGRAMME
    <br/>WFP Emergency Preparedness & Support Response Division
    <br/>Contact: <a href="mailto:hq.gis@wfp.org">HQ Geospatial Support Unit</a></i></p>
    """ % (y,o,r) #add 'a' in the modulo

    alert_fn = alert_csv_path + alert_data_name
    mapview_fn = map_view_path + map_view_name
    s1_fn = s1_txt_path + s1_name
    
    alerts = MIMEBase('application', "octet-stream")
    mapview = MIMEBase('application', "octet-stream")
    s1view = MIMEBase('application', "octet-stream")
    
    alerts.set_payload(open(alert_fn, "rb").read())
    mapview.set_payload(open(mapview_fn, "rb").read())
    s1view.set_payload(open(s1_fn, "rb").read())
    
    encoders.encode_base64(alerts)
    encoders.encode_base64(mapview)
    encoders.encode_base64(s1view)
    
    alerts.add_header('Content-Disposition', 'attachment; filename=%s'%alert_data_name)
    mapview.add_header('Content-Disposition', 'attachment; filename=%s'%map_view_name)
    s1view.add_header('Content-Disposition', 'attachment; filename=%s'%s1_name)
    
    msg.attach(alerts)
    msg.attach(mapview)
    msg.attach(s1view)
    
    msg.attach(MIMEText(body, 'html'))
    
#     s = smtplib.SMTP('smtp.gmail.com', 587) #Anywhere
#     s.starttls() #Anywhere
    
    s = smtplib.SMTP_SSL('smtp.gmail.com', 465) # WFP domain
    
    s.ehlo()
    s.login(smtp_un, smtp_pw)
    sender = 'wfp.hq.dbtrack@gmail.com'
    
    #recipients = ['michaelandrew.manalili@gmail.com']#,'michael.manalili@wfp.org','abdel-lathif.younous@wfp.org','rohini.swaminathan@wfp.org']
    recipients = ['rohini.swaminathan@wfp.org','sirio.modugno@wfp.org','stefano.cairo@wfp.org','michael.manalili@wfp.org','abdel-lathif.younous@wfp.org','michaelandrew.manalili@gmail.com']
            
    if gbl_alert.empty == True:
        pass
        print('DataFrame is Empty, skipping email broadcast...')
        
    else:
        s.sendmail(sender, recipients, str(msg))
        s.quit()
        print('ALERT sent to subscribers!')

#Flood: https://mw1.google.com/crisisresponse/icons/un-ocha/disaster_flood_100px_icon.png
#http://www.gdacs.org/flooddetection/DATA/AMSR2/
#http://www.gdacs.org/flooddetection/DATA/AMSR2/SignalTiffs/2019/
#http://www.gdacs.org/flooddetection/DATA/AMSR2/MagTiffs/2019/

# if mode == 'prod':
#             mailServer = smtplib.SMTP("smtp.gmail.com", 587)
#         else:### in virtual machines security rules don't allow anymore using  (old version) TLS connection, but only SSL. In production, SSL version is good enough
#             mailServer = smtplib.SMTP_SSL("smtp.gmail.com", 465)

### Add the below to the message to include GIS focal points
# <br /><b>GIS Focal point for coordinating the event:</b>
# <br />
# </b>%s


# In[8]:


iso_global = pd.read_csv(root + '/wfp' + '/wfp_country.csv', index_col=0)
iso_global_list = iso_global['ISO3'].tolist()


# In[9]:


# Add population to:
#Forecast + XDAYS + Population (Abdel to advise on this number)


# In[10]:


GBL_adm0 = ee.FeatureCollection("users/wfphqgis/BND/wfp_bnd_inform2019") #iso3
GBL_adm2 = ee.FeatureCollection("users/wfphqgis/BND/admin_2_gaul_2015") #admX_name
#GBL_dfo = ee.FeatureCollection("users/wfphqgis/FLOOD/DFO_HistoricalFloodEvents")
GBL_dfo = ee.FeatureCollection("users/wfphqgis/FLOOD/DFO_HistoricalFloodEvents_M")

wfp_adm0 = GBL_adm0.filter(ee.Filter.inList
                           ('iso3',['DRC','NIC','HND','GTM','CAF','IRQ','LBY','BGD','MYM','COL','CIV',
                                    'BFA','SSD','YEM','SYR','CAM','NIG','DRC','MOZ','PHL','LBN','LKA',
                                    'ZWE','BEN','TGO','GHA','TUR','IND','MLI','PAN','SOM','TCD','NGA',
                                    'IDN','GIN','PAK','TZA','IRN','TJK','AFG','IRQ','SLE','VEN','VNM',
                                    'CMR','GMB','GNB','LBR','MRT','NER','SDN','MRT','MMR','LAO']))


# L3
# ['SSD','YEM','SYR','CAM','NIG','DRC','MOZ']
# L2
# ['MLI','CAF','IRQ','LBY','BGD','MYM','COL', 'BFA']
# M
# ['DRC','NIC','HND','GTM']
# Other
# ['DRC','NIC','HND','GTM','MLI','CAF','IRQ','LBY','BGD','MYM','COL',
# 'BFA','SSD','YEM','SYR','CAM','NIG','DRC','MOZ','PHL','LBN','LKA',
# 'ZWE','BEN','TGO','GHA','TUR','IND','MLI','PAN','SOM','TCD','NGA',
# 'IDN','GIN','PAK','TZA','IRN','TJK','AFG','IRQ','SLE']

# Should work globally..
# wfp_adm0 = GBL_adm0.filter(ee.Filter.inList
#                            ('iso3',iso_global_list))

# Togo, Benin, Ghana, Mali, Codivor, SSD, CAR, PANAMA

dfo = GBL_dfo.filterBounds(wfp_adm0)
wfp_aoi = GBL_adm2.filterBounds(dfo)
dfodata = eeconvert.fcToGdf(GBL_dfo)


# ## GEE data calls

# In[11]:


GSMaP = ee.ImageCollection('JAXA/GPM_L3/GSMaP/v6/operational')
precipitation = GSMaP.select('hourlyPrecipRate')

worldpop = ee.ImageCollection('WorldPop/POP').filter(ee.Filter.equals('year', 2015)).filter(ee.Filter.equals('UNadj', 'yes'))
wpop2015 = worldpop.select('population').reduce(ee.Reducer.max()).clip(wfp_adm0)

ciesn2020 = ee.ImageCollection('CIESIN/GPWv4/unwpp-adjusted-population-count')
ciesn2020_pop = ciesn2020.select('population-count').reduce(ee.Reducer.max()).clip(wfp_adm0)

jrc20yrp = ee.Image('users/wfphqgis/HAZ/floodMapGL_rp20y')
jrc_1m_remap = jrc20yrp.lte(14)
jrc_1m_depth = jrc_1m_remap.remap([1,2,3,4,5,6,7,8,9,10,11,12,13,14],
                                  [1,1,1,1,1,1,1,1,1,1,1,1,1,1])

hectares = ee.Image.pixelArea().divide(10000)

hectares_jrc = jrc20yrp.multiply(hectares)

#pop_haz_area = wpop2015.multiply(jrc_1m_depth)
pop_haz_area = ciesn2020_pop.multiply(jrc_1m_depth)

#Change
pxr = ee.Image('users/wfphqgis/CLIM/Global_IDF')
pxr_3days = pxr.select('b12') #band 12 for 72hours (3Days) band 14 is 120hours (5days) and 
pxr_3days_vector = spatial_red(pxr_3days,ee.Reducer.mode(),30000)

#Uncomment for GFS data (update datetime is not yet fixed)
#gfs_fc = ee.ImageCollection('NOAA/GFS0P25').select("total_precipitation_surface").filter(ee.Filter.eq('creation_time',ee.Date(0).update(2019,7,21,6,0,0)
#                                                                                                      .millis())).filter(ee.Filter.eq('forecast_hours',72)).sum()
#gfs_fc72h_res = gfs_fc.resample('bilinear').reproject(crs=gfs_fc.projection(), scale=10000).clip(wfp_adm0)

#vam = ee.Image('users/wfphqgis/CLIM/wld_cli_rainfall_threshold_q96_25yr_chirps_wfp')
#vam_vector = spatial_red(pxr_3days,ee.Reducer.mode(),5000)


# In[12]:


#GEE Datetime for GPM accumulated Precipitation
day0 = datetime.today()
hrs24 = day0 - timedelta(days=1)
hrs48 = day0 - timedelta(days=2)
hrs72 = day0 - timedelta(days=3)
day7 = day0 - timedelta(days=7)

#Temporal Reducer
gpm_precip_24h = precipitation.filter(ee.Filter.date(hrs24, day0))
gpm_precip_48h = precipitation.filter(ee.Filter.date(hrs48, day0))
gpm_precip_72h = precipitation.filter(ee.Filter.date(hrs72, day0))
gpm_precip_day7 = precipitation.filter(ee.Filter.date(day7, day0))


# In[ ]:


#Changed to WorldClimV2
band_month = datetime.today().strftime("%m")
#band_month = strptime(run_date,'%b').tm_mon
s_bm = str(band_month)
dff_wc_rename = get_wc()
#dff_ecmwf_rename = get_ecmwf()
dff_sm2_rename = get_sm2map()
print('Climate Data loaded..')


# In[ ]:


#Spatial Reducer
gpm72h = gpm_precip_72h.sum().rename('72h').clip(wfp_aoi)
gsmap_vector_mean = spatial_red(gpm72h,ee.Reducer.mode(),500)
wpop_vector = spatial_red(wpop2015,ee.Reducer.sum(),95)
jrc_vector = spatial_red(jrc_1m_depth,ee.Reducer.sum(),1000) #1km resolution all models
exp_vector = spatial_red(pop_haz_area,ee.Reducer.sum(),1000) #95 for WorldPop and 1000 for CIESN

#Change
ciesn_vector = spatial_red(ciesn2020_pop,ee.Reducer.sum(),1000)
#gfs_fc72h_vector = spatial_red(gfs_fc72h_res,ee.Reducer.mode(),10000)


# In[ ]:


print('Computing Global parameters. Please wait...')
df_precip_m = eeconvert.fcToGdf(gsmap_vector_mean)
df_pop = eeconvert.fcToGdf(wpop_vector)
df_jrc = eeconvert.fcToGdf(jrc_vector)
df_exp = eeconvert.fcToGdf(exp_vector)

#Change
df_pxr = eeconvert.fcToGdf(pxr_3days_vector)
df_pop1km = eeconvert.fcToGdf(ciesn_vector)
#df_gfs = eeconvert.fcToGdf(gfs_fc72h_vector)


# In[ ]:


#ECMWF datetime format for FTP (Temporary Solution while migrating to Windows)
forecast = datetime.today()
fc_fmt = forecast.strftime("%m%d")
days3 = forecast + timedelta(days=3)
days5 = forecast + timedelta(days=5)
days7 = forecast + timedelta(days=7)

days3F = days3.strftime("%m%d")
days5F = days5.strftime("%m%d")
days7F = days7.strftime("%m%d")

ecmwf_3D_fc = fc_fmt + days3F + '00' + 'TP' + '.tiff'
ecmwf_5D_fc = fc_fmt + days5F + '00' + 'TP' + '.tiff'
ecmwf_7D_fc = fc_fmt + days7F + '00' + 'TP' + '.tiff'

config.read(os.path.join(root, "config.txt"))
ftp_user = config.get('gisftp', 'ftp_un')
ftp_pw = config.get('gisftp', 'ftp_pw')
ftp_gis =  config.get('gisftp', 'ftp_url')
ftp = FTP(ftp_gis)
ftp.login(user=ftp_user, passwd = ftp_pw)
save2local = ftp.cwd("/ECMWF_processed")

def grabFile(fn):
    filename = fn
    localfile = open(process_ecmwf + filename, 'wb')
    return ftp.retrbinary('RETR ' + filename, localfile.write, 1024)

grabFile(ecmwf_3D_fc)
grabFile(ecmwf_5D_fc)
grabFile(ecmwf_7D_fc)

ftp.quit()


# ### Alert data generation

# In[ ]:


dff_precip = df_precip_m.drop(columns = ['adm0_code', 'adm1_code','disp_area','status'])
dff_p_rename = dff_precip.rename(columns={'mode':'3D_acc_precip'})
dff_pop = df_pop.drop(columns = ['adm0_name', 'adm1_name','adm2_name', 'adm0_code', 'adm1_code','disp_area','geometry','status'])
dff_pop_rename = dff_pop.rename(columns={'sum':'Pop_adm2_WPOP2015'})

dff_jrc = df_jrc.drop(columns = ['adm0_name', 'adm1_name','adm2_name', 'adm0_code', 'adm1_code','disp_area','geometry','status'])
dff_jrc_rename = dff_jrc.rename(columns={'sum':'JRC_Flood_Haz'})

dff_exp = df_exp.drop(columns = ['adm0_name', 'adm1_name','adm2_name', 'adm0_code', 'adm1_code','disp_area','geometry','status'])
dff_exp_rename = dff_exp.rename(columns={'sum':'Exp_Pop_1m'})

#Change PXR RFIDF
dff_pxr = df_pxr.drop(columns = ['adm0_name', 'adm1_name','adm2_name', 'adm0_code', 'adm1_code','disp_area','geometry','status'])
dff_pxr_rename = dff_pxr.rename(columns={'mode':'MaxTP_3D'})

dff_pop1km = df_pop1km.drop(columns = ['adm0_name', 'adm1_name','adm2_name', 'adm0_code', 'adm1_code','disp_area','geometry','status'])
dff_pop1km_rename = dff_pop1km.rename(columns={'sum':'Pop_adm2_CSN2019'})

#dff_gfs = df_gfs.drop(columns = ['adm0_name', 'adm1_name','adm2_name', 'adm0_code', 'adm1_code','disp_area','geometry','status'])
#df_gfs_rename = dff_gfs.rename(columns={'mode':'GFS_3D'})

#Merge all DF
pre_merged = dff_p_rename.merge(dff_sm2_rename, on='adm2_code').merge(dff_pop_rename,on='adm2_code').merge(dff_jrc_rename,on='adm2_code').merge(dff_exp_rename,on='adm2_code').merge(dff_pxr_rename,on='adm2_code').merge(dff_pop1km_rename,on='adm2_code')#.merge(df_gfs_rename,on='adm2_code')
pre_merged['MaxTP_3D'] = pre_merged['MaxTP_3D'] * 1000
pre_merged = pre_merged.round(3)



### For local testing not connected to WFP domain (should have data locally)
#ecmwf_F_3D = '/Users/michael/GEO/adam-floods-alert/processing/process_ecmwf/0830090200TP.tiff'
#ecmwf_F_5D = '/Users/michael/GEO/adam-floods-alert//processing/process_ecmwf/0830090400TP.tiff'

### Uncomment for NRT
ecmwf_F_3D = process_ecmwf + ecmwf_3D_fc
ecmwf_F_5D = process_ecmwf + ecmwf_5D_fc
ecmwf_F_3D_stats = zonal_stats(pre_merged, ecmwf_F_3D, prefix='ecmwf3D_',geojson_out=True)
ecmwf_F_5D_stats = zonal_stats(pre_merged, ecmwf_F_5D, prefix='ecmwf5D_',geojson_out=True)
ecmwf_F_3D_gdf = GeoDataFrame.from_features(ecmwf_F_3D_stats)
ecmwf_F_5D_gdf = GeoDataFrame.from_features(ecmwf_F_5D_stats)

ecmwf_F_3D_gdf['ecmwf3D_mean'] = ecmwf_F_3D_gdf['ecmwf3D_mean'] * 1000 # converts TP values from meters to mm
ecmwf_F_5D_gdf['ecmwf5D_mean'] = ecmwf_F_5D_gdf['ecmwf5D_mean'] * 1000 # converts TP values from meters to mm


dfNew = ecmwf_F_5D_gdf.merge(ecmwf_F_3D_gdf, left_index=True, right_index=True,
                 how='outer', suffixes=('', '_y'))
dfNew.drop(list(dfNew.filter(regex='_y$')), axis=1, inplace=True)
to_drop = ['ecmwf5D_count','ecmwf5D_max','ecmwf5D_min','ecmwf3D_count','ecmwf3D_max','ecmwf3D_min']
gbl_alert = dfNew.drop(columns=to_drop).round(2)
gbl_alert['uuid'] = [uuid.uuid4() for _ in range(len(gbl_alert.index))]
print('Global Alerts created...')
#gbl_alert = gbl_alert.loc[(gbl_alert['3D_acc_precip'] >= 50)]


# In[ ]:


### Before 19 June 2019
#map_me = alert.loc[(alert['status'] == 'red') | (alert['P100'] == True) | (alert['ecmwf3D_mean'] > 50)]

### Deployed version 19 June 2019
#map_me = alert.loc[((alert['status'] == 'green') & (alert['P60'] == True)) | ((alert['status'] == 'red') & (alert['ecmwf3D_mean'] > alert['Mean_Prec']*0.6))] 

#map_me = alert.loc[(alert['3D_acc_precip'] > alert['MaxInt_3D_x']) | alert['P80'] == True]
#map_me.head()


# In[ ]:


# amsr2 = get_amsr2()
# amsr2_stats = zonal_stats(alert, amsr2, prefix='amsr2',geojson_out=True)
# zonal_ppt_gdf = GeoDataFrame.from_features(amsr2_stats)
# zonal_ppt_gdf = zonal_ppt_gdf.loc[(zonal_ppt_gdf['amsr2max'] > 10000)]
# zonal_ppt_gdf.head()


# In[ ]:


def func1(x):
    if x >= 50000:
        return "red"
    elif x <= 10000:
        return "green"
    else:
        return "orange"
    
def func2(x):
    if x >= 500000:
        return "red"
    elif x <= 10000:
        return "green"
    else:
        return "orange"

# def func3(df):
#     if df['nowcast_mod'] == True & df['pop_status'] == 'red':
#         return "green"
#     elif df['nowcast_high'] == True & df['pop_status'] == 'red':
#         return "orange"
#     elif df['forecast_high'] == True & df['pop_status'] == 'red':
#         return "red"
#     else:
#         pass
# gbl_alert['status'] = gbl_alert.apply(func3)

#ABDELS COMMENTS
#map_me = alert.loc[((alert['3D_acc_precip'] + alert['ecmwf3D_mean']) > alert['Mean_Prec'])] 

#Conditions goes here. Needs more work here
#gbl_alert = gbl_alert.loc[gbl_alert['3D_acc_precip'] >=50]
gbl_alert['pop_status'] = gbl_alert['Pop_adm2_CSN2019'].apply(func2)
gbl_alert['exp_status'] = gbl_alert['Exp_Pop_1m'].apply(func1)

#gbl_alert['P60'] = (gbl_alert['3D_acc_precip'] > gbl_alert['Mean_Prec']*0.60)
#gbl_alert['P80'] = (gbl_alert['3D_acc_precip'] > gbl_alert['Mean_Prec']*0.80)
#gbl_alert['P100'] = (gbl_alert['3D_acc_precip'] > gbl_alert['Mean_Prec'])
#alert['ecmwf'] = alert['ECMWF_Mean_Prec'].apply(lambda x: x*1000) #check with abdel TP values

#alert = alert[['adm0_name','adm1_name','adm2_name','adm2_code','status','MaxTP_3D','Mean_Prec','3D_acc_precip','ecmwf3D_mean','ecmwf5D_mean','Est_pop_adm2','JRC_Flood_Haz','Exp_Pop_1m','P60','P80','P100']]

#Condition 1
#gbl_alert['nowcast_low'] = (gbl_alert['3D_acc_precip'] > gbl_alert['MaxTP_3D']).astype(int)
#Condition 2
#gbl_alert['nowcast_mod'] = ((gbl_alert['3D_acc_precip'] + gbl_alert['ecmwf3D_mean']) > gbl_alert['MaxTP_3D']).astype(int)
gbl_alert['nowcast_mod'] = (gbl_alert['3D_acc_precip'] > gbl_alert['MaxTP_3D']).astype(int)
#Condition 3
gbl_alert['nowcast_high'] = (gbl_alert['3D_acc_precip'] > gbl_alert['Mean_Prec']).astype(int)
#Condition 3
gbl_alert['forecast_high'] = ((gbl_alert['3D_acc_precip'] + gbl_alert['ecmwf3D_mean'] > (gbl_alert['Mean_Prec'] + gbl_alert['Mean_Prec']*0.60)) | (gbl_alert['3D_acc_precip'] + gbl_alert['ecmwf5D_mean'] > (gbl_alert['Mean_Prec'] + gbl_alert['Mean_Prec']*0.60))).astype(int)
#gbl_alert = gbl_alert.loc[(gbl_alert['nowcast_mod'] == True & gbl_alert['pop_status'] == 'red') | (gbl_alert['nowcast_high'] == True) | (gbl_alert['forecast_high'] == True)]

#People Living only
# gbl_alert = gbl_alert.loc[((gbl_alert['nowcast_mod'] == True) & (gbl_alert['pop_status'] == 'red') | 
#                            (gbl_alert['nowcast_high'] == True) & (gbl_alert['pop_status'] == 'red') | 
#                            (gbl_alert['forecast_high'] == True) & (gbl_alert['pop_status'] == 'red'))]

##Exposure
gbl_alert.loc[(gbl_alert['3D_acc_precip'] >= 50), 'status'] = 'white'
gbl_alert.loc[(gbl_alert['nowcast_mod'] == True) & (gbl_alert['exp_status'] == 'red'),'status']  = 'yellow'
gbl_alert.loc[(gbl_alert['nowcast_high'] == True) & (gbl_alert['exp_status'] == 'red'), 'status'] = 'orange'
gbl_alert.loc[(gbl_alert['forecast_high'] == True) & (gbl_alert['exp_status'] == 'red'), 'status'] = 'red'
final_df = gbl_alert.loc[(gbl_alert['status'] == 'yellow') | (gbl_alert['status'] == 'orange')| (gbl_alert['status'] == 'red')| (gbl_alert['status'] == 'white')]


# In[ ]:


summary = final_df.groupby('status')['status'].value_counts()
sss = summary.values.tolist()

try:
    w = ''
    y = ''
    o = ''
    r = ''
    w = summary[0]
    y = summary[1]
    o = summary[2]
    r = summary[3]
except Exception:
    pass 

### GIS Officers in the area    
# recip = pd.read_excel(root + '/wfp' + '/wfp_gis.xlsx', index_col=0, sheet_name='GIS_officer')
# input_countries = map_me['adm0_name'].to_list()
# countries = {}
# for country in pycountry.countries:
#     countries[country.name] = country.alpha_3
# c = [countries.get(c, 'Unknown code') for c in input_countries]
# codes = list(set(c))
# codes_df = pd.DataFrame(codes,columns=['ISO3'])

# hoomans = []

# for i in codes:
#     gis_officers = recip.loc[recip['ISO3']==i]    
#     gis_officers_ls = gis_officers.values.tolist()
#     #gis_officer_dict = dict(gis_officers_ls)
#     hoomans.append(gis_officers_ls)   
# list2 = [x for x in hoomans if x != []]
# l = [i.strip('[]') if type(i) == str else str(i) for i in list2]
# a = "\n".join(l)
# a


geo = pre_merged.merge(final_df, on='adm2_code')
geo.drop(list(geo.filter(regex='_y$')), axis=1, inplace=True)
geo = geo.rename(columns={'geometry_x':'geometry'})
#geo.tail()


# ### Save WHITE alerts to DB

# In[ ]:


config.read(os.path.join(root, "config.txt"))
dbserver_out = config.get('db_out','dbserver')
dbname_out = config.get('db_out','dbname')
dbuser_out = config.get('db_out','dbuser')
dbpword_out = config.get('db_out','dbpword')
dbport_out = config.get('db_out','dbport')
engine = create_engine('postgresql+psycopg2://' + dbuser_out + ':' + dbpword_out + '@' + dbserver_out + ':' + dbport_out + '/' + dbname_out)

try:
    df = pd.DataFrame(final_df)
    df.head(0).to_sql(whitelist_data, engine, if_exists='replace',index=False) #truncates the table
    conn = engine.raw_connection()
    cur = conn.cursor()
    output = io.StringIO()
    df.to_csv(output, sep='\t', header=False, index=False)
    output.seek(0)
    contents = output.getvalue()
    cur.copy_from(output, whitelist_data, null="") # null values become ''
    conn.commit()

    sqlstr = "ALTER TABLE {table_name} ALTER COLUMN geometry TYPE geometry;".format(** {
                'table_name': whitelist_data})

    cur.execute(sqlstr)
    conn.commit()
    cur.close()
    print('White Alerts Sent to DB')
except Exception:
    pass 


# ### Save RED alerts to flood.db

# In[ ]:


sensing_date=('NOW-11DAYS',date_cut) #XDAYS XMONTH XWEEK


# In[ ]:


config.read(os.path.join(root, "config.txt"))
dbserver_al = config.get('db_alert','dbserver_a')
dbname_al = config.get('db_alert','dbname_a')
dbuser_al = config.get('db_alert','dbuser_a')
dbpword_al = config.get('db_alert','dbpword_a')
dbport_al = config.get('db_alert','dbport_a')
alert_engine = create_engine('postgresql+psycopg2://' + dbuser_al + ':' + dbpword_al + '@' + dbserver_al + ':' + dbport_al + '/' + dbname_al)

try:
    if final_df.empty == False: 
        df_alert = pd.DataFrame(final_df)
        df_alert.head(0).to_sql(alert_data, alert_engine, if_exists='replace',index=False) #truncates the table
        conn2 = alert_engine.raw_connection()
        cur2 = conn2.cursor()
        output2 = io.StringIO()
        df_alert.to_csv(output2, sep='\t', header=False, index=False)
        output2.seek(0)
        contents2 = output2.getvalue()
        cur2.copy_from(output2, alert_data, null="") # null values become ''
        conn2.commit()

        sqlstr2 = "ALTER TABLE {table_name} ALTER COLUMN geometry TYPE geometry;".format(** {
                    'table_name': alert_data})

        cur2.execute(sqlstr2)
        conn2.commit()
        cur2.close()
        print('Orange Alerts Sent to DB')
        
    else:
        print('No Alerts for today. Check White alerts if necessary!')
except Exception:
    pass 

try:
    if final_df.empty == True:
        pass
    else:
        bbox = final_df.envelope
        bbox_gdf = gpd.GeoDataFrame(gpd.GeoSeries(bbox), columns=['geometry'])
        c = pd.concat([bbox_gdf,final_df],axis=1)
        x = bbox.to_json()
        m = folium.Map(tiles='cartodbpositron') #cartodbpositron #stamentoner #cartodbdark_matter
        folium.GeoJson(x).add_to(m)
        m.fit_bounds(m.get_bounds())
        map_view_name = 'Map_view_' + date_cut + '.html'
        view = m.save(map_view_path + map_view_name)
        alert_geo = pd.DataFrame(final_df)
        alert_content = alert_geo.drop(columns = ['geometry'])
        att_alert = alert_content.to_csv(alert_csv_path + alert_data_name)
        
        recip = pd.read_excel(root + '/wfp' + '/wfp_gis.xlsx', index_col=0, sheet_name='GIS_officer')
        input_countries = final_df['adm0_name'].to_list()
        countries = {}
        for country in pycountry.countries:
            countries[country.name] = country.alpha_3
        c = [countries.get(c, 'Unknown code') for c in input_countries]
        codes = list(set(c))
        codes_df = pd.DataFrame(codes,columns=['ISO3'])

        hoomans = []
        summary = final_df.groupby('status')['status'].value_counts()
        sss = summary.values.tolist()
        
        #grouped = final_df.groupby('status')
        #grouped = final_df.groupby('nowcast_high')#['nowcast_high'].value_counts()
        #s_nc = grouped.values.tolist()
        #summary_fc = gbl_alert.groupby('forecast_high')#['forecast_high'].value_counts()
        #s_fc = summary_fc.values.tolist()
        #for name,group in grouped:

        for index, row in final_df.iterrows():
            if row.status == 'orange':
                get_s1_info(index)
            if row.status == 'red':
                get_s1_info(index)
        
        for i in codes:
            gis_officers = recip.loc[recip['ISO3']==i]    
            gis_officers_ls = gis_officers.values.tolist()
            #gis_officer_dict = dict(gis_officers_ls)
            hoomans.append(gis_officers_ls)   
        list2 = [x for x in hoomans if x != []]
        l = [i.strip('[]') if type(i) == str else str(i) for i in list2]
        a = "\n".join(l)
        job()
                
except Exception:
    pass

print("The script took {0} seconds to complete...".format(time.time() - execTime))


# ### End of script..


# ## Working scratch do not delete!

# In[ ]:


### SCRATCH
# conditions = [(df['column_1'] > 5),
#               (df['column_1'] <= 5) & (df['column_1'] > 0),
#               (df['column_1'] == 0)]

# choices = ['high','low','null']

# df['column_2'] = np.select(conditions, choices, default='null')

# summary_nc = gbl_alert.groupby('nowcast_high')['nowcast_high'].value_counts()
# s_nc = summary_nc.values.tolist()
# summary_fc = gbl_alert.groupby('forecast_high')['forecast_high'].value_counts()
# s_fc = summary_fc.values.tolist()
# try:
#     nc = ''
#     fc = ''
#     nc = s_nc[1]
#     fc = s_fc[1]
# except Exception:
#     pass


# In[ ]:


### This works OK. 
# for name,group in grouped_nc:
#     for i in codes:
#         gis_officers = recip.loc[recip['ISO3']==i]    
#         gis_officers_ls = gis_officers.values.tolist()
#         #gis_officer_dict = dict(gis_officers_ls)
#         hoomans.append(gis_officers_ls)   
#     list2 = [x for x in hoomans if x != []]
#     l = [i.strip('[]') if type(i) == str else str(i) for i in list2]
#     a = "\n".join(l)
#     job()


# In[ ]:


# i = 0
# while i == 0:
#     if name == 'red':
#         i = 0
#         print('i am red')
#     elif name == 'orange':
#         i = 0
#         print('i am orange')
#     else:
#         i += 1
#         print('i am green')
                


# In[ ]:


# for name,group in grouped:
#     if name == 'green':
#         for i in codes:
#             gis_officers = recip.loc[recip['ISO3']==i]    
#             gis_officers_ls = gis_officers.values.tolist()
#             #gis_officer_dict = dict(gis_officers_ls)
#             hoomans.append(gis_officers_ls)   
#         list2 = [x for x in hoomans if x != []]
#         l = [i.strip('[]') if type(i) == str else str(i) for i in list2]
#         a = "\n".join(l)
#         job()
#     else:
#         pass


# ### Email works via wfp domain

# In[ ]:


# alert_geo = pd.DataFrame(map_me)
# alert_content = alert_geo.drop(columns = ['geometry'])
# att_alert = alert_content.to_csv(alert_csv_path + alert_data_name)

# #map_view_name = 'Map_view_' + date_cut + '.html'

# job()

# # try:
# #     job()
# # except Exception:
# #     pass 

# execTime = time.time()
# print("The script took {0} seconds to complete...".format(time.time() - execTime))


# ### Request HTA Map tiles

# ### Get HTA Tiles of the alerted Area (wld_grid_100ka1_wfp) @ 43

# In[ ]:


# config.read(os.path.join(root, "config.txt"))
# dbserver_hta = config.get('hta','dbserver_hta')
# dbname_hta = config.get('hta','dbname_hta')
# dbuser_hta = config.get('hta','dbuser_hta')
# dbpword_hta = config.get('hta','dbpword_hta')
# dbport_hta = config.get('hta','dbport_hta')
# hta_engine = create_engine('postgresql+psycopg2://' + dbuser_hta + ':' + dbpword_hta + '@' + dbserver_hta + ':' + dbport_hta + '/' + dbname_hta)
# conn3 = hta_engine.raw_connection()
# cur3 = conn3.cursor()
# sql = """  SELECT adm0_name,iso3,map_code,map_done,url,shape as geometry
#          FROM wfp_pub.wfp.wld_grid_100ka1_wfp where map_done = 'yes'  """
# remote_df = gpd.GeoDataFrame.from_postgis(sql,conn3,geom_col='geometry')
# conn3.close()

# # import folium
# # m = folium.Map(tiles='cartodbdark_matter') #cartodbpositron #stamentoner #cartodbdark_matter
# # folium.GeoJson(remote_df).add_to(m)
# # m.fit_bounds(m.get_bounds())
# # m


# In[ ]:


# conn3 = hta_engine.raw_connection()
# cur3 = conn3.cursor()
# hta_sql = """  SELECT adm0_name,iso3,map_code,map_done,url,shape as geometry
#           FROM wfp_pub.wfp.wld_grid_100ka1_wfp where map_done != 'yes'  """
# remote_hta_gdf = gpd.GeoDataFrame.from_postgis(hta_sql,conn3,geom_col='geometry')
# request_hta = gpd.overlay(remote_hta_gdf, alert, how='intersection')
# #conn3.close()

# config.read(os.path.join(root, "config.txt"))
# dbserver_hta = config.get('hta','dbserver_hta')
# dbname_hta = config.get('hta','dbname_hta')
# dbuser_hta = config.get('hta','dbuser_hta')
# dbpword_hta = config.get('hta','dbpword_hta')
# dbport_hta = config.get('hta','dbport_hta')
# hta_engine = create_engine('postgresql+psycopg2://' + dbuser_hta + ':' + dbpword_hta + '@' + dbserver_hta + ':' + dbport_hta + '/' + dbname_hta)
# conn3 = hta_engine.raw_connection()
# cur3 = conn3.cursor()
# sql = """  SELECT adm0_name,iso3,map_code,map_done,url,shape as geometry
#          FROM wfp_pub.wfp.wld_grid_100ka1_wfp where map_done = 'yes'  """
# remote_df = gpd.GeoDataFrame.from_postgis(sql,conn3,geom_col='geometry')
# conn3.close()


# In[ ]:


# ## Optional Plot ECMWF raster forecast
# raster = rasterio.open(ecmwf_F_3D)
# array = raster.read()
# from IPython.display import Image
# %matplotlib notebook
# %matplotlib inline
# show(array)


# ### ECMWF Local Compute Processing - Work In Progress (xarray/dask)

# In[ ]:


# out_file = '/Users/michael/GEO/process_cdsapi/'

# c = cdsapi.Client()

# c.retrieve(
#     'reanalysis-era5-single-levels-monthly-means',
#     {
#         'format':'netcdf',
#         'product_type':'monthly_averaged_reanalysis',
#         'variable':'total_precipitation',
#         'year':[
#             '1979','1980','1981',
#             '1982','1983','1984',
#             '1985','1986','1987',
#             '1988','1989','1990',
#             '1991','1992','1993',
#             '1994','1995','1996',
#             '1997','1998','1999',
#             '2000','2001','2002',
#             '2003','2004','2005',
#             '2006','2007','2008',
#             '2009','2010','2011',
#             '2012','2013','2014',
#             '2015','2016','2017',
#             '2018'
#         ],
#         'month':[
#             '01','02','03',
#             '04','05','06',
#             '07','08','09',
#             '10','11','12'
#         ],
#         'time':'00:00'
#     },
#     outfile + 'monthly_average_reanalysis.nc')

# tot_precip = '/Users/michael/GEO/process_cdsapi/total_precip.tif'
# cdsapi_stats = zonal_stats(alert, tot_precip, prefix='cds_x',geojson_out=True)
# zonal_ppt_gdf = GeoDataFrame.from_features(cdsapi_stats)

### Reference: https://github.com/jwagemann/seasonal_forecasts/blob/master/Workflow_seasonal_fc_processing.ipynb
### Reference: https://annefou.github.io/metos_python/07-LargeFiles/


# In[ ]:


#WorldClim V1
#worldclim = ee.ImageCollection('WORLDCLIM/V1/MONTHLY')
#precip = worldclim.select('prec')
#collectionList = precip.toList(precip.size())


# In[ ]:


#import hxl 
#source = hxl.data('http://wfp.org/dataset.xlsx')
#from hdx.location.country import Country

