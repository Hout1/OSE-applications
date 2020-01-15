# In[2]:


import feedparser
import pandas as pd
#import matplotlib.pyplot as plt
import geopandas as gpd
from geopandas import GeoDataFrame
from shapely.geometry import Point
import fiona
import folium
import os
from sqlalchemy import create_engine
import psycopg2 
from io import BytesIO as StringIO
from datetime import datetime, timedelta
import ssl
if hasattr(ssl, '_create_unverified_context'):
    ssl._create_default_https_context = ssl._create_unverified_context
    
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


# In[3]:


#Get RSS
feed = 'https://emergency.copernicus.eu/mapping/list-of-components/EMSR348/feed'
gfds = 'http://www.gdacs.org/flooddetection/data.aspx?from=20190323&to=20190323&type=rss&alertlevel=red&datatype=1DAYS'
EMS = 'http://emergency.copernicus.eu/mapping/activations-rapid/feed' #Filter flood #filter disastertype (flood)
gdacs = 'http://www.gdacs.org/xml/rss.xml'
eo = 'https://earthobservatory.nasa.gov/feeds/earth-observatory.rss'

GDACS_RSS = 'http://www.gdacs.org/xml/rss.xml'
RF_RSS ='https://reliefweb.int/disasters/rss.xml'
EMS_RSS = 'http://emergency.copernicus.eu/mapping/activations-rapid/feed'
GFDS_RSS = 'http://www.gdacs.org/flooddetection/data.aspx?from=20190323&to=20190323&type=rss&alertlevel=red&datatype=1DAYS'


# ## GDACS

# In[4]:


GDACS_URL = [GDACS_RSS]

gdacs_feeds = []
gdacs_iso3 = []
gdacs_country = []
gdacs_etype = []
gdacs_epid = []
gdacs_link = []
gdacs_guid = []
gdacs_pubdate = []
gdacs_bbox = []
gdacs_coords = []
gdacs_alevel = []
gdacs_ascore = []

wfpiso3 = ['TUR','CHL', 'IDN', 'PHL', 'SDN', 'SSD', 'DRC', 'BGD', 'CAF', 
           'BWA', 'TCD', 'KHM', 'IND', 'MYS','ZWE','MDG','MOZ']

def get_gdacs():
    for url in GDACS_URL:
        gdacs_feeds.append(feedparser.parse(url))

    for feed in gdacs_feeds:
        for post in feed.entries:
            gdacs_iso3.append(post.gdacs_iso3)
            gdacs_country.append(post.gdacs_country)
            gdacs_etype.append(post.gdacs_eventtype)
            gdacs_epid.append(post.gdacs_episodeid)
            gdacs_link.append(post.link)
            gdacs_guid.append(post.guid)
            gdacs_pubdate.append(post.published)
            gdacs_bbox.append(post.gdacs_bbox)
            gdacs_coords.append(post.where)
            gdacs_alevel.append(post.gdacs_episodealertlevel)
            gdacs_ascore.append(post.gdacs_alertscore)
            
        guiddf = pd.DataFrame(gdacs_guid)       #0
        pubdatedf = pd.DataFrame(gdacs_pubdate)
        isodf = pd.DataFrame(gdacs_iso3)        #2   
        countrydf = pd.DataFrame(gdacs_country) #3
        etypedf = pd.DataFrame(gdacs_etype)     #4
        epidf = pd.DataFrame(gdacs_epid)        #5
        linkdf = pd.DataFrame(gdacs_link)        #6
        bboxdf = pd.DataFrame(gdacs_bbox)       #7
        coordf = pd.DataFrame(gdacs_coords)     #8 #9
        aleveldf = pd.DataFrame(gdacs_alevel)   #10
        ascoredf = pd.DataFrame(gdacs_ascore)   #11
        
        frames = [guiddf,etypedf,pubdatedf, coordf,linkdf,] #, isodf, countrydf, etypedf, epidf, evidf, bboxdf,coordf, aleveldf, ascoredf]
        df = pd.concat(frames, axis=1, ignore_index = True, names=[frames])
        #return df[df[2].isin(wfpiso3)] #!IMPORTANT!#
        
        tc = df.loc[df[1] == 'TC']
        fl = df.loc[df[1] == 'FL']
        
        x = pd.concat([fl,tc],ignore_index=True, sort=False)
        return x


# In[5]:


EMS_URL = [EMS_RSS]

ems_guid = []
ems_title = []
ems_cat = []
ems_geo = []
ems_pub = []
ems_efeeds = []
ems_link = []

wfpiso3 = ['TUR','CHL', 'IDN', 'PHL', 'SDN', 'SSD', 'DRC', 'BGD', 'CAF', 
           'BWA', 'TCD', 'KHM', 'IND', 'MYS','ZWE','MDG','MOZ']

def get_ems():
    for url in EMS_URL:
        ems_efeeds.append(feedparser.parse(url))

    for feed in ems_efeeds:
        for post in feed.entries:
            ems_guid.append(post.guid)
            ems_title.append(post.title)
            ems_cat.append(post.category)
            ems_geo.append(post.where)
            ems_pub.append(post.published)
            ems_link.append(post.link)
            
        guidf = pd.DataFrame(ems_guid)   #0
        titledf = pd.DataFrame(ems_title)  #1
        catdf= pd.DataFrame(ems_cat)     #2
        geodf= pd.DataFrame(ems_geo)     #3
        pubdf= pd.DataFrame(ems_pub)   #4
        linkdf= pd.DataFrame(ems_link)
        
        frames = [guidf,catdf,pubdf,geodf,linkdf] #titledf, catdf, , 
        df = pd.concat(frames, axis=1, ignore_index = True, names=[frames])
        #return df[df[2].isin(wfpiso3)] #!IMPORTANT!#
        
        st = df.loc[df[1] == 'Storm']
        fl = df.loc[df[1] == 'Flood']
        
        x = pd.concat([fl,st],ignore_index=True, sort=False)
        return x


# In[6]:


RF_URL = [RF_RSS]
rw_feeds = []
rw_pub = []
rw_title = []
rw_category = []
rw_desc = []
rw_geo = []

def get_rw():
    for url in RF_URL:
        rw_feeds.append(feedparser.parse(url))

    for feed in rw_feeds:
        for post in feed.entries:
            rw_pub.append(post.published)
            rw_title.append(post.title)
            rw_category.append(post.category)
            rw_desc.append(post.description)
            
        titledf = pd.DataFrame(rw_title)       #0
        pubdatedf = pd.DataFrame(rw_pub) #1
        catdf = pd.DataFrame(rw_category)        #2   
        descdf = pd.DataFrame(rw_desc)
        geodf = pd.DataFrame(rw_geo)
        
        frames = [titledf,pubdatedf,geodf] #catdf, descdf
        df = pd.concat(frames, axis=1, ignore_index = True, names=[frames])
        #return df[df[2].isin(wfpiso3)] #!IMPORTANT!#
        return df


# In[7]:


GFDS_URL = [GFDS_RSS]

gfds_desc = []
gfds_title = []
gfds_country = []
gfds_geo = []
gfds_magnitude = []
gfds_pub = []
gfds_feeds = []

wfpiso3 = ['TUR','CHL', 'IDN', 'PHL', 'SDN', 'SSD', 'DRC', 'BGD', 'CAF', 
           'BWA', 'TCD', 'KHM', 'IND', 'MYS','ZWE','MDG','MOZ']

def get_gfds():
    for url in GFDS_URL:
        gfds_feeds.append(feedparser.parse(url))

    for feed in gfds_feeds:
        for post in feed.entries:
            gfds_desc.append(post.description)
            gfds_title.append(post.title)
            gfds_country.append(post.country)
            gfds_geo.append(post.where)
            gfds_magnitude.append(post.magnitude)
            #pub.append(post.record_date)
            
        descdf = pd.DataFrame(gfds_desc)   #0
        titledf = pd.DataFrame(gfds_title)  #1
        countrydf= pd.DataFrame(gfds_country)     #2
        geodf= pd.DataFrame(gfds_geo)     #3
        magdf = pd.DataFrame(gfds_magnitude)
        #pubdf= pd.DataFrame(pub)     #4
        
        frames = [descdf,geodf]#titledf,countrydf,, magdf
        df = pd.concat(frames, axis=1, ignore_index = True, names=[frames])
        #return df[df[2].isin(wfpiso3)] #!IMPORTANT!#
        #return df.where[df['category'] == 'FL']
        return df


# ## GET DATA

# In[8]:


gdacs_data = get_gdacs()
ems_data = get_ems()
gfds_data = get_gfds()
rw_data = get_rw()

#reliefweb_f = rw_data.loc[rw_data[2] == 'Flash Flood']
#gdacs_fl = gdacs_data.loc[gdacs_data[4] == 'FL']
#gdacs_tc = gdacs_data.loc[gdacs_data[4] == 'TC']
#ems_f = ems_data.loc[ems_data[2] == 'Flood']
#ems_s = ems_data.loc[ems_data[2] == 'Storm']


# In[9]:


stacked = gdacs_data.append(ems_data)

geometry = [Point(xy) for xy in zip(stacked[3])]
stacked = stacked.drop([3], axis=1)
crs = {'init': 'epsg:4326'}
stacked = GeoDataFrame(stacked, crs=crs, geometry=geometry)

stacked[2]= pd.to_datetime(stacked[2])

stacked.sort_values(by=[2],ascending=False )
#out = stacked.to_csv('/Users/michael/Downloads/test.csv')

# In[10]:


#filt = stacked[(stacked[2] >= '2019-04-01')]


# In[11]:


import os
from configparser import ConfigParser
root = os.path.dirname(os.path.realpath("__file__"))
adam_fl = os.path.dirname(root)

#data_dir = "/Users/michael/GEO/process_amsr2"
#print("I am dir:",adam_fl)
#print("I am root:", root)


config = ConfigParser()
config.read(os.path.join(adam_fl, "config.txt"))

from sqlalchemy import create_engine
import psycopg2 
config.read(os.path.join(root, "config.txt"))
dbserver_out = config.get('db_out','dbserver')
dbname_out = config.get('db_out','dbname')
dbuser_out = config.get('db_out','dbuser')
dbpword_out = config.get('db_out','dbpword')
dbport_out = config.get('db_out','dbport')
engine = create_engine('postgresql+psycopg2://' + dbuser_out + ':' + dbpword_out + '@' + dbserver_out + ':' + dbport_out + '/' + dbname_out)

# In[14]:


now = datetime.today()
d = str(now)
date = d[0:10]
date_cut = date.replace('-', '')
data = 'georss' + date_cut

georss = pd.DataFrame(stacked)

georss.head(0).to_sql(data, engine, if_exists='replace',index=False) #truncates the table
conn = engine.raw_connection()
cur = conn.cursor()
output = StringIO()
georss.to_csv(output, sep='\t', header=False, index=False)
output.seek(0)
contents = output.getvalue()
cur.copy_from(output, data, null="") # null values become ''
conn.commit()

sqlstr = "ALTER TABLE {table_name} ALTER COLUMN geometry TYPE geometry;".format(** {
            'table_name': data})

cur.execute(sqlstr)

conn.commit()
cur.close()

print("Events saved to DB!") 


# In[13]:


#aoi = gdacsjson.loc[gdacsjson['alertlevelepisode'] == 'Green']
        
# m = folium.Map(tiles='cartodbdark_matter') #cartodbpositron #stamentoner #cartodbdark_matter
# folium.GeoJson(stacked).add_to(m)
# m.fit_bounds(m.get_bounds())
# m


# In[14]:


# root = "http://www.gdacs.org/datareport/resources/"
# eventTC = "TC/"
# eventDR = "DR/"
# eventEQ = "EQ/"
# eventID = "1000533"
# objID = "/geojson_"
# objNum = "_1"
# objForm = ".geojson"

# tceventtest = root + eventTC + eventID + objID + eventID + objNum + objForm

# print(tceventtest)

