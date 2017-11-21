import ESMF
from esmf_utils import grid_create, grid_create_periodic
import numpy as np
import numpy.ma as ma
import math as math
import scipy.io.netcdf as nio
import scipy.io as io
from scipy.interpolate import griddata
from scipy import interpolate
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from shapely.geometry.polygon import LinearRing
import os
from mpl_toolkits.basemap import Basemap
import sys
import dist
from scipy import stats
import calendar
# Some locally installed stuff
#import h5py
import gsw
import seawater as sw
from netCDF4 import Dataset
import math

def read_sections(filename):
  '''Read the section indeces defined on NorESM grid. Note that these indeces are defined for Fortran 1 based system so one needs to substract one to get the correct python 0 based indeces'''
  f=open(filename)
  data=f.readlines()
  f.close()
  datanames=[]
  c=0
  for j in range(len(data)+1):
    if j==len(data) or data[j][0:4]=='Name':
     #save the old data
      if (j!=0):
        exec(datanames[c]+'=np.ones((len(dummy),len(dummy[0])))')
        for k in range(len(dummy)):
          exec(datanames[c]+'[k,:]=dummy[k]')
        c=c+1
      if (j!=len(data)):
        #intialize the new variables
        exec('datanames.append("'+data[j][6:-1].strip()+'")')
        dummy=[]
    else:
      dummy.append(np.array(data[j].strip().split(),dtype=int))
  #
  #Put everything into a dictionary and return that
  output={}
  for dataname in datanames:
    exec('output.update({"'+dataname+'": '+dataname+'})')
  #
  return output

def read_section_data(fpath,model,ens,files):
  """Read section transports from CMIP5 models. MFO is the standard output."""
  filenames=[]
  for fil in files:
      filenames.append(fpath+model+'/'+ens+'/'+fil)
  #Read the passage names
  f=nio.netcdf_file(filenames[0])
  passages=f.variables['passage']
  dummy=[]
  #print model
  for passage in passages:
    dummy.append(passage.tostring().replace('\0','').strip())
  #
  #print dummy
  f.close()
  #read the files and construct time and data arrays
  #this is needed because some of the models have more than one file per ensemble
  time_size=[] #length of each file (how many time steps)
  time_size.append(0)
  for c in range(len(filenames)):
    exec('f'+str(c)+'=nio.netcdf_file(filenames[c])')
    exec('time_size.append(len(f'+str(c)+'.variables["time"][:]))')
  #After checking how many files there are create arrays for the mfo and time
  mfo=np.ones((sum(time_size),len(dummy)))*np.nan
  timeaxis=np.ones(sum(time_size))*np.nan
  time_size=np.cumsum(time_size)
  #read the data
  for c in range(len(filenames)):
    exec('mfo[time_size[c]:time_size[c+1],:]=f'+str(c)+'.variables["mfo"][:,:]')
    exec('timeaxis[time_size[c]:time_size[c+1]]=f'+str(c)+'.variables["time"][:].squeeze()')
    exec('f'+str(c)+'.close()')
  #Create the output dictionary
  output={}
  #output.update({"passage_names": dummy})
  output.update({"time": timeaxis})
  c=0
  for passage in dummy:
    exec('output.update({"'+passage+'": mfo[:,c]})')
    c=c+1
  #
  return output, dummy

def section_data():
  #Read the mfo data, start by checking which models and ensembles are available
  fpath='/Data/skd/share/ModData1/CMIP5/ocean/rcp85/mfo/mon/'
  models=os.listdir(fpath)
  ens_num=np.zeros(len(models))
  c=0
  for model in models:
    exec('ens_num[c]=len(os.listdir("'+fpath+model+'/"))')
    c=c+1
  #use the read_section_data for the actual reading
  data={}
  p=0
  for model in models:
   for j in range(int(ens_num[p])):
    exec('mfo_'+model.replace('-','_')+'_r'+str(j+1)+', passage_names=read_section_data(fpath,model,"r'+str(j+1)+'i1p1",os.listdir("'+fpath+model+'/r'+str(j+1)+'i1p1/"))')
    exec('data.update({"mfo_'+model.replace('-','_')+'_r'+str(j+1)+'":mfo_'+model.replace('-','_')+'_r'+str(j+1)+'})')
   p=p+1
  #
  data.update({"passage_names": passage_names})
  data.update({"models": models})
  #
  return data 
 
def PAGO_sections(sections,modelname,fname):
  #READ IN SECTIONS DEFINED IN PAGO AND PUT THEM INTO A NICE FORMAT
  #sections=['brw','brn','be1','be2','fst','lan','nrs','BER','baf','dso','ifo','fso']
  #note the corrections, somehow PAGO assumes that v to be on the north face so have to change v or u points up or down depending on the direction of the turn
  model_sections=io.loadmat("/Data/skd/users/anu074/CMIP_PAGO/"+modelname+fname)
  veci=[]; vecj=[]; udir=[]; vdir=[];
  for s, section in enumerate(sections):
      j=np.where(model_sections["MODEL_sections"]["name"][0][:]==section)[0][0]
      dir1=model_sections["MODEL_sections"][0][j][3].copy()
      dir2=model_sections["MODEL_sections"][0][j][4][0].copy()
      jj=[]; ii=[]; carryover1=False; carryover2=False
      ii.append(int(model_sections["MODEL_sections"][0][j][1][0][0].copy()))
      jj.append(int(model_sections["MODEL_sections"][0][j][2][0][0].copy()))
      for k in range(1,len(model_sections["MODEL_sections"][0][j][1][0])):
        i1=int(model_sections["MODEL_sections"][0][j][1][0][k].copy())
        j1=int(model_sections["MODEL_sections"][0][j][2][0][k].copy())
        if (dir1[k]!=dir1[k-1] and dir1[k-1]=='W' and ((int(dir2[k-1])==1 and int(dir2[k])==-1) or (int(dir2[k-1])==-1 and int(dir2[k])==1)) ) or (j1<jj[k-1] and carryover1): #from west face to north face
          ii.append(i1)
          jj.append(j1+1) #v should be read from north face
          carryover1=True
          if carryover2:
             carryover2=False
        elif (dir1[k]!=dir1[k-1] and dir1[k-1]=='N' and ((int(dir2[k-1])==1 and int(dir2[k])==1) or (int(dir2[k-1])==-1 and int(dir2[k])==-1)) ) or (j1>jj[k-1] and carryover2): #from north face to west face
          ii.append(i1)
          jj.append(j1-1) #u should be read from the cell below the previous one
          carryover2=True
          if carryover1:
             carryover1=False
        else:
          ii.append(i1)
          jj.append(j1)
      veci.append(ii)#(model_sections["MODEL_sections"][0][j][1][0].copy())
      vecj.append(jj)#(model_sections["MODEL_sections"][0][j][2][0].copy())
      #dir1=model_sections["MODEL_sections"][0][j][3].copy()
      #dir2=model_sections["MODEL_sections"][0][j][4][0].copy()
      diru=dir2.copy(); dirv=dir2.copy()
      diru[list(np.where(dir1=='N')[0])]=0; udir.append(diru);
      dirv[list(np.where(dir1=='W')[0])]=0; vdir.append(dirv);
  model={"names": sections, "veci": veci, "vecj": vecj, "udir": udir, "vdir":vdir}
  
  return model

def fix_PAGO(s_v,s_u,data):
  """ Since the PAGO follows the vertices there can be a checkboard pattern in the velocity. This function is attempt to fix that by adding one velocity to the other"""
  if ma.sum(abs(s_v))>=ma.sum(abs(s_u)):
      ind1=ma.where(s_v)[0]
      ind2=ma.where(s_u)[0]
  else:
      ind1=ma.where(s_u)[0]
      ind2=ma.where(s_v)[0]
  if len(data.shape)>1:
    for i in ind2: #add to the 2 closest points - almost work, maybe should be added to all the points in the vicinity
      i1=ma.where(abs(ind1-i)==ma.min(abs(ind1-i)))[0]
      i1=ind1[i1]
      if len(i1)==1:
        data[:,i1[0]]=ma.masked_array(data[:,i1[0]].data+data[:,i].data,mask=data[:,i1[0]].mask)
      else:
        data[:,i1[0]]=ma.masked_array(data[:,i1[0]].data+.5*data[:,i].data,mask=data[:,i1[0]].mask)
        data[:,i1[1]]=ma.masked_array(data[:,i1[1]].data+.5*data[:,i].data,mask=data[:,i1[1]].mask)
    dataout=data[:,ind1]
  else:
    for i in ind2: #add to the 2 closest points - almost work, maybe should be added to all the points in the vicinity
      i1=ma.where(abs(ind1-i)==ma.min(abs(ind1-i)))[0]
      i1=ind1[i1]
      if len(i1)==1:
        data[i1[0]]=ma.masked_array(data.data[i1[0]]+data.data[i],mask=data.mask[i1[0]])
      else:
        data[i1[0]]=ma.masked_array(data.data[i1[0]]+.5*data.data[i],mask=data.mask[i1[0]])
        data[i1[1]]=ma.masked_array(data.data[i1[1]]+.5*data.data[i],mask=data.mask[i1[1]])
    dataout=data[ind1]
  
  return ind1,dataout

def find_bottompatch(x,y,maxd):
  #Note that hear we assume that we want to have a patch starting from
  #lower left corner and moving over with maximum depth maxd to the lower rigth corner
  #From there we move to the upper left corner with depths at y
  #This works for interpolated data
  b_patch=np.zeros((len(x)*2,2))
  for k in range(len(x)):
   b_patch[k,1]=maxd; b_patch[k+len(x),1]=y[len(x)-k-1];
   b_patch[k,0]=x[k]; b_patch[k+len(x),0]=x[len(x)-k-1];
  #
  #If we want follow the grid lines then we do the following
  b_patch2=np.zeros((len(x)*2+4,2))
  #bottom part
  b_patch2[0:2,1]=maxd; b_patch2[0,0]=min(x); b_patch2[1,0]=max(x);
  #surface part
  imid=int(len(x)/2.); imid2=imid*2
  b_patch2[2:imid2+2:2,0]=x[imid-1:-1][::-1]; #from left to middle
  b_patch2[3:imid2+2:2,0]=x[imid-1:-1][::-1]; #from left to middle
  b_patch2[imid2+3:-1:2,0]=x[:imid][::-1]; #from rigth to middle
  b_patch2[imid2+4:-1:2,0]=x[:imid][::-1]; #from rigth to middle
  b_patch2[imid2+2,0]=x[imid-1]; b_patch2[imid*2+2]=x[imid-1];
  b_patch2[2,0]=x[-1]; b_patch[-1,0]=x[0]
  #
  b_patch2[2,1]=y[-1]; b_patch2[-1,1]=y[0] #First surface points
  b_patch2[3:imid2+3:2,1]=y[imid-1:-1][::-1]; #from left to middle
  b_patch2[4:imid2+4:2,1]=y[imid-1:-1][::-1]; #from left to middle
  b_patch2[imid2+2:-2:2,1]=y[:imid][::-1]; #from rigth to middle
  b_patch2[imid2+3:-2:2,1]=y[:imid][::-1]; #from rigth to middle
  
  return b_patch, b_patch2


def point_inside_polygon(x,y,poly):
    """ As it says, find out if points x,y are inside polygon poly """
    #Determine if a point is inside a given polygon or not
    #Polygon is a list of (x,y) pairs.
    #Note that this does not work in lon,lat coordinates -> use basemap to convert lon,lat to some map projection coordinates
    #Adapted from http://www.ariel.com.au/a/python-point-int-poly.html
    n = len(poly)
    inside =False
    
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y
    
    return inside

#def insphpoly(lon,lat,lonv,latv,lon0,lat0):
#    """ INSPHPOLY True for points inside or on a polygonal region. INSIDE = INSPHPOLY(LON,LAT,LONV,LATV,LON0,LAT0) returns a matrix INSIDE the size of LON and LAT.  INSIDE(p,q) = 1 if the point (LON(p,q), LAT(p,q)) is either strictly inside or on the edge of the spherical polygonal region whose vertices are specified by the vectors LONV and LATV; otherwise INSIDE(p,q) = 0. All positions are stereographically projected onto a plane with (LON0, LAT0) in origo before polygon testing is done."""
#  #
#  #rad=pi/180;
#  #w=np.tan((90-lat0)*rad/2).*np.exp(np.i*lon0*rad); 
#  #z=np.tan((90-lat)*rad/2).*np.exp(np.i*lon*rad);
#  #z=(z-w)./(np.conj(w).*z+1);
#  #zv=np.tan((90-latv)*rad/2).*np.exp(np.i*lonv*rad);
#  #zv=np.(zv-w)./(conj(w).*zv+1);
#  #
#  #inside=inpolygon(real(z),imag(z),real(zv),imag(zv));
#  #return inside

def basin_coords(x,y,basin_name):
    """ Determine the basins based on the basin name and create the polygon """
    if basin_name=='nansen':
       #xp=[20.5, 20.5, 70, 100, 80, 50, 20.5]
       #yp=[82, 84, 85.5, 84, 83.5, 83.5, 82]
        xp=[350.,20.,20.,70.,110.,115.,125.,115.,70.,20.,350.,350.]
        yp=[81.,81.,81.,82.,80.,77.,77.,82.5,85.5,85.,85.,81.]
    elif basin_name=='amundsen':
       #xp=[330,160,110,135,125.5,100,80.5,0.5,330]
       #yp=[86.5, 89.5,88,84,83.5,87,88,86.5,86.5]
       xp=[350.,10.,70.,90.,115.,125.,126.,140.,140.,150.,140.,300.,300.,350.]
       yp=[85.5,85.5,85.5,85.,82.5,80.,78.,79.5,82.5,85.,87.5,87.5,85.5,85.5]
    elif basin_name=='makarov':
       #xp=[165, 165, 280, 280, 165]
       #yp=[86, 88, 88, 86, 86]
       xp=[300.,60.,145.,180.,180.,240.,285,300.]
       yp=[85.5,89.5,80.,79.,85.,86.,86.,85.5]
    elif basin_name=='canada':
       #xp=[220, 215, 210, 205, 208, 210, 225, 230, 230, 225, 225,220.5,220]
       #yp=[73, 72, 73, 74, 77, 81, 82, 82, 78, 77, 73,73,73]
       #xp=[200.,210.,220.,227.,227.,240.,250.,210.,185.,185.,205,205.,200.]
       #yp=[73.,71.5,71.5,72.5,77.5,80.,82.5,84.,84.,80.,80.,75.,73.]
       xp=[-103.84775688995,-132.426110480092,-146.32851987161,-157.129445685198,-165.289764266082,-174.948859566318,-169.908721619095,-153.99152252417,-153.654033773183,-155.991549118887,-157.264115189183,-158.99374858648,-156.200224800736,-149.147691533142,-146.758355108426,-144.575059890815,-142.121921555256,-140.683340702781,-138.653559249327,-137.204255767801,-134.56804180152,-129.946421652944,-127.799470427715,-127.884540272314,-127.978284839156,-127.025389763871,-120.96578199099,-110.556210660899,-109.124896536309,-108.367743820184,-105.94954447772,-103.84775688995]
       yp=[83.1681125452104,83.7662334620532,83.8235854688569,83.8506192123663,83.5387740470199,81.0961814080998,79.63568188243,77.4297453757112,76.9425777887249,76.1151804105261,74.9362263809805,73.5789329669126,72.8750577746576,71.4351508667292,71.0889175071879,71.0180792528114,70.7030506843364,70.6869725167624,70.4111375153443,70.4316926024079,70.6888027417022,71.8499446192718,72.9555979284964,73.7334233896998,74.5127217931511,75.3222606599783,77.8864195154148,80.0342652352252,80.3858803567825,81.84979009224,82.139759675712,83.1681125452104]
    elif basin_name=='eurasia':
       xp=[27.2474628708497,36.5447032085722,48.5214896162041,60.7001380803934,77.9434238049383,94.8206611280039,107.876624193844,114.76253701556,121.721644200227,134.412781688966,140.650432725167,141.964704250335,145.155248340791,149.246374638137,132.143273251575,149.76109018589,-79.6880622455017,-55.18591213471,-43.4591200137502,-22.0049915144465,-20.6758818057775,-17.5930251981611,-9.73692516655174,5.28638472739742,27.2474628708497]
       yp=[81.2047365447009,81.5007859895957,82.2436675361916,82.721919243158,82.4532187163616,81.6658524321853,80.1201616988957,78.868040305759,77.3745806589471,78.343810419039,79.2198008825623,81.0361747596064,82.8850251870597,84.6970166179212,88.170272215848,88.5112003265041,89.0201530656119,86.82607324433,85.3307236364029,85.7137862915936,84.4424860812103,83.6275436815664,82.8313390219846,82.2610471668164,81.2047365447009]
    elif basin_name=='arctic':
       xp=[20.0,66.0,67.0,59.0,54,63.0,95.0,150.0,-175.0,-164,-110.0,-72.0] # ,75,190,235,300,345,20]
       yp=[80.0,80.0,76.5,75.0,72,69.0,66.0,66.0 , 66.5,65.5, 66.0,   80.0]  #70,68,70,83,80,80]
    elif basin_name=='arctic_mediterranean':
       xp=[-32,-20,-14,-7,-1.5,14,20,30,63.0,95.0,150.0,-175.0,-164,-110.0,-72.0,-40.]
       yp=[ 68, 65, 65,62,60.5,61,68,66,69.0,66.0,66.0 ,  66.5,65.5,  66.0,  80.0,80.]
    elif basin_name=='eurasian_shelves':
       xp=[60,60,100,120,130,160,180,180,140,120,90,60]
       yp=[70,80,80,77,77,77,72,70,70,70,70,70]
    elif basin_name=='fram_strait':
       xp=[0.0, 10., 10., 0.0, 0.0]
       yp=[78., 78., 80., 80., 78.]
    elif basin_name=='AW_inflow1':
       xp=[20., 30., 30., 20., 20.]
       yp=[81., 81., 85., 85., 81.]
    elif basin_name=='AW_inflow2':
       xp=[80., 90., 90., 80., 80.]
       yp=[83., 83., 85., 85., 83.]
    elif basin_name=='mediterranean':
       xp=[-5,-5,12,56,59,40,33,22,-5,-5]
       yp=[36,40,50,50,33,32,29,29,31,35]
    elif basin_name=='siberia':
       xp=[30,180,180,30,30]
       yp=[60,60,70,70,60]
    elif basin_name=='namerica':
       xp=[180,290,290,180,180]
       yp=[60,60,70,70,60]
    elif basin_name=='icemelt':
       xp=[93,100,110,115,120,82,65,93] #
       yp=[82,82,81,81,85,86,85,82]
    elif basin_name=='nordic_seas':
       xp=[-32,-20,-14,-7,-1.5,14,20,17,-32,-32]
       yp=[68,65,65,62,60.5,61,69,79,79,68]
    elif basin_name=='nordic_seas_e': #Norwegian Sea
       xp=[-14,-7,-1.5,14,20,17,7,8,-8,-2,-6.5,-14]
       yp=[65,62,60.5,61,69,79,78,74,71,68,65,65]
    elif basin_name=='nordic_seas_w': #Greenland and Iceland Seas
       xp=[-32,-20,-14,-6.5,-2,-8,8,7,17,-32,-32]
       yp=[68,65,65,65,68,71,74,78,79,79,68]
    elif basin_name=='north_atlantic': #Region in the north_atlantic
       xp=[-60,-45,-45,-60,-60]
       yp=[37,37,41,41,37]
    elif basin_name=='subpolar_gyre': #Region in the subpolar_gyre
       xp=[-50,-60,-45,-40,-30,-30,-40,-50]
       yp=[ 50, 60, 59, 60, 65, 50, 50, 50]
    elif basin_name=='natlantic':
       xp=[-10.0,-20.0,-30.0,-40.0,-50.0,-60.0,-76.0,-77.0,-78.0,-82.0,-91.0,-101.0,-130.0,-164.0,-175.0,150.0,95.0,36.0,20.0,11.0, 9.0, 9.0,-5.5,-5.5,-10.0]
       yp=[  9.0,  9.0,  9.0,  9.0,  9.0,  9.0,  7.0,  8.5,  9.5,  9.0, 17.0,  21.0,  65.6,  65.5,  66.5, 66.0,66.0,62.0,68.0,60.0,55.0,49.0,36.0,35.0, 10.0]
    elif basin_name=='atlantic':
       xp=[-65.0,-76.0,-77.0,-78.0,-82.0,-91.0,-101.0,-130.0,-164.0,-175.0,150.0,95.0,36.0,20.0,11.0, 9.0, 9.0,-5.5,-5.5,25.0, 19.5, 10.0,  0.0,-13.0,-20.0,-30.0,-40.0,-50.0,-60.0,-70.0,-60.0,-65.0]
       yp=[  9.0,  7.0,  8.5,  9.5,  9.0, 17.0,  21.0,  65.6,  65.5,  66.5, 66.0,66.0,62.0,68.0,60.0,55.0,49.0,36.0,35.0, 0.0,-34.5,-40.0,-45.0,-50.0,-50.0,-50.0,-50.0,-50.0,-50.0,-50.0,-20.0,  9.0]   
    
    return xp,yp

def basin_mask(x,y,basin_name,is_basin=True,xp=None,yp=None):
    """ Find out whether lon, lat coordinated are inside a polygon """
    #if the basin polygon is not given find it based on the name
    if is_basin:
      xp,yp=basin_coords(x,y,basin_name)
    m=Basemap(projection='npstere', boundinglat=20, lon_0=0)
    xp,yp=m(xp,yp)
    poly=[]
    for j in range(len(xp)):       
       poly.append((xp[j],yp[j]))
    
    if(len(x.shape)==2):
      inside=np.ones(x.shape)*np.nan
      for i in range(x.shape[0]):
        for j in range(x.shape[1]):
          xx,yy=m(x[i,j],y[i,j])
          inside[i,j]=point_inside_polygon(xx,yy,poly)
    elif(len(x.shape)==1):
      inside=np.ones((y.shape[0],x.shape[0]))*np.nan
      for j in range(y.shape[0]):
        for i in range(x.shape[0]):
          xx,yy=m(x[i],y[j])
          inside[j,i]=point_inside_polygon(xx,yy,poly)
    
    return inside

def woa_basins(lon,lat,basin_name):
    """ READ WOA BASINS AND DEFINE BASIN MASK"""
    data=np.loadtxt('WOA_ocean_basins.txt',delimiter=',',skiprows=2,usecols=(0,1,2))
    lonwoa=data[:,1]; latwoa=data[:,0]; datawoa=data[:,2]
    mask=griddata((lonwoa,latwoa),datawoa,(lon,lat),method='nearest',fill_value=99)
    if basin_name=='arctic':
      mask[mask!=11]=0; mask[mask!=0]=1
    
    return mask

def drainage_basin(x,y,regime):
    data=np.loadtxt('Budgeteers_Arctic_Ocean_Drainage_Basin_Mask_plus_Bugeteers_Ocean_Mask_v1.0.txt',skiprows=6)
    lat=np.loadtxt('latitude_budgeteers.txt',skiprows=6)
    lon=np.loadtxt('longitude_budgeteers.txt',skiprows=6)
    data[data==-9999]=0.0
    if regime=='lnd':
      data[data==2]=0.0
    elif regime=='ocn':
      data[data==1]=0.0
      data[data==2]=1
    jj,ii=np.where(data<=2)
    mask=griddata((lon[jj,ii],lat[jj,ii]),data[jj,ii],(x,y),method='nearest',fill_value=0.0)
    
    return mask

def NorESM_masks(basin_name,iinds,jinds,barents=False,secfile='secindex.dat',grid='bipolar',lon=None,lat=None, baltic=True):
    '''NorESM_masks(basin,iinds,jinds,barents)
       Basin is either 'arctic' or 'nordic seas'. Barents is a flag, whether or not to include it to the 'arctic'. For the most parts we are following the gates defined in the secindex file. Note that in that file the indices are in Fortran 1 based format, that is why there is -1 in this function. Note that the difference to the above functions is that this is very NorESM specific and wont work with any other model while the functions above are able to figure out the points inside polygon in any model.'''
    inds_all=read_sections(secfile)
    coords=[]
    if basin_name=='atlantic':
     #Atlantic+Arctic
     if grid in ['bipolar']:
       coords.extend([[300,0],[295,46],[293,50],[293,220],[282,220],[282,225],[280,225],[280,243],[273,250],[267,250],[267,333],[70,333],[70,350],[60,350],[48,340],[48,335],[48,315],[32,300],[32,289],[32,250],[54,250],[54,0]])
     elif grid in ['tripolar']:
       lons=[-109,-100,-91,-84,-84,-81,-80,-77,-59,-59,20,20,-5.5,-5.5,8.5,94]
       lats=[63,21,17,14,11,9.5,9.5,8,-5,-34,-34,10,34,39,51,64]
       for j in range(len(lons)):
         j_out,i_out=lonlat_index(lon,lat,lons[j],lats[j],inc=1)
         coords.extend([[i_out,j_out]])
       coords.extend((inds_all.get('bering_strait')[:,:2]-1).tolist())
       lons=[-160,-115]
       lats=[66,66]
       for j in range(len(lons)):
         j_out,i_out=lonlat_index(lon,lat,lons[j],lats[j],inc=0.5)
         coords.extend([[i_out,j_out]])
       coords.extend([[coords[-1][0],lat.shape[0]]])
       coords.extend([[coords[0][0],lat.shape[0]]])
    if basin_name=='nordic_seas':
     if grid in ['bipolar']:
       coords.extend([[6,384]]) #first point
     elif grid in ['tripolar']:
       j_out,i_out=lonlat_index(lon,lat,-40,70,inc=0.5)
       coords.extend([[i_out,j_out]])
     coords.extend((inds_all.get('denmark_strait')[:,:2]-1).tolist())
     coords.extend((inds_all.get('iceland_faroe_channel')[:,:2]-1).tolist())
     coords.extend((inds_all.get('faroe_scotland_channel')[:,:2]-1).tolist())
     if grid in ['bipolar']:
       coords.extend([[36,340], [47,340]]) #from scotland to norway
     elif grid in ['tripolar']: #this is a bit different, take the North Sea and Baltic into account
       if baltic:
         coords.extend((inds_all.get('english_channel')[:,:2]-1).tolist())
         coords.extend([[coords[-1][0],coords[-1][1]-1]])
         j_out,i_out=lonlat_index(lon,lat,42,48.5,inc=0.5)
         coords.extend([[i_out,j_out]])
         j_out,i_out=lonlat_index(lon,lat,42,58.5,inc=0.5)
         coords.extend([[i_out,j_out]])
       else:
         coords.extend([[coords[-1][0],coords[-1][1]-1]])
         j_out,i_out=lonlat_index(lon,lat,15,62,inc=0.5)
         coords.extend([[i_out,j_out]])
     coords.extend((inds_all.get('barents_opening')[::-1,:2]-1).tolist())
     coords.extend((inds_all.get('fram_strait')[::-1,:2]-1).tolist())
     if grid in ['bipolar']:
       coords.extend([[115,384]]) #last point
     elif grid in ['tripolar']:
       j_out,i_out=lonlat_index(lon,lat,-40,75,inc=0.5)
       coords.extend([[i_out,j_out]])
    if basin_name=='arctic':
     if grid in ['bipolar']:
       coords.extend([[112,384]])
     elif grid in ['tripolar']:
       coords.extend((inds_all.get('canadian_archipelago')[22:,:2]-1).tolist())
       j_out,i_out=lonlat_index(lon,lat,-37.5,79.5,inc=0.5)
       coords.extend([[i_out,j_out]])
     coords.extend((inds_all.get('fram_strait')[:,:2]-1).tolist())
     if barents:
       #include barents sea
       coords.extend((inds_all.get('barents_opening')[:,:2]-1).tolist())
       if grid in ['bipolar']:
         coords.extend([[68,333]]) #
       elif grid in ['tripolar']:
         j_out,i_out=lonlat_index(lon,lat,42,58.5,inc=0.5)
         coords.extend([[i_out,j_out]])
     else:
       #don't include barents sea 
       coords.extend([[96,369],[99,369],[99,367],[103,367],[115,361],[117,359],[114,352],[104,352],[104,351],[101,351],[101,350],[98,350],[98,349],[96,349],[94,344],[94,333]])
     coords.extend((inds_all.get('bering_strait')[:,:2]-1).tolist())
     if grid in ['bipolar']:
       coords.extend([[245,333]])
       coords.extend([[245,355]])
     elif grid in ['tripolar']:
       #add points between Bering Strait and CAA
       j_out,i_out=lonlat_index(lon,lat,-155.5,66.0,inc=0.5)#,58.5,inc=0.5)
       coords.extend([[i_out,j_out]])
       j_out,i_out=lonlat_index(lon,lat,-127,66.0,inc=0.5); coords.extend([[i_out,j_out]])
       j_out,i_out=lonlat_index(lon,lat,-125,69.5,inc=0.5); coords.extend([[i_out,j_out]])
     if grid in ['bipolar']:
       coords.extend((inds_all.get('canadian_archipelago')[:2,:2]-1).tolist())
       coords.extend([[226,359],[232,362],[234,363],[236,366],[236,374],[204,375]])
       coords.extend((inds_all.get('canadian_archipelago')[2:,:2]-1).tolist())
       coords.extend([[199,384]])
     elif grid in ['tripolar']:
       coords.extend((inds_all.get('canadian_archipelago')[:22,:2]-1).tolist())
       coords.extend([[coords[-1][0],lat.shape[0]]])
       coords.extend([[coords[0][0],lat.shape[0]]])
       #j_out,i_out=lonlat_index(lon,lat,-37.5,79.5,inc=0.5)
       #coords.extend([[i_out,j_out]])
    #Then the actual mask
    ring=LinearRing(coords)
    area = Polygon(ring)
    if grid in ['bipolar']:
      mask=np.zeros((384,320))
    elif grid in ['tripolar']:
      mask=np.zeros(lat.shape)
    for j in range(len(iinds)):
       mask[jinds[j],iinds[j]]=area.contains(Point(iinds[j],jinds[j]))
    
    return mask

def lonlat_index(lon,lat,lon_p,lat_p,inc=0.5):
   """Find the closest - or at least pretty close indeces from 2D lon,lat arrays that correspond to point lon_p,lat_p. 'inc' should be in order of the grid size"""
   j1,i1=ma.where(ma.logical_and(lat<lat_p+inc,lat>lat_p-inc))
   j2=ma.where(ma.logical_and(lon[j1,i1]<lon_p+inc,lon[j1,i1]>lon_p-inc))[0]
   j_out,i_out=ma.where(ma.logical_and(lat==lat[j1,i1][j2][0],lon==lon[j1,i1][j2][0]))
   
   return j_out,i_out

def lonlatfix(lon,lat):
    """Make Longitudes to be in between -180 and 180, and output 2D lon, lat"""
    norm_grid=False
    if ma.min(lon)<-180:
      lon[lon<0.0]=lon[lon<0.0]+360
      lon[lon>180]=lon[lon>180]-360
    elif ma.min(lon)>=0.:
      lon[lon>180.0]=lon[lon>180.0]-360
      norm_grid=True
    if len(lon.shape)<2:
      lon,lat=np.meshgrid(lon, lat)
      norm_grid=False
    if ma.max(lon)>180 or ma.min(lon)<-180:
      print 'something is wrong max(lon)='+str(ma.max(lon))+', min(lon)='+ma.min(lon)
    #
    return lon, lat, norm_grid


def enable_global(tlon,tlat,data):
  """Fix NorESM/CCSM4 (only works with these models) data in such a way that it can to be plotted on a global projection on its native grid"""
  tlon = np.where(np.greater_equal(tlon,min(tlon[:,0])),tlon-360,tlon)
  tlon=tlon+abs(ma.max(tlon)); tlon=tlon+360
  # stack grids side-by-side (in longitiudinal direction), so
  # any range of longitudes may be plotted on a world map.
  tlon = np.concatenate((tlon,tlon+360),1)
  tlat = np.concatenate((tlat,tlat),1)
  data = ma.concatenate((data,data),1)
  tlon = tlon-360.
  
  return tlon, tlat, data

def NorESM_add_cyclic(lon,lat,data,dim=1):
    '''Add a cyclic point to NorESM data. Note that here we assume that data has the same size as lon and lat (so 2D). dim is the dimension in which the cyclic point is to be added. For most cases this would be the default option dim=1, but for generality the possibility is given for arbitrary dimension. If the variable is 3D with time as first variable then dim should most likely be 2'''
    #create the new shape
    dims=list(data.shape)
    dims[dim]=dims[dim]+1
    dims=tuple(dims)
    #creat the new arrays and assign the data - copy the first column/row to the end of the array
    for v in ['lon','lat','data']:
     if len(dims)<3 or v in ['lon','lat']:
      if len(dims)>2:
        exec('d'+v+'=ma.zeros((dims[1],dims[2]))')
      else:
        exec('d'+v+'=ma.zeros(dims)')
      if dim==1:
        exec('d'+v+'[:,:-1]='+v+'')
        exec('d'+v+'[:,-1]='+v+'[:,0]')
      elif dim==0:
        exec('d'+v+'[:-1,:]='+v+'')
        exec('d'+v+'[-1,:]='+v+'[0,:]')
     elif len(dims)>2 and v in ['data']:
      exec('d'+v+'=ma.zeros(dims)')
      if dim==2:
        exec('d'+v+'[:,:,:-1]='+v+'')
        exec('d'+v+'[:,:,-1]='+v+'[:,:,0]')
      elif dim==1:
        exec('d'+v+'[:,:-1,:]='+v+'')
        exec('d'+v+'[:,-1,:]='+v+'[:,0,:]')
    
    return dlon, dlat, ddata
 
#def make_dist_grid(lon,lat):
# ''' THIS IS TOO SLOW, THINK ABOUT AWAY TO AVOID THE LOOPS '''
#    if len(lon.shape)>1:
#      dx=np.zeros(lon.shape)
#      dy=np.zeros(lon.shape)
#      for i in range(1,lon.shape[0]):
#        for j in range(1,lon.shape[1]):       
#         d=gsw.earth.distance([0,lon[i,j]],[0,lat[i,j]])
#         dy[i,j]=d*np.sin(np.radians(lat[i,j]));dx[i,j]=d*np.cos(np.radians(lon[i,j]))

def across_line(lon,lat,u,v,lon_line,lat_line,FirstTime=True,iind=None,jind=None,dxout=None,dyout=None,dxin=None,dyin=None):
     """Line is a list of lon,lat points defining the line, with a straight line only start and end points are iven. U and V should be the eastward and northward transport (velocity) components. One can use the vecrotc to rotate the components if they are initially on a model grid."""
     #
     #First we find a smaller area around the line to speed up the interpolation
     if FirstTime:
       dx=2.;dy=2.;
       xp=[min(lon_line)-dx,max(lon_line)+dx,max(lon_line)+dx,min(lon_line)-dx,min(lon_line)-dx]
       yp=[min(lat_line)-dy,min(lat_line)-dy,max(lat_line)+dy,max(lat_line)+dy,min(lat_line)-dy]
       inside=basin_mask(lon,lat,'',is_basin=False,xp=xp,yp=yp)
       jind,iind=np.where(inside!=0.0)
       #find a equidistant map projection
       #m=Basemap(projection='aeqd',lon_0=min(lon_line),lat_0=np.min(lat_line))
       #xin,yin=m(lon,lat)
       #xout,yout=m(lon_line,lat_line)
       #Make a distance 'grid', check how far everything is from the equator (0,0)
       lo=lon[jind,iind]; la=lat[jind,iind]
       dxin=np.zeros(len(lo));dyin=np.zeros(len(la))
       dxout=np.zeros(len(lon_line)); dyout=np.zeros(len(lat_line))
       for j in range(len(dxin)):
         din=gsw.earth.distance([0,lo[j]],[0,la[j]])
         dxin[j]=din*np.cos(np.radians(lo[j]));dyin[j]=din*np.sin(np.radians(la[j]))
       for j in range(len(dxout)):
         dout=gsw.earth.distance([0,lon_line[j]],[0,lat_line[j]])
         dxout[j]=dout*np.cos(np.radians(lon_line[j])); dyout[j]=dout*np.sin(np.radians(lat_line[j]))
     #d=gsw.earth.distance(lon,lat)
     #distx=np.cumsum(d,axis=0);disty=np.cumsum(d,axis=1)
     #dist=np.sqrt(distx**2+disty**2)
     #Then we interpolate to get the northward and eastward components along the line
     #tx=griddata((xin[jind,iind],yin[jind,iind]),u[jind,iind],(xout,yout),method='linear',fill_value=0.0)
     #ty=griddata((xin[jind,iind],yin[jind,iind]),v[jind,iind],(xout,yout),method='linear',fill_value=0.0)
     tx=griddata((dxin,dyin),u[jind,iind],(dxout,dyout),method='linear',fill_value=0.0)
     ty=griddata((dxin,dyin),v[jind,iind],(dxout,dyout),method='linear',fill_value=0.0)
     phi=np.arctan2(np.gradient(lat_line),(np.gradient(lon_line)*np.cos(np.radians(lat_line))))
     #Finally we rotate the components to get the transports accross the line
     txr=tx*np.cos(phi)-ty*np.sin(phi)
     tyr=tx*np.sin(phi)+ty*np.cos(phi)
     #
     return txr,tyr,phi,iind,jind,dxout,dyout,dxin,dyin

def across_line2(x_line,y_line,lon_line,lat_line,u_line,v_line):
    ''' Transport across line without interpolation. Note that now lon,lat,u,v are the model variables along the line given by indices i_line and j_line. Since the line is defined along the indices we can rotate the indices just once i.e. phi is both the angle of the line and angle of the gridcells. The function above is more general and can be used to interpolate transports along any line whether or not it follows the model grid indices. However, this is probably more accurate. '''
    #phi=np.arctan2(np.gradient(lat_line),(np.gradient(lon_line)*np.cos(np.radians(lat_line))))
    phi=np.arctan2(np.gradient(y_line),np.gradient(x_line))
    #t_across=u_line*np.cos(phi)*(lon_line/np.abs(lon_line))+v_line*np.sin(phi)*(lat_line/np.abs(lat_line))
    t_across=u_line*np.cos(phi)+v_line*np.sin(phi)
    t_along=u_line*np.sin(phi)+v_line*np.cos(phi) #this is probably wrong
    
    return t_across,t_along,phi

def vecrotc(lon,lat,u,v,scalar=True):
    """ VECROTC rotate vector components [UR,VR] = VECROT(LON,LAT,U,V) rotate the vector components U and V defined at Arakawa C grid velocity points to zonal (UR) and meridional (VR) components defined at scalar points, or optionally at their respective points. LON and LAT defines the geographic location of the scalar points. Points near the pole singularities are set to NaN."""
    
    # Centered latitude and longitude differences in one direction
    dlat=np.zeros(lat.shape); dlon=np.zeros(lon.shape)
    dlat[:,0]=lat[:,1]-lat[:,0]; dlat[:,1:-1]=(lat[:,2:]-lat[:,0:-2])*.5; dlat[:,-1]=lat[:,-1]-lat[:,-2]
    dlon[:,0]=lon[:,1]-lon[:,0]; dlon[:,1:-1]=(lon[:,2:]-lon[:,0:-2])*.5; dlon[:,-1]=lon[:,-1]-lon[:,-2]
    dlon[dlon>180]=dlon[dlon>180]-360
    dlon[dlon<-180]=dlon[dlon<-180]+360
    dlon[dlon>90]=dlon[dlon>90]-180
    dlon[dlon<-90]=dlon[dlon<-90]+180
    
    # Compute rotation angle
    rad=np.pi/180
    phi=np.arctan2(dlat,(dlon*np.cos(np.radians(lat))))
    #
    us=u.copy();vs=v.copy()
    if scalar==True:
      if len(u.shape)==2:
        # Get velocity components at scalar point
        us[:,:-1]=ma.sum([u[:,:-1],u[:,1:]],0)*.5
        #us[:,-1]=u[:,-1]
        vs[:-1,:]=ma.sum([v[:-1,:],v[1:,:]],0)*.5 
        #vs[-1,:]=v[-1,:];
      elif len(u.shape)==3:
        # Get velocity components at scalar point
        us[:,:,:-1]=ma.sum([u[:,:,:-1],u[:,:,1:]],0)*.5
        #us[:,-1]=u[:,-1]
        vs[:,:-1,:]=ma.sum([v[:,:-1,:],v[:,1:,:]],0)*.5
        #vs[-1,:]=v[-1,:];
    
    # Rotate the vector components
    ur=us*np.cos(phi)-vs*np.sin(phi)
    vr=us*np.sin(phi)+vs*np.cos(phi)
    
    # Set points near pole singularities to NaN
    #ind=find(lat>88|lat<-88);
    #ur(ind)=nan;
    #vr(ind)=nan;
    
    return ur, vr, phi

def noresm2WOA(datain,shift=False, grid='gx1v6',dest='1deg'):
    '''noresm2WOA(datain,shift=False, grid='gx1v6') 
    Interpolate from the 2D datain field (choose grid to be one of the following: 'gx1v6', 'tnx1v1', 'tnx0.25v1') to WOA09 cartesian grid using a predefined weights defined in sparse matrix S'''
    #load the data
    if grid in ['gx1v6']:
      mapping_data=io.loadmat('/Home/siv22/anu074/NorESM/map_noresm_gx1v6_to_woa09_1deg_aave.mat')
      #S=mapping_data['S'].copy()
    #these two don't work, maybe something wrong with the way data is written in mat file
    elif grid in ['tnx1v1']:
      mapping_data=io.loadmat('/Data/skd/users/anu074/norstore/map_noresm_tnx1v1_to_woa09_1deg_aave_v2.mat')
      #S=mapping_data["S"].values()[0][:]
    elif grid in ['tnx0.25v1']:
      if dest in ['1deg']:
        mapping_data=io.loadmat('/Data/skd/users/anu074/norstore/map_noresm_tnx0.25v1_to_woa09_1deg_aave_v2.mat')
      elif dest in ['0.25deg']:
        mapping_data=io.loadmat('/Data/skd/users/anu074/norstore/map_noresm_tnx0.25v1_to_woa09_0_25deg_aave_v2.mat')
    S=mapping_data['S'][:]
    lon_b=mapping_data['lon_b'][:].squeeze()
    lat_b=mapping_data['lat_b'][:].squeeze()
    nx_b=mapping_data['nx_b'][:].squeeze()
    ny_b=mapping_data['ny_b'][:].squeeze()
    #
    s_a=datain.flatten()
    if False:
      ind=ma.where(1-datain.mask.flatten())[0]
      S2=S[:,ind]
      s_a2=s_a[ind]
      s_b=ma.reshape(S2*s_a2,(int(ny_b), int(nx_b)))#(int(nx_b),int(ny_b)))
    else:
      s_b=ma.reshape(S*s_a,(int(ny_b), int(nx_b)))
    #
    if shift:
      j=nx_b/2
      s_c=s_b.copy(); s_c[:,:j]=s_b[:,j:]; s_c[:,j:]=s_b[:,:j]; s_b=s_c
      lon_c=lon_b.copy(); lon_c[:j]=lon_b[j:]; lon_c[j:]=lon_b[:j]; lon_b=lon_c
    
    return lon_b, lat_b, s_b

def noresm2scalar(u,v):
    """Interpolate the velocity components to scalar points on NorESM C grid """
    us=u.copy();vs=v.copy()
    if len(u.shape)==2:
        # Get velocity components at scalar point
        us[:,:-1]=ma.sum([u[:,:-1],u[:,1:]],0)*.5
        vs[:-1,:]=ma.sum([v[:-1,:],v[1:,:]],0)*.5
    elif len(u.shape)==3:
        # Get velocity components at scalar point
        us[:,:,:-1]=ma.sum([u[:,:,:-1],u[:,:,1:]],0)*.5
        vs[:,:-1,:]=ma.sum([v[:,:-1,:],v[:,1:,:]],0)*.5
    
    return us,vs     

def scalar2u(S,dp=None,area=None):
    """ Linearly interpolate scalar S to u point on NorEsm C grid"""
    Sout=S.copy()
    if type(dp)==type(None):
      dp=np.ones(S.shape)
    if type(area)==type(None):
      area=np.ones(S.shape)
    if len(S.shape)==3 and len(area.shape)==2:
      np.tile(area,(S.shape[0],1,1))
    if len(S.shape)==3:
      Sout[:,:,1:]=ma.sum([(area*dp*S)[:,:,:-1],(area*dp*S)[:,:,1:]],0)/ma.sum([(area*dp)[:,:,:-1],(area*dp)[:,:,1:]],0)#0.5*(S[:,:,:-1]+S[:,:,1:])
      Sout[:,:,0]=ma.sum([(area*dp*S)[:,:,-1],(area*dp*S)[:,:,0]],0)/ma.sum([(area*dp)[:,:,-1],(area*dp)[:,:,0]],0) #0.5*(S[:,:,-1]+S[:,:,0])
    elif len(S.shape)==2:
      Sout[:,1:]=ma.sum([(area*S)[:,:-1],(area*S)[:,1:]],0)/ma.sum([area[:,:-1],area[:,1:]],0) #0.5*(S[:,:-1]+S[:,1:])
      Sout[:,0]=ma.sum([(area*S)[:,-1],(area*S)[:,0]],0)/ma.sum([area[:,-1],area[:,0]],0) #0.5*(S[:,-1]+S[:,0])
    #
    return Sout

def scalar2v(S,dp=None,area=None):
    """ Linearly interpolate scalar S to v point on NorEsm C grid"""
    Sout=S.copy()
    if type(dp)==type(None):
      dp=np.ones(S.shape)
    if type(area)==type(None):
      area=np.ones(S.shape)
    if len(S.shape)==3 and len(area.shape)==2:
      np.tile(area,(S.shape[0],1,1))
    if len(S.shape)==3:
      Sout[:,:-1,:]=ma.sum([(area*dp*S)[:,1:,:],(area*dp*S)[:,:-1,:]],0)/ma.sum([(area*dp)[:,1:,:],(area*dp)[:,:-1,:]],0) #0.5*(S[:,1:,:]+S[:,:-1,:])
    elif len(S.shape)==2:
      Sout[:-1,:]=ma.sum([(area*dp*S)[1:,:],(area*dp*S)[:-1,:]],0)/ma.sum([(area*dp)[1:,:],(area*dp)[:-1,:]],0) #0.5*(S[1:,:]+S[:-1,:])
    #
    return Sout

def runningMean(x, N, axis=None):
    """Fast running mean using convolve. x is the variable and N is the window size"""
    if len(x.shape)==1:
       runmean=np.convolve(x, np.ones((N,))/N)[(N-1):]
    elif len(x.shape)==2:
       runmean=np.zeros(x.shape)
       if axis==None: axis=0
       if axis==1: x=x.T; runmean=runmean.T;
       for j in range(x.shape[1]):
           runmean[:,j]=np.convolve(x[:,j].squeeze(), np.ones((N,))/N)[(N-1):]
       if axis==1: runmean=runmean.T
    return runmean

def timeMean(x,year0,inds=None,xtype='monthly',dt=None):
    """Make averages over time up to a year from annual daily data (i.e. there is 365+ entries). If longer time mean is required then calculate monthly means with this script and annual/longer averages from that. It's suggested not to give the indices, but to use 'monthly' flag or dt for all the rest."""
    dim=x.shape
    n1=365
    if type(inds)==type(None):
      if xtype in ['monthly']:
        inds=np.cumsum([0,31,28,31,30,31,30,31,31,30,31,30,31])
      else:
        inds=np.arange(0,n1+1,dt)
        inds[-1]=n1
    n=len(inds)-1
    if dt!=None and dim[0]%int(dt)==0:
      if len(dim)==1:
        y=np.nanmean(np.reshape(x,(dim[0]/dt,dt)),1)
      elif len(dim)==2:
        y=np.nanmean(np.reshape(x,(dim[0]/dt,dt,dim[1])),1)
      elif len(dim)==3:
        y=np.nanmean(np.reshape(x,(dim[0]/dt,dt,dim[1],dim[2])),1)
    else:
      years=int(x.shape[0]/365.) #how many years
      if len(dim)==1: 
        y=np.zeros(years*(len(inds)-1)) #initialize the output
        for year in range(year0,year0+years): #loop over the years
          if year==year0: 
            inds1=inds.copy() #if the first year then just copy the inds
          else:
            inds1=inds+inds1[-1] #if not the first years, then increase the indicaes
          if calendar.isleap(year):
            if n==12: #if monthly then add the leap day to february
              inds1[2:]=inds1[2:]+1
            else: #if not monthly then just use the additional day in the very end
              inds1[-1]=inds1[-1]+1
          for j in range(n-1): #loop over the time discretization (months,weeks,etc)
            y[j]=np.nanmean(x[inds1[j]:inds1[j+1]])
      else:
        y=np.zeros(list([years*(len(inds)-1)])+list(dim[1:]))
        c=0
        for year in range(year0,year0+years): #loop over the years
          if year==year0:
            inds1=inds.copy() #if the first year then just copy the inds
          else:
            inds1=inds+inds1[-1] #if not the first years, then increase the indicaes
          if calendar.isleap(year):
            if n==12: #if monthly then add the leap day to february
              inds1[2:]=inds1[2:]+1
            else: #if not monthly then just use the additional day in the very end
              inds1[-1]=inds1[-1]+1
          for j in range(n): #loop over the time discretization (months,weeks,etc)
            y[c,:]=np.nanmean(x[inds1[j]:inds1[j+1],:],0)
            c=c+1
    
    return y

def AnnualMean(x,w=None):
    """mean=AnnualMean(x,w=None)
       Calculate the annual mean from monthly timeseries, time is assumed to be the zero dimension.
       If months are of different length then use w for weights"""
    dim=x.shape
    ###############
    if len(dim)==1:
      m=ma.zeros(dim[0]/12.)
      if type(w)==type(None):
        for j in range(x.shape[0]/12):
          m[j]=ma.mean(x[j*12:(j+1)*12])
      else:
        for j in range(x.shape[0]/12):
          m[j]=ma.sum(w[j*12:(j+1)*12]*x[j*12:(j+1)*12])/ma.sum(w[j*12:(j+1)*12])
    #################
    elif len(dim)==2:
      m=ma.zeros((dim[0]/12.,dim[1]))
      if type(w)==type(None):
        for j in range(dim[0]/12):
          m[j,:]=ma.mean(x[j*12:(j+1)*12,:],0)
      else:
        w=np.tile(w,(dim[1],1)).T
        for j in range(dim[0]/12):
          m[j,:]=ma.sum(w[j*12:(j+1)*12,:]*x[j*12:(j+1)*12,:],0)/ma.sum(w[j*12:(j+1)*12,:],0)
    #################
    elif len(dim)==3:
      m=ma.zeros((dim[0]/12.,dim[1],dim[2]))
      if type(w)==type(None):
        for j in range(dim[0]/12):
          m[j,:,:]=ma.mean(x[j*12:(j+1)*12,:,:],0)
      else:
        w=np.tile(w,(dim[2],dim[1],1)).T
        for j in range(dim[0]/12):
          m[j,:,:]=ma.sum(w[j*12:(j+1)*12,:,:]*x[j*12:(j+1)*12,:,:],0)/ma.sum(w[j*12:(j+1)*12,:,:],0)
    #################
    elif len(dim)==4:
      m=ma.zeros((dim[0]/12.,dim[1],dim[2],dim[3]))
      if type(w)==type(None):
        for j in range(dim[0]/12):
          m[j,:,:,:]=ma.mean(x[j*12:(j+1)*12,:,:,:],0)
      else:
        w=np.tile(w,(dim[3],dim[2],dim[1],1)).T
        for j in range(dim[0]/12):
          m[j,:,:,:]=ma.sum(w[j*12:(j+1)*12,:,:,:]*x[j*12:(j+1)*12,:,:,:],0)/ma.sum(w[j*12:(j+1)*12,:,:,:],0)
    #
    return m

def centered_diff(t,xaxis):
    """ Use a 1D polynomial method to calculate centered difference on a irregularly spaced grid. t=input data (2D) dx=grid spacing"""
    #Initialize the data
    dx=np.diff(xaxis) #dx
    dxaxis=np.cumsum(dx) #new x-axis
    dtout=np.ones((t.shape[0],t.shape[1]-1))*np.nan #output data
    #
    #first with simple forward derivative, how about polynomial fit?
    dtout[:,0]=(t[:,1]-t[:,0])/dx[0]
    dtout[:,-1]=(t[:,-1]-t[:,-2])/dx[-1]
    #mid-points with second order approximation applying the polynomial fit
    for k in range(1,len(xaxis)-1):
       dtout[:,k]=-dx[k]*t[:,k-1]/(dx[k-1]*(dx[k]+dx[k-1]))+(dx[k]-dx[k-1])*t[:,k-1]/(dx[k]*dx[k-1])+dx[k-1]*t[:,k+1]/(dx[k]*(dx[k]+dx[k-1]))
    #
    return dtout, dxaxis

def interp_lev(data,lev, z):
    """ Use a 1D interpolation to find a certain depth surface from a 3D data: assumes, data is 3D masked array and lev is the z coorinates. z is the depth surface where we want to interpolate. Note that this should work for sigma level data as well, as far as lev=np.cumsum(dz)."""
    if np.any(lev==z):
      z0=np.nonzero(lev==z)[0][0]
      data_out=data[z0,:,:]
    else:
     #construct the 2D land mask - take the first level deeper than the level of interest
     z0=np.nonzero(lev>=z)[0][0]
     lmask=data[z0,:,:].mask.copy()
     #construct the 2D surface
     dummy=ma.array(ma.ones(lmask.shape),mask=lmask)
     dummy=dummy.flatten()
     dummy2=ma.reshape(data.flatten(),(len(lev),-1)) #original data
     indices=ma.where(~lmask.flatten())        #ocean indices
     for ii in indices[0]:
       dummy[ii]=np.interp(z,lev,dummy2[:,ii])
     #
     data_out=ma.reshape(dummy,lmask.shape)
     data_out=ma.masked_where(data_out==0,data_out) #just make sure all the points below surface are masked
     #
    return data_out

def interp_sigmacoord_to_lev(data,dz, z, jinds, iinds, bad_value, save_data=False, case=None, year=None, month=None):
    """ Use a 1D interpolation to find a certain depth surface, ie 2D field, from a 3D data defined in sigma coordinates: assumes, data is 3D array and dz is a corresponding 3D array of layer depths associated with it. z is the depth surface where we want to interpolate. Note that this is probably better for sigma coords than interp_lev because it takes into account the fact that the sigma layers have a certain depth."""
    data_out=ma.zeros((data.shape[-2],data.shape[-1]))
    data_out2=ma.zeros((data.shape[-2],data.shape[-1]))
    weights=np.zeros(dz.shape)
    zaxis=ma.cumsum(dz,0)-dz/2. #center of the cell
    zl_bound=ma.cumsum(dz,0) #Lower bound
    zu_bound=np.zeros(dz.shape)
    zu_bound[1:,:,:]=ma.cumsum(dz,0)[:-1,:,:] #Upper bound
    for h in range(len(jinds)):
       jj=jinds[h]; ii=iinds[h]
       # find all layers larger than certain cutoff thickness given as bad_value
       ind=ma.where(dz[:,jj,ii]>=bad_value)[0] #np.where(dz[:,jj,ii]!=bad_value)[0]
       #Try linear interpolation and just 'nearest neighbour' ie. find in which layer the level would be
       if (len(ind)>1 and zl_bound[:,jj,ii][-1]>z):
         #First simple linear interpolation
         #data_out[jj,ii]=np.interp(z,zaxis,data[ind,jj,ii].squeeze()) #np.interp(z,zaxis,(data[ind,jj,ii]/dz[ind,jj,ii]).squeeze())
         data_out[jj,ii]=np.interp(z,zaxis[ind,jj,ii],data[ind,jj,ii].squeeze())
         #Then do the optimal linear interpolation: find the layer in which the depth can be found and assing this depth to the value
         kk=0
         #make the kk to be the index in ind, this is needed because kk+1 and kk-1 need to be reasonable
         while ind[kk]<ind[-1] and (zl_bound[ind[kk],jj,ii]<z):
         #while kk<ind[-1] and (zl_bound[kk,jj,ii]<z):
           kk=kk+1
         #
         if ind[-1]>ind[kk] and abs(z-zaxis[ind[kk+1],jj,ii])<abs(z-zaxis[ind[kk-1],jj,ii]):
           #if level is above the last cell and the cell below is closer than the cell below
           weights[ind[kk],jj,ii]=1-abs(z-zaxis[ind[kk],jj,ii])/(abs(z-zaxis[ind[kk+1],jj,ii])+abs(z-zaxis[ind[kk],jj,ii]))
           weights[ind[kk+1],jj,ii]=1-abs(z-zaxis[ind[kk+1],jj,ii])/(abs(z-zaxis[ind[kk+1],jj,ii])+abs(z-zaxis[ind[kk],jj,ii]))
         elif zaxis[ind[kk],jj,ii]<z:
           #now the level below is closer than the level above, but check that we are not between bottom and the last cell midpoint
           weights[ind[kk],jj,ii]=1-abs(z-zaxis[ind[kk],jj,ii])/(abs(z-zaxis[ind[kk-1],jj,ii])+abs(z-zaxis[ind[kk],jj,ii]))
           weights[ind[kk-1],jj,ii]=1-abs(z-zaxis[ind[kk-1],jj,ii])/(abs(z-zaxis[ind[kk-1],jj,ii])+abs(z-zaxis[ind[kk],jj,ii]))
         elif zaxis[ind[kk],jj,ii]>z:
           #if we are between the midpoint of the last cell and bottom then just use the value of the last cell
           weights[ind[kk],jj,ii]=1
         #The output is just the weight times the data
         data_out2[jj,ii]=ma.sum(data[:,jj,ii]*weights[:,jj,ii])
         #data[ind[kk],jj,ii].squeeze() #data[ind[kk],jj,ii]/dz[ind[kk],jj,ii] #divide by layer thickness
       #If only one layer exist then can't do linear interpolation but check if level is inside this layer
       elif len(ind)==1 and (zl_bound[:,jj,ii][0]<z and zu_bound[:,jj,ii][0]>z):
         data_out2[jj,ii]=data[ind[0],jj,ii].squeeze() #data[ind[0],jj,ii]/dz[ind[0],jj,ii]
         data_out[jj,ii]=data[ind[0],jj,ii].squeeze()
         weights[ind[0],jj,ii]=1.
    if save_data:
         #Save the weights for later use as well
         exec('np.savez("/Data/skd/users/anu074/sigma_to_lev_weights/weights_'+case+'_'+str(z)+'_'+str(year)+'-'+str(month)+'.npz", weights=weights)')
     #
    return data_out, data_out2

def climatology(S):
    '''Calculate annual climatology from monthly data. Here we assume that the first dimension of S is time.'''
    tdim, jdim, idim=S.shape
    years=tdim/12
    C=ma.zeros((12,jdim,idim))
    for m in range(12):
      C[m,:,:]=ma.mean(S[m::12,:,:],0)
    
    return C

def remove_clim(S):
    """ Remove monthly climatology from the data. Here we assume that the first dimension of S is time with monthly timesteps. This function uses the climatology function to calculate the monthly climatology. """
    C=climatology(S)
    Sout=S.copy()
    tdim, jdim, idim=S.shape
    years=tdim/12
    for y in range(years):
      Sout[y*12:(y+1)*12,:,:]=S[y*12:(y+1)*12,:,:]-C
    
    return Sout

def zonal_average(data,area,lat,lon,iinds,jinds,mask=None,ave=True,dz=None,dp=None,zlevels=None,zlev=None,dlat=1,lat_out=None):
    """Calculate zonal average for Atlantic. Zonal average will have 1deg resolution in meridional direction. Note that this is rather crude way, we loop over latitudes and pick everything that is +/- .5*dlat degrees around a full degree so no fancy interpolation included. Data can be either 2D or 3D while all the other variables have to be 2D"""
    #create some variables
    if type(lat_out)==type(None):
      lat_out=np.arange(-90,90+dlat,dlat)
    dz_out=None
    if len(data.shape)>2:
      zdim=data.shape[0]
      data_out=np.zeros((zdim,len(lat_out)))
      if zdim==53:
        dz_out=data_out.copy()
    else:
       data_out=np.zeros(len(lat_out))
    #One possibility would be to make a list of longitude limits that correspond to each latitude
    #Another possibility is to define a polygon which we have used here (in NorESM_masks function)
    if type(mask)==type(None): #if there is no mask make a mask
      mask=NorESM_masks('atlantic',iinds,jinds,barents=False)
    #pick up the points that are inside the mask
    #xinds,yinds=ma.where(1-mask)
    a_area=ma.masked_array(data=area,mask=mask)
    if len(data.shape)>2:
      mask2=np.tile(mask,[zdim,1,1]) #here we tile up the mask to be same size as the variable
      a_area2=np.tile(area,[zdim,1,1]) #do the same for area
      mask3=data.mask+mask2; mask3[mask3>1]=1 #create a combined mask
      a_data=ma.masked_array(data=data,mask=mask3) #mask the variable
      a_area2=ma.masked_array(data=a_area2,mask=mask3) #mask the area
    else:
      a_data=ma.masked_array(data=data,mask=mask)
    #loop over every latitude
    if len(data.shape)>2:
      if zdim==53:
        dz=ma.masked_array(data=dz,mask=mask3) #mask the dz as well
        dz_mask=dz.mask; dz_mask[dz<1E-3]=1; dz=ma.masked_array(data=dz,mask=dz_mask) #refine the mask, take out all the cells smaller than 1mm
        a_data=ma.masked_array(data=a_data,mask=dz_mask)
        #-------------------------------------------------------------------#
        if zlevels: #This part calculates the zonal average on depth levels
          a_data2=ma.zeros((70,area.shape[0],area.shape[1]))
          for z,zz in enumerate(zlev): #first interpolate the data to z grid
            out2,out=interp_sigmacoord_to_lev(a_data,dz,zz,jinds,iinds,1E-1)
            a_data2[z,:,:]=out2
          for j,l in enumerate(lat_out):
             ii,jj=ma.where(abs(lat-l)<=0.5) #This is the full zonal band
             if ave:
                data_out[:,j]=ma.sum((a_data*a_area2)[:,ii,jj],1)/ma.sum((a_area2)[:,ii,jj],1)
             else:
                data_out[:,j]=ma.sum(a_data[:,ii,jj],1)
        #-------------------------------------------------------------------#
        else: #This part calculates the zonal average on the sigma levels
          for j,l in enumerate(lat_out):
             ii,jj=ma.where(abs(lat-l)<=dlat) #This is the full zonal band
             dz_out[:,j]=ma.sum((dz*a_area2)[:,ii,jj],1)/ma.sum(a_area2[:,ii,jj],1)
             if ave:
                data_out[:,j]=ma.sum((a_data*dp*a_area2)[:,ii,jj],1)/ma.sum((dp*a_area2)[:,ii,jj],1)
             else:
                data_out[:,j]=ma.sum(a_data[:,ii,jj],1)
        #-------------------------------------------------------------------#
      else:
        for j,l in enumerate(lat_out):
           #ii,jj=ma.where(abs(lat-l)<=0.5) #This is the full zonal band
           ii,jj=ma.where(ma.logical_and((lat-l)>=(-.5*dlat), (lat-l)<(.5*dlat)))
           #Since the latitude is already masked we can just calculate the weighted average. 
           #However, mask area again with data mask to make sure the depth is taken into account
           if ave:
             data_out[:,j]=ma.sum((a_data*a_area2)[:,ii,jj],1)/ma.sum(a_area2[:,ii,jj],1)
           else:
             data_out[:,j]=ma.sum(a_data[:,ii,jj],1) #depending on the variable one might want to calculate just the sum over the latitude band for example on energy
      #----------------------------------------------------------------------#
    else:
       for j,l in enumerate(lat_out):
         #ii,jj=ma.where(abs(lat-l)<=0.5)
         ii,jj=ma.where(ma.logical_and((lat-l)>=(-.5*dlat), (lat-l)<(.5*dlat)))
         if ave:
           data_out[j]=ma.sum((a_data*a_area)[ii,jj])/ma.sum(a_area[ii,jj])
         else:
           data_out[j]=ma.sum(a_data[ii,jj])
    
    mask3=np.zeros(data_out.shape); mask3[ma.where(data_out==0.0)]=1
    data_out=ma.masked_array(data=data_out,mask=mask3)
    
    return data_out, lat_out, dz_out


#for k in range((data.shape[0])):
           #  data_out[k,j]=ma.sum((a_data*a_area)[k,ii,jj])/ma.sum(ma.masked_array(data=a_area,mask=mask3[k,:,:].squeeze())[ii,jj])
def MLD(templvl,salnlvl,zlev,dcrit=0.03):
    '''Calculate MLD assuming density criteria of 0.03 kg/m^3'''
    salnlvl[salnlvl==0]=np.nan; templvl[templvl==0]=np.nan;
    dens=sw.pden(salnlvl, templvl, 0)
    mld=np.zeros(salnlvl.shape[0])
    for k in range(len(mld)):
      mld[k]=np.interp(dcrit, dens[k,:]-dens[k,0], zlev)
    
    return mld

def seasonind(ys,ye,season='winter'):
    '''Define a seasonindex for monthly data '''
    yy=ye-ys+1
    if season=='winter':
          seasonind=np.ones(yy*3); seasonind[0::3]=np.arange(0,12*yy,12); seasonind[1::3]=np.arange(1,12*yy,12); seasonind[2::3]=np.arange(2,12*yy,12);
    elif season=='spring':
          seasonind=np.ones(yy*3); seasonind[0::3]=np.arange(3,12*yy,12); seasonind[1::3]=np.arange(4,12*yy,12); seasonind[2::3]=np.arange(5,12*yy,12);
    elif season=='summer':
          seasonind=np.ones(yy*3); seasonind[0::3]=np.arange(6,12*yy,12); seasonind[1::3]=np.arange(7,12*yy,12); seasonind[2::3]=np.arange(8,12*yy,12);
    elif season=='autumn':
          seasonind=np.ones(yy*3); seasonind[0::3]=np.arange(9,12*yy,12); seasonind[1::3]=np.arange(10,12*yy,12); seasonind[2::3]=np.arange(11,12*yy,12);
     
    return seasonind

########################################################################################
# THE FOLLOWING ARE THE EQUATION OF STATE OF SEAWATER USED IN NORESM1-M (CMIP5 VERSION)

def eosben07_const():
    a11= 9.9985372432159340e+02;
    a12= 1.0380621928183473e+01;
    a13= 1.7073577195684715e+00;  
    a14=-3.6570490496333680e-02;
    a15=-7.3677944503527477e-03;
    a16=-3.5529175999643348e-03;
    a21= 1.0;  
    a22= 1.0316374535350838e-02;
    a23= 8.9521792365142522e-04;
    a24=-2.8438341552142710e-05;
    a25=-1.1887778959461776e-05;
    a26=-4.0163964812921489e-06;
    b11= 1.7083494994335439e-02;
    b12= 7.1567921402953455e-05;
    b13= 1.2821026080049485e-05;
    b21= 1.1995545126831476e-05;
    b22= 5.5234008384648383e-08;
    b23= 8.4310335919950873e-09;
    
    return a11,a12,a13,a14,a15,a16,a21,a22,a23,a24,a25,a26,b11,b12,b13,b21,b22,b23

def p_alpha(p,p0,th,s):
    """ Explanation of the function """
    a11,a12,a13,a14,a15,a16,a21,a22,a23,a24,a25,a26,b11,b12,b13,b21,b22,b23=eosben07_const()
    #
    b1i=1/(b11+b12*th+b13*s)
    a1=(a11+th*(a12+a14*th+a15*s)+s*(a13+a16*s))*b1i
    a2=(a21+th*(a22+a24*th+a25*s)+s*(a23+a26*s))*b1i
    b2=(b21+b22*th+b23*s)*b1i
    #
    r=b2*(p-p0)+(a2-a1*b2)*ma.log((a1+p)/(a1+p0))
    #
    return r

def eosben07_rho(p,th,s):
    """ Explanation of the function """
    a11,a12,a13,a14,a15,a16,a21,a22,a23,a24,a25,a26,b11,b12,b13,b21,b22,b23=eosben07_const()
    #
    r=(a11+th*(a12+a14*th+a15*s)+s*(a13+a16*s)+p*(b11+b12*th+b13*s))/(a21+th*(a22+a24*th+a25*s)+s*(a23+a26*s)+p*(b21+b22*th+b23*s))
    #
    return r

def eosben07_rho_th(p,th,s):
    """ Explanation of the function """
    a11,a12,a13,a14,a15,a16,a21,a22,a23,a24,a25,a26,b11,b12,b13,b21,b22,b23=eosben07_const()
    #
    P=a11+th*(a12+a14*th+a15*s)+s*(a13+a16*s)+p*(b11+b12*th+b13*s)
    Qi=1/(a21+th*(a22+a24*th+a25*s)+s*(a23+a26*s)+p*(b21+b22*th+b23*s))
    #
    r=Qi*(a12+2.0*a14*th+a15*s+b12*p-Qi*P*(a22+2.0*a24*th+a25*s+b22*p))
    #
    return r

def eosben07_rho_s(p,th,s):
    """ Explanation of the function """
    a11,a12,a13,a14,a15,a16,a21,a22,a23,a24,a25,a26,b11,b12,b13,b21,b22,b23=eosben07_const()
    #
    P=a11+th*(a12+a14*th+a15*s)+s*(a13+a16*s)+p*(b11+b12*th+b13*s);
    Qi=1/(a21+th*(a22+a24*th+a25*s)+s*(a23+a26*s)+p*(b21+b22*th+b23*s))
    #
    r=Qi*(a13+a15*th+2.0*a16*s+b13*p-Qi*P*(a23+a25*th+2.0*a26*s+b23*p))
    #
    return r

################################################################################

def read_mimoc():
    """ Read MIMOC climatology on pressure levels. p,lon,lat,saln,temp=read_mimoc()"""
    saln=np.zeros((12,81,341,720))
    temp=np.zeros((12,81,341,720))
    for j in range(12):
      fm=Dataset('/Data/skd/share/ObsData/MIMOC/pressure_grid/MIMOC_Z_GRID_v2.2_PT_S_month'+str(j+1).zfill(2)+'.nc')
      saln[j,:,:,:]=fm.variables['SALINITY'][:].copy().squeeze()
      temp[j,:,:,:]=fm.variables['POTENTIAL_TEMPERATURE'][:].copy().squeeze()
      if j==0:
        lat=fm.variables['LATITUDE'][:].copy()
        lon=fm.variables['LONGITUDE'][:].copy()
        p=fm.variables['PRESSURE'][:].copy()
    
    return p,lon,lat,saln,temp

def read_mimoc_ML(wm=False):
    """ Read MIMOC climatology on pressure levels. p,lon,lat,saln,temp=read_mimoc()"""
    saln=np.zeros((12,341,720))
    temp=np.zeros((12,341,720))
    depth=np.zeros((12,341,720))
    for j in range(12):
      if not wm:
        fm=Dataset('/Data/skd/share/ObsData/MIMOC/mixed_layer/MIMOC_ML_v2.2_PT_S_MLP_month'+str(j+1).zfill(2)+'.nc')
      else:
        fm=Dataset('/Data/skd/share/ObsData/MIMOC/mixed_layer_wm/MIMOC_ML_v2.2wm_PT_S_MLP_month'+str(j+1).zfill(2)+'.nc')
      saln[j,:,:]=fm.variables['SALINITY_MIXED_LAYER'][:].copy().squeeze()
      temp[j,:,:]=fm.variables['POTENTIAL_TEMPERATURE_MIXED_LAYER'][:].copy().squeeze()
      depth[j,:,:]=fm.variables['DEPTH_MIXED_LAYER'][:].copy().squeeze()
      if j==0:
        lat=fm.variables['LATITUDE'][:].copy()
        lon=fm.variables['LONGITUDE'][:].copy()

    return lon,lat,saln,temp,depth

def smooth2D(lon,lat,datain,n=1,use_weights=False,weights_only=False,save_weights=False,use_median=False,save_path='/Home/siv22/anu074/temp_data/JHU/'):
  """2D smoothing of (preferably masked) array datain (should be shape (lat,lon)), will be using halo of n, if n=1 (default) then the each point will be 9 point average. Option to use distance weights"""
  dataout=ma.zeros(datain.shape)
  ymax,xmax=datain.shape
  if ma.is_masked(datain):
    jind,iind=ma.where(1-datain.mask)
    dataout=ma.masked_array(dataout,mask=datain.mask)
  else:
    jind,iind=ma.where(np.ones(datain.shape))
  #weights_out=np.zeros(len(jind))
  weights_out=None #initialize output weights even if one doesn't want them
  for k in range(len(jind)):
    j=jind[k]; i=iind[k] #j=1; i=1
    jind2=[]; iind2=[]; dxy=[]
    c=0
    for ib in range(-n,n+1):
      for jb in range(-n,n+1):
        if ((j+jb)>=ymax or (j+jb)<0):
         jind2.append(j)
        else:
         jind2.append(j+jb)
        if (i+ib)>=xmax: #note that xmin case is automatically covered thanks to python indexing
         iind2.append((i+ib)-xmax)
        elif (i+ib)<0:
         iind2.append(xmax+(i+ib))
        else:
         iind2.append(i+ib)
        if datain.mask[jind2[-1],iind2[-1]]:
         jind2[-1]=j; iind2[-1]=i
        if use_weights:
          dxy.append(dist.distance([lat[j],lon[i]],[lat[jind2[c]],lon[iind2[c]]]))
        c=c+1
    if k%10000.==0:
       print k, c, j, i
    if use_weights:
      if k==0:
        weights_out=np.zeros((len(jind),c,3)) #3-d array with (weights,jind2,iind2)
      dxy=np.array(dxy)
      if ma.sum(dxy)==0:
        weights=np.ones(len(dxy))
        diind=np.argsort(dxy)
      else:
        diind=np.argsort(dxy)
        weights=(float(ma.sum(np.sort(dxy)))-np.sort(dxy))/ma.sum(float(ma.sum(np.sort(dxy)))-np.sort(dxy))
      weights_out[k,:,0]=weights
      weights_out[k,:,1]=np.array(jind2)[diind]
      weights_out[k,:,2]=np.array(iind2)[diind]
    if not weights_only:
      if use_weights:
        dataout[j,i]=ma.sum(datain[jind2[diind],iind2[diind]]*weights)/ma.sum(weights)
      elif use_median:
        dataout[j,i]=ma.median(datain[jind2,iind2]) 
      else:
        dataout[j,i]=ma.mean(datain[jind2,iind2])
  if save_weights:
    np.savez(save_path+str(n)+'_degree_smoothing_weights.npz',weights_out=weights_out,jind=jind,iind=iind)
  #
  return dataout, weights_out

def uv_interp_closest_dist(xlon,ylat,jj,ii,u,v,ulon,ulat,vlon,vlat,umask,vmask,ine,inw):
  """uu,vv=uv_interp_closest_dist(xlon,xlat,jj,ii,u,v,ulon,ulat,vlon,vlat,umask,vmask) returns weighted mean uu,vv of the 4 closest points around. jj,ii are the indices in p points, xlon and xlat are the position vectors and j defines the current position. u & v defined the velocity field that is used for the interpolation."""
  du=[]
  jjii_u=[]
  dv=[]
  jjii_v=[]
  uu=0; vv=0
  #dutot=0; dvtot=0
  #u loop over the surrounding indices, calculates the distance to the current location and saves the indices
  for hh in [jj,jj+1,jj-1]:
   if hh>383: break
   if hh<0: break
   for kk in [ii,inw[hh,ii],ine[hh,ii],ine[hh,ine[hh,ii]]]:
     du.append(dist.distance([xlon,ylat],[ulon[hh,kk],ulat[hh,kk]]))
     jjii_u.append([hh,kk])
  
  #sort the distance list and pick 4 closest
  for k in np.argsort(du)[:4]:
    jj1,ii1=jjii_u[k]
    weigth=(float(ma.sum(ma.sort(du)[:4]))-du[k])/ma.sum(float(ma.sum(np.sort(du)[:4]))-np.sort(du)[:4])
    uu=ma.sum([uu,u[jj1,ii1]*weigth],0)
    #uu=ma.sum([uu,u[jj1,ii1]*du[k]/float(np.sum(np.sort(du)[:4]))],0)
    #dutot=dutot+du[k]*umask[jj1,ii1]
  #uu=uu/dutot
  
  #v loop over the surrounding indices, calculates the distance to the current location and saves the indices
  for hh in [jj,jj+1,jj-1,jj+2]:
   if hh>383: break
   if hh<0: break
   for kk in [ii,inw[hh,ii],ine[hh,ii]]:
     dv.append(dist.distance([xlon,ylat],[vlon[hh,kk],vlat[hh,kk]]))
     jjii_v.append([hh,kk])
  
  #sort the distance list and pick 4 closest - note that the weigth is total_dist-dist/sum(total_dist) giving largest weigth to the closest point
  for k in np.argsort(dv)[:4]:
    jj1,ii1=jjii_v[k]
    weigth=(float(ma.sum(ma.sort(dv)[:4]))-dv[k])/ma.sum(float(ma.sum(np.sort(dv)[:4]))-np.sort(dv)[:4])
    vv=ma.sum([vv,v[jj1,ii1]*weigth],0) #
    #dvtot=dvtot+dv[k]*vmask[jj1,ii1]
  #vv=vv/dvtot
  
  return uu,vv

def grid_angle(lat_vertices,lon_vertices,lat,lon):
    """Calculate the grid orientation, note that this is probably not exactly correct but close enough (will be the mean of the grid vertices, but probably should be the value at the grid center)"""
    dlat1=lat_vertices[:,:,-1]-lat_vertices[:,:,0] 
    dlat2=lat_vertices[:,:,2]-lat_vertices[:,:,1]
    #note that there is a correction for the longitude differences in cases where we cross the 0 longitude
    #note that this assumes longitudes to be from 0-360
    dlon1=lon_vertices[:,:,-1]-lon_vertices[:,:,0]; dlon1[ma.where(dlon1>abs(dlon1-360))]=dlon1[ma.where(dlon1>abs(dlon1-360))]-360
    dlon2=lon_vertices[:,:,2]-lon_vertices[:,:,1];  dlon2[ma.where(dlon2>abs(dlon2-360))]=dlon2[ma.where(dlon2>abs(dlon2-360))]-360 #dlon2[ma.where(dlon2>300)]=dlon2[ma.where(dlon2>300)]-360
    #calculate the angle,  this is simply angle=(dlat/dlon*cos(lat))
    angle1=-1*np.arctan2(dlon1*np.cos(np.pi*lat/180),dlat1)
    angle2=-1*np.arctan2(dlon2*np.cos(np.pi*lat/180),dlat2)
    #return the angles of the two vertices
    return angle1,angle2

def trajectories(uin,vin,lon,lat,ulon,ulat,vlon,vlat,ine,inw,pdy,pdx,angleu=None, anglev=None,dmask=None, smax=10E3, tmax=3600*24*200, dt=3600*24, cartesian=False):
    """Drift around particles along given velocity field. The particles start from points defined by the mask. By default the timestep will be one day and the drifting time will be 200 days, but the script will adjust it's timestep along the way so that the distance between points will be indentical (give by the smax, 10km by default). Note that the angleu and anglev are the angles of the u and v faces, in NorESM (rotated c grid in general) this means that the u angle is the angle of the p points west face similarly, this means that the v angle is the angle of the p points south face, which is the same as u points east face."""
    x=[]
    y=[]
    g=9.81
    r=6371E3
    ymax,xmax=uin.shape
    ymax=ymax-1; xmax=xmax-1
    umask=uin.mask; vmask=vin.mask
    for ll in ['','u','v']:
      exec('ld='+ll+'lon')
      if ma.max(ld)>180:
         exec(ll+'lon[ma.where('+ll+'lon>180)]='+ll+'lon[ma.where('+ll+'lon>180)]-360')
    if len(umask.shape)>2:
       umask=umask[0,:,:]; vmask=vmask[0,:,:] #here we assume that the u and v masks are going to be constant in time
    jinds2,iinds2=ma.where(dmask)
    xtrack=[]; ytrack=[]
    for k in range(len(iinds2)):
     xx=iinds2[k]; yy=jinds2[k]
     print k
     if (1-umask[yy,xx]) and (1-vmask[yy,xx]):
       xlon=[lon[yy,xx]];ylat=[lat[yy,xx]]
       xlon2=[lon[yy,xx]];ylat2=[lat[yy,xx]] #these are just for picking up the monthly location
       xcum=0; ycum=0; tcum=0; j=1; tt=0
       while tcum<tmax:
         if len(uin.shape)>2: #interpolate in time, very crude as it is, but should work
            u=(1-((tcum/3600./24.)%30)/30.)*uin[tt,:,:]+(((tcum/3600./24.)%30)/30.)*uin[tt+1,:,:]
            v=(1-((tcum/3600./24.)%30)/30.)*vin[tt,:,:]+(((tcum/3600./24.)%30)/30.)*vin[tt+1,:,:]
         else:
            u=uin; v=vin
         if j>1: #check if you have crossed to the next cell
           #note that here yy and xx just keep record of p cell index. This is enough since we know the exact location in lon,lat space where we do the distance calculation
           y1=np.trunc(ycum/pdy[yy,xx]); x1=np.trunc(xcum/pdx[yy,xx])
           yy=yy+y1; xx=xx+x1
           #check if the cell is on the border -> x loops over, y stops
           if xx>xmax: #319:
             xx=0
           elif xx<0:
             xx=xmax #319
           if yy>ymax: #383:
             yy=ymax #383
             break
           elif yy<0:
             yy=0
             break
           if umask[yy,xx] or vmask[yy,xx]:
             break
           if not cartesian:
             ud,vd=uv_interp_closest_dist(xlon[j-1],ylat[j-1],yy,xx,u,v,ulon,ulat,vlon,vlat,umask,vmask,ine,inw)
           else:
             ud=u[yy,xx]; vd=v[yy,xx]
           if (np.isnan(ud) or np.isnan(vd)) or (ud==0 or vd==0):
              break
           if abs(y1)>0:
              ycum=0;
           if abs(x1)>0:
              xcum=0;
         else:
           if not cartesian:
             ud,vd=uv_interp_closest_dist(xlon[j-1],ylat[j-1],yy,xx,u,v,ulon,ulat,vlon,vlat,umask,vmask,ine,inw)
           else:
             ud=u[yy,xx]; vd=v[yy,xx]
           if (np.isnan(ud) or np.isnan(vd)) or (ud==0 or vd==0):
              break
         #adjust the timestep so that every segment will have the same length if smax is given
         #This speeds up the code remarkably - very slow points are not integrated 
         if smax!=None:
           dt=smax/np.sqrt(ud**2+vd**2)
           if tcum+dt>tmax: #we limit the time to be exactly tmax so the last timestep may need to be adjusted
              dt=tmax-tcum
           if ((tcum+dt)/3600./24.)/((tt+1)*30.)>=1: #restrict dt in such a way that there will be a one point every full month
              dt=(tt+1)*30.*24*3600-tcum
         #  tcum=tcum+dt #smax/np.sqrt(ud**2+vd**2)
         #else:
         #  tcum=tcum+dt
         tcum=tcum+dt #finally accumulate time
         #accumulate distance
         xcum=xcum+ud*dt;  ycum=ycum+vd*dt;
         #for the lon,lat coords we need to use north south velocities!
         uu=ud*np.cos(angleu[yy,xx])-vd*np.sin(anglev[yy,xx]) #this might not be really correct should be the angle in u and v points
         vv=ud*np.sin(angleu[yy,xx])+vd*np.cos(anglev[yy,xx])
         dx=uu*dt
         dy=vv*dt
         dlat=(dy/r)*(180/np.pi)
         #dlon depends on latitude - so check that the latitude is sensible (here we use the midpoint latitude
         if ylat[j-1]+dlat*.5>90:
            dummy=90-(ylat[j-1]+dlat*.5-90)
            dlon=(dx/(r*np.cos((dummy)*np.pi/180.)))*(180./np.pi)
         else:
            dlon=(dx/(r*np.cos((ylat[j-1]+dlat*.5)*np.pi/180.)))*(180./np.pi)
         if abs(dlon)>360: #make sure dlon is between -360 and 360
            dlon=360*(dlon/360.-np.trunc(dlon/360))
         #do the same check for latitude - check if crossing north pole
         if ylat[j-1]+dlat>90:
            ylat.append(90-(ylat[j-1]+dlat-90))
         else:
            ylat.append(ylat[j-1]+dlat)
         #check the -180/180 limit for longitude
         if xlon[j-1]+dlon<180 and xlon[j-1]+dlon>-180:
            xlon.append(xlon[j-1]+dlon)
         elif xlon[j-1]+dlon>180: #no larger than 180 allowed
            xlon.append(xlon[j-1]+dlon-360)
         elif xlon[j-1]+dlon<-180: #no smaller than -180 allowed
            xlon.append(xlon[j-1]+dlon+360)
         j=j+1;
         if (tcum/3600./24.)/((tt+1)*30.)>=1: #tt is month index, save only the location at the end of each month
            tt=tt+1; xlon2.append(xlon[j-1]); ylat2.append(ylat[j-1]);
       xtrack.append(xlon); ytrack.append(ylat); 
    #save the xlon2,ylat2 locations at the end of the month - this is done because I wanted to pick up a timeseries along a path
    #For plotting it is reasonable to save xlon and ylat, ie locations in given distances.
    #Maybe best to simply return both
    #############################################################################################
    #Main loop ends here
    #Finally convert the list of lists to an array (matrix) filled with nan's where there is no data
    length = len(sorted(xtrack,key=len, reverse=True)[0])
    xtrack2=np.array([xi+[np.nan]*(length-len(xi)) for xi in xtrack])
    ytrack2=np.array([yi+[np.nan]*(length-len(yi)) for yi in ytrack])
    #mask for convenience
    xmask=xtrack2.copy(); xmask[np.where(np.isnan(xtrack2))]=1; xmask[np.where(~np.isnan(xtrack2))]=0
    xtrack2=ma.masked_array(xtrack2,mask=xmask)
    ytrack2=ma.masked_array(ytrack2,mask=xmask)
    
    return xtrack2, ytrack2

def load_CMIP5(var,regime1,dims,models=None,scenario='rcp85',ensembles=None,z=[None],ymax=101):
    """all_dict=load_CMIP5(var,regime1,dims,models=None,scenario=None,ensembles=None,z=None,ymax=101)
      reads the timeseries of given variable. var, regime1, dims, and z are lists of variables, their regimes 
      (ocean, atmos, land), list of dimensions (2D or 3D neglecting time) and list of where depths at which 
      the data will be pulled. z is only relevant if the data is 3D, but one only wants certain depth level. 
      var,regime1,dims and z have to be defined. models, scenario, and ensembles can be left empty. If models
      is not defined the function will attempt to use all the available models that have data for all the variables. 
      Because of this if one wants to all the data for one model the function should be called with one variable only.
      If scenario is not specified or is None it defaults to rcp8.5 and if ensembles is not defined it  will try all the 
      possible ensemble members. ymax limits the length of the timeseries in years (relevant for long runs such as piControl)."""
    if scenario==None:
     scenarios=['piControl','historical','rcp85','rcp26','rcp45','rcp60','lgm']
     scenario=scenarios[2]; #by default it will be rcp85
    if ensembles==None:
     ensembles=['r1i1p1','r2i1p1','r3i1p1','r4i1p1','r5i1p1','r6i1p1','r7i1p1','r8i1p1','r9i1p1','r10i1p1','r11i1p1','r12i1p1','r13i1p1','r14i1p1','r1i1p150']
    regime2=[]
    for r in regime1:
       if r in ['ocean']:
         regime2.append('ModData1')
       elif r in ['atmos']:
         regime2.append('ModData2')
       elif r in ['land','seaice']:
         regime2.append('ModData3')
    if models==None:
      models=os.listdir('/Data/skd/share/'+regime2[0]+'/CMIP5/'+regime1[0]+'/'+scenario+'/'+var[0]+'/mon/')
      if len(var)>=1:
        for c,v in enumerate(var[1:]):
           modelsd=os.listdir('/Data/skd/share/'+regime2[c+1]+'/CMIP5/'+regime1[c+1]+'/'+scenario+'/'+var[c+1]+'/mon/')
           models=list(set(models)&set(modelsd)) #only keep the models that have all the variables
      models.sort() #this should be list of models that have all the variables
      if not models[0]:
         print 'There are no common models for these variables'
         return
    print models
    all_dict={'models':models} #put the list of models in the dictionary
    m_ens=[];tax=[];years=[]
    mm=0
    for m2,model in enumerate(models):
     m=m2-mm
     print model
     ens=[]
     #First read some grid variables 
     #- additionally lon,lat information is dowloaded for each variable separately
     if 'atmos' in regime1:
       if not model in os.listdir('/Data/skd/share/ModData4/CMIP5/fixed/atm/areacella/'):
          print 'There is no areacella for this model'
          all_dict['models'].pop(m); mm=mm+1
          continue
       elif not model in os.listdir('/Data/skd/share/ModData4/CMIP5/fixed/atm/sftlf/'):
          print 'There is no sftlf for this model'
          all_dict['models'].pop(m); mm=mm+1
          continue
       exec('pathg="/Data/skd/share/ModData4/CMIP5/fixed/atm/areacella/'+model+'/"')
       exec('pathga="/Data/skd/share/ModData4/CMIP5/fixed/atm/sftlf/'+model+'/"')
       dirListg=os.listdir(pathg)
       dirListga=os.listdir(pathga)
       dirListg.sort()
       dirListga.sort()
       fg=Dataset(pathg+dirListg[0])
       fga=Dataset(pathga+dirListga[0])
       a=fg.variables["areacella"][:].copy() #area of the atmospheric grid
       if model in ['CNRM-CM5']:
         a=a*100.
       aa=fga.variables["sftlf"][:].copy() #land area of the atmospheric grid [%]
       if ma.max(aa)<99.: #this is needed for MIROC-ESM and MIROC-CHEM, if not in % convert, how difficult is it to follow the protocol!!
          aa=aa*100
       lata=fg.variables["lat"][:].copy()
       lona=fg.variables["lon"][:].copy()
       lona,lata,norm_grida=lonlatfix(lona,lata)
       #put the grid variables into the dictionary
       exec('all_dict["lona'+str(m)+'"]=lona'); exec('all_dict["lata'+str(m)+'"]=lata')
       exec('all_dict["a'+str(m)+'"]=a'); exec('all_dict["aa'+str(m)+'"]=aa');
     if 'ocean' in regime1 or 'seaice' in regime1:
       if not model in os.listdir('/Data/skd/share/ModData4/CMIP5/fixed/ocn/areacello/'):
          print 'There is no areacello for this model'
          all_dict['models'].pop(m); mm=mm+1
          continue
       exec('pathgo="/Data/skd/share/ModData4/CMIP5/fixed/ocn/areacello/'+model+'/"')
       dirListgo=os.listdir(pathgo)
       dirListgo.sort()
       fgo=Dataset(pathgo+dirListgo[0])
       ao=fgo.variables["areacello"][:].copy() #area of the ocean grid
       lato=fgo.variables["lat"][:].copy()
       lono=fgo.variables["lon"][:].copy()
       lono,lato,norm_grido=lonlatfix(lono,lato)
       #put the grid variables into the dictionary
       exec('all_dict["lono'+str(m)+'"]=lono'); exec('all_dict["lato'+str(m)+'"]=lato')
       exec('all_dict["ao'+str(m)+'"]=ao');
       exec('pathdo="/Data/skd/share/ModData4/CMIP5/fixed/ocn/deptho/'+model+'/"')
       if os.path.isdir(pathdo):
         dirListdo=os.listdir(pathdo)
         fg1=Dataset(pathdo+dirListdo[0])
         deptho=fg1.variables['deptho'][:]; pmask=ma.getmask(deptho)
         if not ma.max(pmask): #if the mask was empty
           pmask=deptho.copy();pmask[deptho>0.25]=1; pmask[deptho<1]=0; pmask=1-pmask
       else:
         pmask=np.zeros(lato.shape)
       exec('all_dict["pmask'+str(m)+'"]=ma.masked_array(pmask,mask=pmask)');
     if 'seaice' in regime1:
       if model in ['CanESM2','CSIRO-Mk3-6-0']:
         exec('pathgi="/Data/skd/share/ModData4/CMIP5/fixed/atm/areacella/'+model+'/"')
         exec('pathgia="/Data/skd/share/ModData4/CMIP5/fixed/atm/sftlf/'+model+'/"')
         dirListgi=os.listdir(pathgi)
         dirListgi.sort()
         fgi=Dataset(pathgi+dirListgi[0])
         dirListgia=os.listdir(pathgia)
         dirListgia.sort()
         fgia=Dataset(pathgia+dirListgia[0])
         ai=((100-fgia.variables["sftlf"][:].copy())/100.)*fg.variables["areacella"][:].copy() #area of the atmospheric grid
       else:
         exec('pathgi="/Data/skd/share/ModData4/CMIP5/fixed/ocn/areacello/'+model+'/"')
         dirListgi=os.listdir(pathgi)
         dirListgi.sort()
         fgi=Dataset(pathgi+dirListgi[0])
         ai=fgi.variables["areacello"][:].copy() #area of the ocean grid
       lati=fgi.variables["lat"][:].copy()
       loni=fgi.variables["lon"][:].copy()
       loni,lati,norm_gridi=lonlatfix(loni,lati)
       #put the grid variables into the dictionary
       exec('all_dict["loni'+str(m)+'"]=loni'); exec('all_dict["lati'+str(m)+'"]=lati')
       exec('all_dict["ai'+str(m)+'"]=ai');
     #Check the common ensemble members
     dum=ensembles
     for c,v in enumerate(var):
       dum=list(set(dum)&set(os.listdir('/Data/skd/share/'+regime2[c]+'/CMIP5/'+regime1[c]+'/'+scenario+'/'+var[c]+'/mon/'+model+'/')))
     #loop over the ensemble members
     dum.sort()
     print dum
     for k,ensemble in enumerate(ensembles):
      if ensemble not in dum: #if the ensemble is not in the dum then break
         continue 
      ensd=[]
      for c,v in enumerate(var): #loop over the variables
        if regime1[c]=='ocean': # or regime1[c]=='seaice':
          lon=lono.copy(); lat=lato.copy()
        elif regime1[c]=='seaice':
          lon=loni.copy(); lat=lati.copy()
        elif regime1[c]=='atmos':
          lon=lona.copy(); lat=lata.copy()
        exec('pathd="/Data/skd/share/'+regime2[c]+'/CMIP5/'+regime1[c]+'/'+scenario+'/'+v+'/mon/'+model+'/'+ensemble+'/"')
        if os.path.isdir(pathd): #if ensemble exists then go forward - this should be not needed with the above loop
          exec('path'+v+str(k)+'="/Data/skd/share/'+regime2[c]+'/CMIP5/'+regime1[c]+'/'+scenario+'/'+v+'/mon/'+model+'/'+ensemble+'/"')
          exec('dirList'+v+str(k)+'=os.listdir(path'+v+str(k)+')') #these files are listed in wrong order, sort them first
          exec('dirList'+v+str(k)+'.sort()')
          exec('dirList=dirList'+v+str(k))
          ensd.append(k)
          ind=0
          ind2=0
          #############################################
          # PUT TOGETHER A LIST OF FILES FOR DOWNLOAD #
          #############################################
          #This is based on the CMIP5 naming convention makes life easy 
          #since we don't need to care about the different calendars etc
          #of infividual models- we can just sinmply pick the timeseries 
          #based on the given filenames. Note that for this reason the 
          #selection will depend on the amount of files the data is in
          #if there is only one file the function will download everything
          #no matter the selection
          if scenario in ['historical','lgm']:
            ind2=len(dirList)
            #This is away to limit the historical data to years after 1970 
            #for j in range(len(dirList)):
            #  if (int(dirList[j].split('_')[-1].split('.')[0].split('-')[0][0:4]) < ymin):
            #    ind=ind+1
          if scenario in ['piControl']: #,'lgm']:
            #limit the selection to ymax years (just put a large number if you want everything, 100 years by default)
            for j in range(len(dirList)):
              if (int(dirList[j].split('_')[-1].split('.')[0].split('-')[0][0:4])-int(dirList[0].split('_')[-1].split('.')[0].split('-')[0][0:4])<ymax):
                ind2=ind2+1;
          if 'rcp' in scenario:
            #If we are using RCP's limit the selection to years before 2101 or load until the year ymax
            for j in range(len(dirList)):
              #print j, dirList[j]
              if (int(dirList[j].split('_')[-1].split('.')[0].split('-')[0][0:4]) < 2101 and int(dirList[j].split('_')[-1].split('.')[0].split('-')[0][0:4]) < 2005+ymax):
                ind2=ind2+1
                year1=int(dirList[j].split('_')[-1].split('.')[0].split('-')[1][0:4])
          #
          year0=int(dirList[ind].split('_')[-1].split('.')[0].split('-')[0][:4]) #This is the first year
          month0=(int(dirList[ind].split('_')[-1].split('.')[0].split('-')[0][4:])-1)/12. #This is the first month (year fraction), assuming month always starts from its first day
          year1=int(dirList[ind2-1].split('_')[-1].split('.')[0].split('-')[-1][:4]) #This is the last year
          month1=int(dirList[ind2-1].split('_')[-1].split('.')[0].split('-')[-1][4:])/12. #This is the last month (year fraction), assuming month always ends to its last day
          year=year0
          years.append(np.asarray([year0,year1]))
          taxis=np.arange(year0+month0,year1+month1,1./12)
          tsum=0
          # --- -------------------------------------------------------------------------------- --- #
          # --- LOOP OVER THE FILES - SOMETIMES THERE IS MORE THAN ONE FILE - SOMETIMES ONLY ONE --- #
          # --- -------------------------------------------------------------------------------- --- #
          #INTIALIZE VARIABLES
          #then the timeseries of a full 2D field (could be for example temperature at certain layer)
          if dims[c]==2 and v not in ['msftmyz','mfo','hfbasin']:
            exec(v+'_'+str(m)+'_'+str(k)+'=np.zeros((len(taxis),lon.shape[0],lon.shape[1]))')
            exec('f_d=Dataset(path'+v+str(k)+'+dirList'+v+str(k)+'[0])')
            #if model in ['CanESM2','CSIRO-Mk3-6-0']:
            #  exec(v+'_'+str(m)+'_'+str(k)+'=np.zeros((len(taxis),f_d.dimensions["lat"].size,f_d.dimensions["lon"].size))')
            exec('lon_'+str(m)+'_'+v+',lat_'+str(m)+'_'+v+',dbool=lonlatfix(f_d.variables["lon"][:],f_d.variables["lat"][:])')
            exec('all_dict["lat_'+str(m)+'_'+v+'"]=lat_'+str(m)+'_'+v+'')
            exec('all_dict["lon_'+str(m)+'_'+v+'"]=lon_'+str(m)+'_'+v+'')
            if v in ['thetao','uo','vo','hfx','hfy']: #calculate also the relevant angles
              if 'lat_vertices' in list(f_d.variables):
                lon_vertices,lat_vertices,dbool=lonlatfix(f_d.variables["lon_vertices"][:],f_d.variables["lat_vertices"][:])
                exec('angle_'+str(m)+'_'+v+'1,angle_'+str(m)+'_'+v+'2=grid_angle(lat_vertices,lon_vertices,lat_'+str(m)+'_'+v+',lon_'+str(m)+'_'+v+')')
              elif 'lat_bnds' in list(f_d.variables): #in the case the variable is lat_bnds the grid is cartesian and the angle is 0
                exec('angle_'+str(m)+'_'+v+'1=np.zeros((lon_'+str(m)+'_'+v+'.shape[0],lon_'+str(m)+'_'+v+'.shape[1]))')
                exec('angle_'+str(m)+'_'+v+'2=np.zeros((lon_'+str(m)+'_'+v+'.shape[0],lon_'+str(m)+'_'+v+'.shape[1]))')
              exec('all_dict["angle_'+str(m)+'_'+v+'1"]=angle_'+str(m)+'_'+v+'1')
              exec('all_dict["angle_'+str(m)+'_'+v+'2"]=angle_'+str(m)+'_'+v+'2')
            f_d.close()
          #then the timeseries of full 3D field
          elif dims[c]==3 and v not in ['msftmyz','mfo','hfbasin']:
            exec('f_d=Dataset(path'+v+str(k)+'+dirList'+v+str(k)+'[0])')
            if regime1[c]=='ocean':
              lev=f_d.variables['lev'][:].copy()
            elif regime1[c]=='atmos':
              lev=f_d.variables['plev'][:].copy()
            exec(v+'_'+str(m)+'_'+str(k)+'=np.zeros((len(taxis),len(lev),lon.shape[0],lon.shape[1]))')
            exec('lon_'+str(m)+'_'+v+',lat_'+str(m)+'_'+v+',dbool=lonlatfix(f_d.variables["lon"][:],f_d.variables["lat"][:])')
            exec('all_dict["lat_'+str(m)+'_'+v+'"]=lat_'+str(m)+'_'+v+'')
            exec('all_dict["lon_'+str(m)+'_'+v+'"]=lon_'+str(m)+'_'+v+'')
            exec('all_dict["lev_'+str(m)+'_'+v+'"]=lev')
            exec('all_dict["lev_bnds_'+str(m)+'_'+v+'"]=f_d.variables["lev_bnds"][:].copy()')
            if v in ['uo','vo','umo','vmo']: #in case the variable is a vector calculate also the relevant angles in order to be able to calculate north-south vel components
              if 'lat_vertices' in list(f_d.variables):
                lon_vertices,lat_vertices,dbool=lonlatfix(f_d.variables["lon_vertices"][:],f_d.variables["lat_vertices"][:])
                exec('angle_'+str(m)+'_'+v+'1,angle_'+str(m)+'_'+v+'2=grid_angle(lat_vertices,lon_vertices,lat_'+str(m)+'_'+v+',lon_'+str(m)+'_'+v+')')
                exec('all_dict["lon_vertices_'+str(m)+'_'+v+'"]=lon_vertices')
                exec('all_dict["lat_vertices_'+str(m)+'_'+v+'"]=lat_vertices')
              elif 'lat_bnds' in list(f_d.variables): #in the case the variable is lat_bnds the grid is cartesian and the angle is 0
                exec('angle_'+str(m)+'_'+v+'1=np.zeros((lon_'+str(m)+'_'+v+'.shape[0],lon_'+str(m)+'_'+v+'.shape[1]))')
                exec('angle_'+str(m)+'_'+v+'2=np.zeros((lon_'+str(m)+'_'+v+'.shape[0],lon_'+str(m)+'_'+v+'.shape[1]))')
              exec('all_dict["angle_'+str(m)+'_'+v+'1"]=angle_'+str(m)+'_'+v+'1')
              exec('all_dict["angle_'+str(m)+'_'+v+'2"]=angle_'+str(m)+'_'+v+'2')
            f_d.close()
          #SPECIAL CASES - OVERTURNING STREAMFUNCTION AND THE TRANSPORTS ACROSS STANDARD SECTIONS
          #These are hardcoded at the moment but one could probably do something bit smarter
          #overturning stream function
          elif v in ['msftmyz']:
            exec('f_amoc=Dataset(path'+v+str(k)+'+dirList'+v+str(k)+'[0])')
            amoc_lat=f_amoc.variables["lat"][:].copy().squeeze()
            lev2=f_amoc.variables["lev"][:].copy().squeeze()
            exec(v+'_'+str(m)+'_'+str(k)+'=np.zeros((len(taxis),3,len(lev2),len(amoc_lat)))')
            exec('all_dict["amoc_lat_'+str(m)+'"]=amoc_lat')
            exec('all_dict["amoc_lev_'+str(m)+'"]=lev2')
            f_amoc.close()
          elif v in ['hfbasin']:
            exec('f_hfbasin=Dataset(path'+v+str(k)+'+dirList'+v+str(k)+'[0])')
            hfbasin_lat=f_hfbasin.variables["lat"][:].copy().squeeze()
            region=f_hfbasin.variables["region"][:]
            time_bnds=f_hfbasin.variables["time_bnds"][:].copy().squeeze()
            exec(v+'_'+str(m)+'_'+str(k)+'=np.zeros((len(taxis),len(region),len(hfbasin_lat)))')
            exec('all_dict["region_'+str(m)+'_'+str(k)+'"]=region')
            exec('all_dict["lat_hfbasin_'+str(m)+'_'+str(k)+'"]=hfbasin_lat')
            exec('all_dict["time_bnds_'+str(m)+'_'+str(k)+'"]=time_bnds')
            f_hfbasin.close()
          elif v in ['mfo']:
            exec('f_mfo=Dataset(path'+v+str(k)+'+dirList'+v+str(k)+'[0])')
            passages=f_mfo.variables["passage"][:]
            exec(v+'_'+str(m)+'_'+str(k)+'=np.zeros((len(taxis),len(passages)))')
            exec('all_dict["passages_'+str(m)+'_'+str(k)+'"]=passages')
            f_mfo.close()
          # --- ############################# --- #
          # --- # HERE STARTS THE MAIN LOOP # --- #
          # --- ############################# --- #
          cc=0 #counter for index
          for i in range(ind,ind2): #MAIN LOOP (files ie ~time)
             print str(i)
             exec('f_'+v+str(i)+'=Dataset(path'+v+str(k)+'+dirList'+v+str(k)+'[i])') #open the file
             exec('data=f_'+v+str(i)+'.variables["'+v+'"][:].squeeze()') #load the variable
             #save 2D field
             if dims[c]==2:
               if z[c]!=None and len(data.shape)>3: #if the data is in fact 3D figure out which level you want
                  if i==ind: #pickup the actual level - note that this is not interpolation - we just check in which level the given depth (height?) is and simply take that level
                   if regime1[c] in ['ocean']:
                    exec('lev=f_'+v+str(i)+'.variables["lev"][:]')
                    exec('lb=f_'+v+str(i)+'.variables["lev_bnds"][:]')
                    zz=np.where(lb[:,0]<=z[c])[0][-1] #the last level where the upper bound is smaller than the given depth
                   elif regime1[c] in ['atmos']: zz=z[c]
                  data=data[:,zz,:,:].squeeze()
               if i==ind:
                 #this is a hack to account for some models that output their u,v on different sized grid than the rest of the atmospheric output
                 nt1,ny1,nx1=data.shape
                 exec('nt2,ny2,nx2='+v+'_'+str(m)+'_'+str(k)+'.shape')
                 ny3=min(ny1,ny2);nx3=min(nx1,nx2)
               #here we save the field 
               exec(v+'_'+str(m)+'_'+str(k)+'[cc:cc+data.shape[0],:ny3,:nx3]=data')
             #save 3D field
             elif dims[c]==3:
               if i==ind:
                 nt1,nz1,ny1,nx1=data.shape
                 exec('nt2,nz2,ny2,nx2='+v+'_'+str(m)+'_'+str(k)+'.shape')
                 nz3=min(nz1,nz2);ny3=min(ny1,ny2);nx3=min(nx1,nx2)
               exec(v+'_'+str(m)+'_'+str(k)+'[cc:cc+data.shape[0],:nz3,:ny3,:nx3]=data')
             elif v in ['msftmyz','hfbasin']:
               exec(v+'_'+str(m)+'_'+str(k)+'[cc:cc+data.shape[0],:,:]=data')
             elif v in ['mfo']:
               exec(v+'_'+str(m)+'_'+str(k)+'[cc:cc+data.shape[0],:]=data')
             cc=cc+data.shape[0];
             #finally check that the next file does not overlap, for example HadGEM-ES has 1 month overlap in the turn of the 22st century - probably a bug?
             if i<ind2-1 and int(dirList[i].split('_')[-1].split('.')[0].split('-')[-1][:4])+int(dirList[i].split('_')[-1].split('.')[0].split('-')[-1][4:])/12.==int(dirList[i+1].split('_')[-1].split('.')[0].split('-')[0][:4])+int(dirList[i+1].split('_')[-1].split('.')[0].split('-')[0][4:])/12.:
                cc=cc-1
          exec('dmask='+v+'_'+str(m)+'_'+str(k)+'.copy()')
          dmask[dmask<1E20]=0; dmask[dmask>=1E20]=1
          exec('all_dict["'+v+'_'+str(m)+'_'+str(k)+'"]=ma.masked_array('+v+'_'+str(m)+'_'+str(k)+',mask=dmask)') #put the data into the dict
          tax.append(taxis)
      ens.append(ensd)
     m_ens.append(ens)
    all_dict['ens']=m_ens #m_ens is a list of a list defining which ensemble members were available
    all_dict['tax']=tax #list of time axis
    all_dict['years']=years #list of bounding years
    #Return the dictionary with all the data
    return all_dict

def increase_resolution(t1,s1,sig1,dz1,c=35, cn=15):
  """ t2,s2,sig2,dz2,zaxis2=increase_resolution(t1,s1,sig1,dz1,c=35, n=15) 
      This is a way to increase vertical resolution in MICOM - usefull if you want to create a new startup file.
      Creates new fields t2,s2,sig2,dz2 by adding n=cn-1 new layers between levels c,c+cn.
      The addition is done in sigma space and the thickness of the new layer is 1/2 of the thickness of the
      layer above and below. NOTE that cn must be 2>=2.
  """
  n=cn-1 #this is the number of new levels
  c2=c
  zn,yn,xn=s1.shape
  s2=ma.zeros((zn+n,s1.shape[1],s1.shape[2]))
  t2=ma.zeros((zn+n,t1.shape[1],t1.shape[2]))
  dz2=ma.zeros((zn+n,t1.shape[1],t1.shape[2]))
  sig2=ma.zeros((zn+n,t1.shape[1],t1.shape[2])) 
  #
  #sig11=np.zeros(zn)
  #for j in range(zn):
  #  sig11[j]=ma.max(sig1[j,:,:]) 
  #
  #sig21=np.zeros(zn+n)
  #sig21[:c]=sig11[:c]; sig21[c+2*n:]=sig11[c+n:]
  #for j in range(c,c+n):
  #  sig21[c2]=sig11[j]
  #  sig21[c2+1]=sig11[j]+(sig11[j+1]-sig11[j])/2.
  #  c2=c2+2
  #c2=c
  # Interpolate in sigma space, just linear interpolation, dz will be divided between the two new layers, have to figure out what to do with t and s.
  # Probably best to give the mean s value and iterate (if no inversion function available?) to find the temperature
  dz1_dum=dz1.copy()
  sig2[:c,:,:]=sig1[:c,:,:]; sig2[c+2*n:,:,:]=sig1[c+n:,:,:]
  s2[:c,:,:]=s1[:c,:,:]; s2[c+2*n:,:,:]=s1[c+n:,:,:]
  t2[:c,:,:]=t1[:c,:,:]; t2[c+2*n:,:,:]=t1[c+n:,:,:]
  dz2[:c,:,:]=dz1[:c,:,:]; dz2[c+2*n:,:,:]=dz1[c+n:,:,:]
  #Here is the idea: only insert a new layer with a thickness if both the layer above and layer below exist
  #otherwise just fill in a reasonable value
  #Keep track of the interface hights and adjust them accordingly.
  zaxis=ma.cumsum(dz1,0)-dz1/2. #center of the cell
  zl_bound1=ma.cumsum(dz1,0) #Lower bound
  zu_bound1=np.zeros(dz1.shape)
  zu_bound1[1:,:,:]=ma.cumsum(dz1,0)[:-1,:,:] #Upper bound  
  #upper bounds above and below the region with higher res
  zu_bound2=np.zeros(sig2.shape)
  zu_bound2[:c,:,:]=zu_bound1[:c,:,:];
  zu_bound2[c+2*n:,:,:]=zu_bound1[c+n:,:,:]
  #lower bounds above and below the region with higher res
  zl_bound2=np.zeros(sig2.shape); 
  zl_bound2[:c,:,:]=zl_bound1[:c,:,:]; 
  zl_bound2[c+2*n:,:,:]=zl_bound1[c+n:,:,:]
  #The high res area
  #initialize the old cells with the initial thickness
  zu_bound2[c:c+2*n:2,:,:]=zu_bound1[c:c+n,:,:];
  zl_bound2[c:c+2*n:2,:,:]=zl_bound1[c:c+n,:,:]  
  #initialize the new cells with zero thickness
  zu_bound2[c+1:c+2*n:2,:,:]=zl_bound1[c:c+n,:,:];
  zl_bound2[c+1:c+2*n:2,:,:]=zl_bound1[c:c+n,:,:];
  #ocean points
  jinds,iinds=np.where(dz1[0,:,:]!=0)
  #
  c2=c
  sig2[c2:c2+2*n:2,:,:]=sig1.data[c:c+n,:,:]
  s2[c2:c2+2*n:2,:,:]=s1.data[c:c+n,:,:]
  t2[c2:c2+2*n:2,:,:]=t1.data[c:c+n,:,:]
  sig2[c2+1:c2+2*n:2,:,:]=sig1.data[c:c+n,:,:]+(sig1.data[c+1:c+n+1,:,:]-sig1.data[c:c+n,:,:])/2.
  s2[c2+1:c2+2*n:2,:,:]=s1.data[c:c+n,:,:]+(s1.data[c+1:c+n+1,:,:]-s1.data[c:c+n,:,:])/2.
  t2[c2+1:c2+2*n:2,:,:]=t1.data[c:c+n,:,:]+(t1.data[c+1:c+n+1,:,:]-t1.data[c:c+n,:,:])/2.
  #loop over the ocean points
  for h in range(len(jinds)):
    yy=jinds[h]; xx=iinds[h]
    c2=c; zinds=ma.where(~dz1.mask[:,yy,xx])[0]
    #
    #sig2[c2:c2+2*n:2,yy,xx]=sig1.data[c:c+n,yy,xx]
    #s2[c2:c2+2*n:2,yy,xx]=s1.data[c:c+n,yy,xx]
    #sig2[c2+1:c2+2*n:2,yy,xx]=sig1.data[c:c+n,yy,xx]+(sig1.data[c+1:c+n+1,yy,xx]-sig1.data[c:c+n,yy,xx])/2.
    #s2[c2+1:c2+2*n:2,yy,xx]=s1.data[c:c+n,yy,xx]+(s1.data[c+1:c+n+1,yy,xx]-s1.data[c:c+n,yy,xx])/2.
    #if there are no active layers or only one then continue the loop
    #if not bool(set(range(c,c+n)) & set(zinds)) or len(set(range(c,c+n)) & set(zinds))<2:
    # continue
    #else:
    if bool(set(range(c,c+n)) & set(zinds)) or len(set(range(c,c+n)) & set(zinds))>1:
     #otherwise go on
     for j in range(c,c+n):
       if (j in list(set(range(c,c+n)) & set(zinds))) and (j+1 in list(set(range(c,c+n)) & set(zinds))): #do stuff only if the layer exists
          #the upper bound is the same as the bottom bound above
          zu_bound2[c2,yy,xx]=zl_bound2[c2-1,yy,xx]
          #rise the lower bound of the original cell by 1/4th of the grid cell depth
          zl_bound2[c2,yy,xx]=zl_bound1[j,yy,xx]-dz1[j,yy,xx]/4.
          #the lower bound of the new grid cell is 1/4 of the depth of the cell below from the original lower bound
          zl_bound2[c2+1,yy,xx]=zl_bound1[j,yy,xx]+dz1[j+1,yy,xx]/4.
          #upper bounds are copies of the lower bounds
          zu_bound2[c2+1,yy,xx]=zl_bound2[c2,yy,xx].copy()
          #adjust the lower bound whenever there is a discontinuity
          if not (j+2 in list(set(range(c,c+n)) & set(zinds))):
             zu_bound2[c2+2,yy,xx]=zl_bound2[c2+1,yy,xx].copy()
       c2=c2+2 #counter for the new layer
  dz2=zl_bound2-zu_bound2
  zaxis2=ma.cumsum(dz2,0)-dz2/2.
  mask=dz2.copy(); mask[mask!=0]=1; mask=1-mask
  dz2=ma.masked_array(dz2,mask=mask)
  sig2=ma.masked_array(sig2,mask=mask)
  t2=ma.masked_array(t2,mask=mask)
  s2=ma.masked_array(s2,mask=mask)
  #
  return t2,s2,sig2,dz2,zaxis2


def latitude_line(lat0, lat):
    """Define the indices which mark a latitude """
    iind=[]
    jind=[]
    sum=0
    #lat0=80
    i=0 #keeps track of the model i index
    i2=0 #keeps track of the length of the jind and iind
    i3=0 #helper index for looping backwards
    maxy=lat.shape[0]
    maxx=lat.shape[1]-1
    keep_looping=True
    backwards=False
    bipolar=False
    if len(np.where(np.diff(lat[-1,:])==0)[0])==1:
       bipolar=True
    #for i in range(320):
    while keep_looping: #normally loop over i indices from 0:320
      if not backwards and (lat0<max(lat[:,i]) and lat0>=min(lat[:,i])):
        #if the latitude is available append the index, this is the normal situation
        ind=np.where(lat0<=lat[:,i])[0][0] #(np.logical_and((lat-l)>=(-.5*dlat), (lat-l)<(.5*dlat)))
        jind.append(ind)
        iind.append(i)
        i=i+1; i2=i2+1; i3=i3+1
      elif len(jind)>0 and bipolar: #not (lat0<ma.max(lat[:,i:]) and lat0>=ma.min(lat[:,i:])): 
        #if the latitude doesn't exist and some indices are already there (situation close to north pole in in bipolar grid)
        #Also check that the latitude doesn't exist in the rest of the matrix (which cab be the case for the tripolar setup)
        #Then loop backwards
        if (lat0<max(lat[:,i-1]) and lat0>=min(lat[:,i-1])):
          #ind=np.round(np.interp(lat0, lat[jind[i3-1]:,i-1], np.arange(jind[i3-1],maxy)))
          ind=np.where(lat0<=lat[:,i-1])[0][-1]
          jind.append(ind)
          iind.append(i-1)
          i2=i2+1; i3=i3-1
        else:
          keep_looping=False
          #fill in the the list if needed
          if jind[-1]-jind[0]>1:
            kk=jind[-1]-jind[0]
            for k in range(kk):
              jind.append(jind[-1]-1)
              iind.append(iind[-1])
        i=i-1;
        backwards=True
      else:
        i=i+1;
      if i>maxx or i<0:
        keep_looping=False
    #
    return iind, jind

def heat_trasport(iind,jind,xtransport,ytransport):
    """ calculate the heat transport accross a given line. 
        calculate first iind and jiind. Note that this will work
        in a cartesian grid and on a NorESM type of C grid."""
    #looks already pretty good some things should be still figured out
    #First cell
    sumtot=ytransport[:,jind[0],iind[0]]
    if jind[1]>jind[0]:
          #if the next step is up right then add the transport from the cell to the right
          sumtot=ma.sum([sumtot,-1*xtransport[:,jj,ii+1]],0)
    #Last cell
    if iind[-1]==xtransport.shape[-1]-1:
    #if normal case with increasing indices
      if jind[-1]==jind[0]:
        sumtot=ma.sum([sumtot, ytransport[:,jind[-1],iind[-1]]],0)
      elif jind[-1]>jind[0]:
        sumtot=ma.sum([sumtot, ytransport[:,jind[-1],iind[-1]]+xtransport[:,jind[0],iind[0]]],0)
      elif jind[-1]<jind[0]:
        sumtot=ma.sum([sumtot, ytransport[:,jind[-1],iind[-1]]-xtransport[:,jind[0],iind[0]]],0)
    #if a tripolar grid
    elif iind[-1]>iind[-2] and jind[-1]>jind[-2]:
      sumtot=ma.sum([sumtot, ytransport[:,jind[-1],iind[-1]]-xtransport[:,jind[-1],iind[-1]]],0)
    ##########################
    # - LOOP OVER THE REST - #
    ##########################
    for j in range(1,len(jind)-1):
      #note that the last point is the copy of the first in case of bibolar
      jj=jind[j]; ii=iind[j]
      ##################################
      #Straight Line in X
      if jind[j-1]==jj and iind[j-1]<ii:
        #add the transport from the cell below
        sumtot=ma.sum([sumtot, ytransport[:,jj,ii]],0)
        if jind[j+1]>jj:
          #if the cell is last one in a strike of a cells before a step upwardright
          sumtot=ma.sum([sumtot, -1*xtransport[:,jj,ii+1]],0)
      ###################################
      #Straight backward line in x
      elif jind[j-1]==jj and iind[j-1]>ii and jj+1<ytransport.shape[1]:
        #add the transport from the cell above
        sumtot=ma.sum([sumtot, -1*ytransport[:,jj+1,ii]],0)
        if jind[j+1]<jj and iind[j+1]<ii:
          #if the cell is last one in a strike of a cells before a step downleft add the positive of xtransport
          sumtot=ma.sum([sumtot, xtransport[:,jj,ii-1]],0)
      ###################################
      #Straight line in y downwards
      if jind[j-1]>jj and iind[j-1]==ii:
         sumtot=ma.sum([sumtot, xtransport[:,jj,ii]],0)
         if iind[j+1]>ii:
           #if the cell is last one in a strike of a cells before a step right add the ytransport from below
           sumtot=ma.sum([sumtot, ytransport[:,jj,ii]],0)
      ###################################
      #Straight line in y upwards
      if jind[j-1]<jj and iind[j-1]==ii:
         sumtot=ma.sum([sumtot, -1*xtransport[:,jj,ii+1]],0)
         if iind[j+1]<ii and jj+1<xtransport.shape[-2]:
           #if the cell is last one in a strike of a cells before a step left add the ytransport from above
           sumtot=ma.sum([sumtot, -1*ytransport[:,jj+1,ii]],0)
      ###################################
      #Step down-right
      elif jind[j-1]>jj and iind[j-1]<ii:
        #add transport from the cell to the left
        sumtot=ma.sum([sumtot,xtransport[:,jj,ii]],0)
        if iind[j+1]!=ii:
          #and if the next move is away from this point ie the next cell is not the cell below
          #then add also the transport from below
          sumtot=ma.sum([sumtot,ytransport[:,jj,ii]],0)
      ####################################
      #Step upright
      elif jind[j-1]<jj and iind[j-1]<ii:
        #Add the ytransport from cell below
        sumtot=ma.sum([sumtot,ytransport[:,jj,ii]],0)
        if jind[j+1]!=jj:
          #and if the next step is not next to it then negative of the x transport from the cell to the right
          sumtot=ma.sum([sumtot,-1*xtransport[:,jj,ii+1]],0)
          if iind[j+1]<ii:
          #if the next step is step up-left (ie you're in the turning point to backward stepping)
            sumtot=ma.sum([sumtot,-1*ytransport[:,jj+1,ii]],0)
      #####################################
      #Step up-left (backwards up)
      elif jind[j-1]<jj and iind[j-1]>ii:
        #add x transport from the cell to the right 
        sumtot=ma.sum([sumtot,-1*xtransport[:,jj,ii+1]],0)
        if iind[j+1]<ii and jj+1<ytransport.shape[1]:
          #if the next step is not directly above add the transport from the cell above
          sumtot=ma.sum([sumtot,-1*ytransport[:,jj+1,ii]],0)
          if jind[j+1]<jj:
            #and if the next step is down left then add transport from the cell to the left
            sumtot=ma.sum([sumtot,xtransport[:,jj,ii]],0)
      ######################################
      #Step down-left (backwards down)
      elif jind[j-1]>jj and iind[j-1]>ii:
        #add y transport from above
        sumtot=ma.sum([sumtot,-1*ytransport[:,jj+1,ii]],0)
        if jind[j+1]<jj:
        #and if the next cell is not the cell to the left add x transport from the cell to the left
          sumtot=ma.sum([sumtot,xtransport[:,jj,ii]],0)
    #
    return sumtot

def tot_transport(lat,xtransport,ytransport,dlat=1,modeltype='NorESM'):
    """calculate the northward heat transport"""
    lat_new=np.arange(-90,90+dlat,dlat)#int(np.ceil(np.max(abs(np.diff(lat,axis=0))))))
    total=np.zeros((ytransport.shape[0], len(lat_new)))
    iinds=[]; jinds=[]
    for j,lat0 in enumerate(lat_new):
      iind,jind=latitude_line(lat0, lat)
      if modeltype=='NEMO': iind=list(np.asarray(iind)-1); jind=list(np.asarray(jind)-1);
      iinds.append(iind); jinds.append(jind)
      if len(iind)>0:
        sumtot=heat_trasport(iind,jind,xtransport,ytransport)
        total[:,j]=sumtot
      else:
        total[:,j]=0
    return total,lat_new

def linear_trend(data,monthly=False):
    """trends,r,p=linear_trend(data,monthly=False)
       return a linear trend along axis 0. Default is to assume annual data and calculate annual trend, but for monthly data the function can also return a monthly trends (monthly=True)."""
    dims=data.shape
    ##########################
    if monthly:
       months=12
    else:
       months=1
    ##########################
    if len(dims)==1:
       trends=np.zeros(months)
       r=trends.copy(); p=trends.copy()
       for j in range(months):
         slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(dims[0]/months),data[j::months])
         trends[j]=slope; r[j]=r_value; p[j]=p_value
    elif len(dims)==2:
       trends=np.zeros((months,dims[1]))
       r=trends.copy(); p=trends.copy()
       for j in range(months):
         for i in range(dims[1]):
           slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(dims[0]/months),data[j::months,i])
           trends[j,i]=slope; r[j,i]=r_value; p[j,i]=p_value
       trends=trends.squeeze(); r=r.squeeze(); p=p.squeeze()
    elif len(dims)==3:
       trends=np.zeros((months,dims[1],dims[2]))
       r=trends.copy(); p=trends.copy()
       for j in range(months):
         for i in range(dims[1]):
           for k in range(dims[2]):
             slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(dims[0]/months),data[j::months,i,k])
             trends[j,i,k]=slope; r[j,i,k]=r_value; p[j,i,k]=p_value
       trends=trends.squeeze(); r=r.squeeze(); p=p.squeeze()
    elif len(dims)==4:
       trends=np.zeros((months,dims[1],dims[2],dims[3]))
       r=trends.copy(); p=trends.copy()
       for j in range(months):
         for i in range(dims[1]):
           for k in range(dims[2]):
             for n in range(dims[3]):
               slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(dims[0]/months),data[j::months,i,k,n])
               trends[j,i,k,n]=slope; r[j,i,k,n]=r_value; p[j,i,k,n]=p_value
       trends=trends.squeeze(); r=r.squeeze(); p=p.squeeze()
    
    return trends, r, p

def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''
    
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    return segments

def create_ESMF_mask(adstr,cutoff=0.9,var='sst',infile='sst.day.mean.2000.v2.nc'):
    '''
    Create proper mask for regridded fields
    '''
    inpath='/datascope/hainegroup/anummel1/Projects/MicroInv/'+var+'_data/annual_files/'
    file1=Dataset(inpath+infile)
    sla=file1.variables[var][0,:,:].squeeze() #sla_in[n,:,:] #file1.variables['sla'][n,:,:].squeeze()
    srcmask=sla.mask
    srcgrid=grid_create_periodic([1440,720])
    if adstr in ['0.5deg']:
      dstgrid=grid_create_periodic([720,360]) #0.5 deg
    elif adstr in ['1deg']:
      dstgrid=grid_create_periodic([360,180]) #1 deg
    elif adstr in ['2deg']:
      dstgrid=grid_create_periodic([180,90])  #2 deg
    elif adstr in ['3deg']:
      dstgrid=grid_create_periodic([120,60])  #3 deg
    elif adstr in ['4deg']:
      dstgrid=grid_create_periodic([90,45])   #4 deg
    #
    file1=Dataset(inpath+infile)
    src_field_mask=ESMF.Field(srcgrid, 'mask',staggerloc=ESMF.StaggerLoc.CENTER)
    #
    src_field_mask.data[720:,:]=srcmask.T[:720,:]
    src_field_mask.data[:720,:]=srcmask.T[720:,:]
    dst_field_mask = ESMF.Field(dstgrid, 'sla_1_mask',staggerloc=ESMF.StaggerLoc.CENTER)
    regridSrc2Dst_mask = ESMF.Regrid(src_field_mask, dst_field_mask, regrid_method=ESMF.RegridMethod.CONSERVE,unmapped_action=ESMF.UnmappedAction.ERROR)
    dstfield_mask = regridSrc2Dst_mask(src_field_mask, dst_field_mask)
    dstfield_mask2=dstfield_mask.data
    dstfield_mask2[np.where(dstfield_mask2<0.9)]=0
    dstfield_mask2=np.ceil(dstfield_mask2)
    #
    return dstfield_mask2.T, dstgrid.coords[0][0][:,0].squeeze(), dstgrid.coords[0][1][0,:].squeeze()


def geographic_midpoint(lat,lon,w=None):
    '''Geographic Midpoint'''
    if w is None:
      if ma.is_masked(lat):
        w=ma.masked_array(np.ones(lat.shape),mask=lat.mask)
      else:
        w=np.ones(lat.shape)
    x=ma.sum(np.cos(lat*np.pi/180.)*np.cos(lon*np.pi/180.)*w,0)/ma.sum(w,0)
    y=ma.sum(np.cos(lat*np.pi/180.)*np.sin(lon*np.pi/180.)*w,0)/ma.sum(w,0)
    z=ma.sum(np.sin(lat*np.pi/180.)*w,0)/ma.sum(w,0)
    #
    lon_out=np.arctan2(y,x)*180./np.pi
    lat_out=np.arctan2(z,np.sqrt(x*x+y*y))*180./np.pi
    #
    return lat_out,lon_out
