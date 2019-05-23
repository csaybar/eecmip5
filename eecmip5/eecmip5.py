# -*- coding: utf-8 -*-

"""Main module."""

#!/usr/bin/env python
import ee
import pandas as pd
from datetime import datetime
import os, sys, re

class CMIP5(object):
  def __init__(self,variable,model,geometry,download,folder):

    self.ee_connection()        
    self.download = download
    self.variable = variable
    self.model = model 
    self.geometry = self.eegeometry(geometry)    
    self.folder = folder
    self.units()
        
  def eegeometry(self,geom):
    '''
    Create an ee.Geometry.Rectangle from a list

    Parameters:
    -----------

    - geom: The minimum and maximum corners of the rectangle, as a list of
    two points each in the format of GeoJSON 'Point' coordinates, or a
    list of two ee.Geometry describing a point, or a list of four
    numbers in the order xMin, yMin, xMax, yMax.

    Returns:
    --------
    An ee.Geometry.Rectangle object
    '''        
    return ee.Geometry.Rectangle(geom[0],geom[1],geom[2],geom[3])
      
  def units(self):
    '''
    Change PP units of kg/(m^2*s) to mm and T of Kelvin to Celsius
    '''    
    if self.variable is 'pr':
      self.mltpy = 86400
      self.add = 0
    elif self.variable is 'tasmin' or self.variable is 'tasmax':
      self.mltpy = 1
      self.add = -273.15
    else:
      raise ValueError( self.variable + ' does not supported!')
            
  def ee_connection(self):
    '''
    Verify if exist a connection with the Earth Engine server
    '''
    try:
      ee.Initialize()
      print('The Earth Engine package initialized successfully!')
    except ee.EEException as error:
      print('The running was stopped:')
      print(error)
      
  def accumulate(self,image,img):
    '''
    Iterator function
    '''
    name_image = image.get('system:index')
    image = image.select([0],[name_image])
    cumm = ee.Image(img).addBands(image)
    return cumm
    
class historical(CMIP5):
  def __init__(self,variable,model,geometry,
               download = 'drive',
               historical_start = '1950-01-01',
               historical_end = '2006-01-01',
               folder = 'CMIP5'):
    
    CMIP5.__init__(self,variable,model,geometry,download,folder)                
    self.historical_start = datetime.strptime(historical_start, '%Y-%m-%d')
    self.historical_end = datetime.strptime(historical_end, '%Y-%m-%d')
    self.scenario = 'historical'
            
  def daily_to_annual(self):
  
    CMIP5 = ee.ImageCollection('NASA/NEX-GDDP') \
              .filterDate(self.historical_start,self.historical_end) \
              .select(self.variable) \
              .filter(ee.Filter.eq('scenario', self.scenario)) \
              .filter(ee.Filter.eq('model', self.model))

    years = ee.List.sequence(self.historical_start.year,self.historical_end.year-1)

    
    if self.variable is 'pr':
      byyear = ee.ImageCollection.fromImages(years.map(lambda x: CMIP5.filter(ee.Filter.calendarRange(x,x,'year'))\
                 .sum()\
                 .set('year',x)))      
    else:
      byyear = ee.ImageCollection.fromImages(years.map(lambda x: CMIP5.filter(ee.Filter.calendarRange(x,x,'year'))\
               .mean()\
               .set('year',x)))
    
    timeseries = ee.FeatureCollection(byyear.map(lambda img: img.multiply(self.mltpy)\
                                                                .add(self.add)\
                                                                .reduceRegions(self.geometry,ee.Reducer.mean(),500)\
                                                                .first()\
                                                                .set('year',img.get('year'))))
    if self.download == 'URL':
      print(timeseries.getDownloadURL('csv',['year','mean']))                
    elif self.download == 'drive':
      desc_name = '%s%s%s%s' % (self.variable,self.scenario,self.model,
                                datetime.now().strftime("%Y%m%d%H%M%S"))
      my_task = ee.batch.Export.table.toDrive(
                        collection = timeseries,
                        fileFormat='csv',
                        folder = self.folder,
                        description = desc_name,
                        selectors=['year','mean']).start()
      return my_task
    else:
      raise ValueError( self.download + ' does not supported!')    

  def download_data(self):
    time_series = pd.date_range(self.historical_start, 
                                self.historical_end,
                                freq='YS')
    time_series = time_series[0:(len(time_series)-1)]

    for time in time_series:
      str_date = str(time)[0:10]
      end_date = str(time +1)[0:10]

      CMIP5 = ee.ImageCollection('NASA/NEX-GDDP') \
                .filterDate(str_date,end_date) \
                .select(self.variable) \
                .filter(ee.Filter.eq('scenario', self.scenario)) \
                .filter(ee.Filter.eq('model', self.model))                       

      col_band = CMIP5.map(lambda img: img.multiply(self.mltpy)\
                                          .add(self.add)\
                                          .set('system:time_start', img.get('system:time_start'))\
                                          .set('system:index', img.get('system:index')))

      #  ImageCollection to List           
      col_list = col_band.toList(col_band.size())

      #  Define the initial value for iterate.
      base = ee.Image(col_list.get(0))
      base_name = base.get('system:index')
      base = base.select([0], [base_name])

      #  Eliminate the image 'base'.
      new_col = ee.ImageCollection(col_list.splice(0,1))
      img_cummulative = ee.Image(new_col.iterate(self.accumulate,base))

      print("Downloading data for the date: " + str_date)

      if self.download == 'URL':
        print(img_cummulative.clip(self.geometry)\
                             .getDownloadURL({'scale':25000,'crs':'EPSG:4326'}))

      elif self.download == 'drive':
        ee.batch.Export.image.toDrive(image = img_cummulative.clip(self.geometry),
                folder = self.folder,
                fileNamePrefix = '%s%s%s%s' % (self.model,self.variable,self.scenario,str_date),
                scale = 25000).start()
      else:
        raise ValueError( self.download + ' does not supported!')
        

class rcp45(CMIP5):
  def __init__(self,variable,model,geometry,
               download = 'drive',
               scenario_start = '2006-01-01',
               scenario_end = '2100-12-31',
               folder = 'CMIP5'):

    CMIP5.__init__(self,variable,model,geometry,
                   download,folder)
    
    self.scenario = 'rcp45'
    self.scenario_start = datetime.strptime(scenario_start, '%Y-%m-%d')
    self.scenario_end = datetime.strptime(scenario_end, '%Y-%m-%d')
    
  def daily_to_annual(self):    
    CMIP5 = ee.ImageCollection('NASA/NEX-GDDP') \
              .filterDate(self.scenario_start,self.scenario_end) \
              .select(self.variable) \
              .filter(ee.Filter.eq('scenario', self.scenario)) \
              .filter(ee.Filter.eq('model', self.model))
    
    years = ee.List.sequence(self.scenario_start.year,self.scenario_end.year-1)
    
    if self.variable is 'pr':
      byyear = ee.ImageCollection.fromImages(years.map(lambda x: CMIP5.filter(ee.Filter.calendarRange(x,x,'year'))\
                 .sum()\
                 .set('year',x)))      
    else:
      byyear = ee.ImageCollection.fromImages(years.map(lambda x: CMIP5.filter(ee.Filter.calendarRange(x,x,'year'))\
               .mean()\
               .set('year',x)))    
      
    timeseries = ee.FeatureCollection(byyear.map(lambda img: img.multiply(self.mltpy)\
                                                                .add(self.add)\
                                                                .reduceRegions(self.geometry,ee.Reducer.mean(),500)\
                                                                .first()\
                                                                .set('year',img.get('year'))))                    
    if self.download == 'URL':
      print(timeseries.getDownloadURL('csv',['year','mean']))
                
    elif self.download == 'drive':
      desc_name = '%s%s%s%s' % (self.variable,self.scenario,self.model,
                                      datetime.now().strftime("%Y%m%d%H%M%S"))
      my_task = ee.batch.Export.table.toDrive(
                          collection = timeseries,
                          fileFormat='csv',
                          folder = self.folder,
                          description = desc_name,
                          selectors=['year','mean']).start()                
      return my_task
    else:
      raise ValueError( self.download + ' does not supported!')
      
  def download_data(self):
    
    time_series = pd.date_range(self.scenario_start, 
                                self.scenario_end,
                                freq='YS')
    
    time_series = time_series[0:(len(time_series)-1)]
    url_list = []

    for time in time_series:

      str_date = str(time)[0:10]
      end_date = str(time +1)[0:10]
                
      CMIP5 = ee.ImageCollection('NASA/NEX-GDDP') \
                .filterDate(str_date,end_date) \
                .select(self.variable) \
                .filter(ee.Filter.eq('scenario', self.scenario)) \
                .filter(ee.Filter.eq('model', self.model))                       
      
      col_band = CMIP5.map(lambda img: img.multiply(self.mltpy)\
                                          .add(self.add)\
                                          .set('system:time_start', img.get('system:time_start'))\
                                          .set('system:index', img.get('system:index')))

      #  ImageCollection to List           
      col_list = col_band.toList(col_band.size())

      #  Define the initial value for iterate.
      base = ee.Image(col_list.get(0))
      base_name = base.get('system:index')
      base = base.select([0], [base_name])
    
      #  Eliminate the image 'base'.
      new_col = ee.ImageCollection(col_list.splice(0,1))

      img_cummulative = ee.Image(new_col.iterate(self.accumulate,base))

      print("Downloading data for the date: " + str_date)
                
      if self.download == 'URL':
        list_download = img_cummulative.clip(self.geometry)\
                                                 .getDownloadURL({'scale':25000,'crs':'EPSG:4326'})
        url_list.append(list_download)

      elif self.download == 'drive':
        ee.batch.Export.image.toDrive(
            image = img_cummulative.clip(self.geometry),
            folder = self.folder,
            fileNamePrefix = '%s%s%s%s' % (self.model,self.variable,self.scenario,str_date),
            scale = 25000).start()
      else:
        raise ValueError( self.download + ' does not supported!')
        
      if self.download == 'URL':
        return url_list


class rcp85(CMIP5):
  def __init__(self,variable,model,geometry,
               download = 'drive',
               scenario_start = '2006-01-01',
               scenario_end = '2100-12-31',
               folder = 'CMIP5'):

    CMIP5.__init__(self,variable,model,geometry,download,folder)
                
    self.scenario = 'rcp85'
    self.scenario_start = datetime.strptime(scenario_start, '%Y-%m-%d')
    self.scenario_end = datetime.strptime(scenario_end, '%Y-%m-%d')
    
  def daily_to_annual(self):
    CMIP5 = ee.ImageCollection('NASA/NEX-GDDP') \
              .filterDate(self.scenario_start,self.scenario_end) \
              .select(self.variable) \
              .filter(ee.Filter.eq('scenario', self.scenario)) \
              .filter(ee.Filter.eq('model', self.model))        
    years = ee.List.sequence(self.scenario_start.year,self.scenario_end.year-1)
    
    if self.variable is 'pr':
      byyear = ee.ImageCollection.fromImages(years.map(lambda x: CMIP5.filter(ee.Filter.calendarRange(x,x,'year'))\
                 .sum()\
                 .set('year',x)))      
    else:
      byyear = ee.ImageCollection.fromImages(years.map(lambda x: CMIP5.filter(ee.Filter.calendarRange(x,x,'year'))\
               .mean()\
               .set('year',x)))
            
    timeseries = ee.FeatureCollection(byyear.map(lambda img: img.multiply(self.mltpy)\
                                                                .add(self.add)\
                                                                .reduceRegions(self.geometry,ee.Reducer.mean(),500)\
                                                                .first()\
                                                                .set('year',img.get('year'))))
                    
    if self.download == 'URL':
      print(timeseries.getDownloadURL('csv',['year','mean']))
                
    elif self.download == 'drive':
      desc_name = '%s%s%s%s' % (self.variable,self.scenario,self.model,
                                datetime.now().strftime("%Y%m%d%H%M%S"))
      my_task = ee.batch.Export.table.toDrive(
                            collection = timeseries,
                            fileFormat='csv',
                            folder = self.folder,
                            description = desc_name,
                            selectors=['year','mean']).start()
      return my_task
    else:
      raise ValueError( self.download + ' does not supported!')

    def download_data(self):

        time_series = pd.date_range(self.scenario_start, 
                                        self.scenario_end,
                                        freq='YS')
                                        
        time_series = time_series[0:(len(time_series)-1)]
        url_list = []

        for time in time_series:

            str_date = str(time)[0:10]
            end_date = str(time +1)[0:10]
                
            CMIP5 = ee.ImageCollection('NASA/NEX-GDDP') \
                      .filterDate(str_date,end_date) \
                      .select(self.variable) \
                      .filter(ee.Filter.eq('scenario', self.scenario)) \
                      .filter(ee.Filter.eq('model', self.model))                       

            col_band = CMIP5.map(lambda img: img.multiply(self.mltpy)\
                                                .add(self.add)\
                                                .set('system:time_start', img.get('system:time_start'))\
                                                .set('system:index', img.get('system:index')))

            #  ImageCollection to List           
            col_list = col_band.toList(col_band.size())

            #  Define the initial value for iterate.
            base = ee.Image(col_list.get(0))
            base_name = base.get('system:index')
            base = base.select([0], [base_name])
    
            #  Eliminate the image 'base'.
            new_col = ee.ImageCollection(col_list.splice(0,1))

            img_cummulative = ee.Image(new_col.iterate(self.accumulate,base))

            print("Downloading data for the date: " + str_date)
                
            if self.download == 'URL':
                list_download = img_cummulative.clip(self.geometry)\
                                               .getDownloadURL({
                                                    'scale':25000,
                                                    'crs':'EPSG:4326'})
                url_list.append(list_download)

            elif self.download == 'drive':
                ee.batch.Export.image.toDrive(
                        image = img_cummulative.clip(self.geometry),
                        folder = self.folder,
                        fileNamePrefix = '%s%s%s%s' % (self.model,self.variable,self.scenario,str_date),
                        scale = 25000).start()
            else:
                raise ValueError( self.download + ' does not supported!')
        
        if self.download == 'URL':
            return url_list
