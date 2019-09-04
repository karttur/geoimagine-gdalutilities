'''
Created on 22 Nov 2018

@author: thomasgumbricht
'''

__version__ = '0.3.1'
VERSION = tuple( int(x) for x in __version__.split('.') )
metadataD = { 'name':'gdalutilities', 
             'author':'Thomas Gumbricht', 
             'author_email':'thomas.gumbricht@gmail.com',
             'title':'GDAL utilities', 'label':'Intarface to GDAL utilities.',
             'prerequisites':'GDAL utilities must be available via the operating system terminal',
             'image':'avg-trmm-3b43v7-precip_3B43_trmm_2001-2016_A'}