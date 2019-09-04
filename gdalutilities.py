'''
Created on 22 Nov 2018

@author: thomasgumbricht
'''

import os
from sys import exit
import geoimagine.gis.mj_gis_v80 as mj_gis 

GDALpath = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs'

class ProcessGdalUtilities:
    '''class for GDAL processing'''   
    def __init__(self, process, session, verbose):
        self.session = session
        self.verbose = verbose
        self.process = process
        if len(self.process.srcLayerD) == len(self.process.dstLayerD) == 1 and list(self.process.srcLayerD.keys())[0] != list(self.process.dstLayerD.keys())[0]:
            self._LoopWithRegionChange()
        else:
            SNULLE

    def _LoopWithRegionChange(self):
        #differet src and dst regions
        srclocus = list(self.process.srcLayerD.keys())[0]
        dstlocus = list(self.process.dstLayerD.keys())[0]
        #Get the region for the dstLocus from the db
        queryD = {}
        queryD['regionid'] = {'val':dstlocus, 'op':'=' }
        paramL = ['epsg', 'ullat', 'ullon', 'urlat', 'urlon', 'lrlat', 'lrlon', 'lllat', 'lllon', 'minx', 'miny', 'maxx', 'maxy']

        coords = self.session._SelectDefRegionExtent(queryD,paramL)
        if coords == None:
            exitstr = 'The dst_region %(regionid)s does not exist, cannot translate using gdal utlities' %queryD
            exit(exitstr)
        self.extentD = dict(zip(paramL,coords))
        for datum in self.process.dstLayerD[dstlocus]:
            cont = True
            for dstcomp in self.process.dstLayerD[dstlocus][datum]:
                if not self.process.dstLayerD[dstlocus][datum][dstcomp]._Exists() or self.process.overwrite:
                    dstLayer = self.process.dstLayerD[dstlocus][datum][dstcomp]
                else:
                    self.session._InsertLayer(self.process.dstLayerD[dstlocus][datum][dstcomp], self.process.overwrite, self.process.delete)
                    cont = False
            if cont:
                
                for srccomp in self.process.dstLayerD[dstlocus][datum]:
                    srcLayer = self.process.srcLayerD[srclocus][datum][srccomp]   
                    self.srcMeta = mj_gis.RasterGetLayerMeta(self.process.srcLayerD[srclocus][datum][srccomp].FPN)
  
                self.ulx = self.extentD['ullon']
                self.uly = self.extentD['ullat']
                self.lrx = self.extentD['lrlon']
                self.lry = self.extentD['lrlat']
  
                if self.process.proc.processid in ['gdal_translateancillary']:
                    self._Gdal_Translate(srcLayer,dstLayer,True)
                else:
                    exitStr = 'EXITING, unknown process in ProcessGdalUtilities (gdalutilities.py: %(s)s)' %{'s':self.process.proc.processid}
                    exit(exitStr)
                self.session._InsertLayer(dstLayer, self.process.proc.overwrite, self.process.proc.delete)
      
    def _Gdal_Translate(self, srcLayer, dstLayer, projwin):
        '''
        '''
        tr = True
        scale = True
        outsize = True
        resampling = True
        
        if self.process.params.xsize == self.process.params.ysize == 0:
            outsize = False
        else:
            tr = False
            
        if self.process.params.xres == self.process.params.yres == 0:
            tr = False
        if self.process.params.src_min == self.process.params.dst_min and self.process.params.src_max == self.process.params.dst_max:
            scale = False 
        if resampling == 'nearest':
            resampling = False
 
        gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_translate'
        
        #if srcLayer.comp.hdrfiletype.lower() not in ['tif','.tif']:
        #    gdalcmd = '%(cmd)s -of' %{'cmd':gdalcmd}
        gdalcmd = '%(cmd)s -ot %(ot)s' %{'cmd':gdalcmd, 'ot':dstLayer.comp.celltype}
        if outsize:
            gdalcmd = '%(cmd)s -outsize %(xsize)s %(ysize)s' %{'cmd':gdalcmd, 'xsize':self.process.params.xsize, 'ysize':self.process.params.ysize}

        if tr:
            gdalcmd = '%(cmd)s -tr %(xres)f %(yres)f' %{'cmd':gdalcmd, 'xres':self.xres, 'yres':self.yres}
        if projwin:
            gdalcmd = '%(cmd)s -projwin %(ulx)f %(uly)f %(lrx)f %(lry)f' %{'cmd':gdalcmd, 'ulx':self.ulx, 'uly':self.uly, 'lrx':self.lrx, 'lry':self.lry}

        #resamplig reqires GDAL 2, change
        if resampling:
            gdalcmd = '%(cmd)s -r %(resampling)s' %{'cmd':gdalcmd, 'resampling':self.process.params.resampling}

        if scale:
            gdalcmd = '%(cmd)s -scale %(src_min)d %(src_max)d %(dst_min)d %(dst_max)d' %{'cmd':gdalcmd, 
                        'src_min':self.process.params.src_min, 'src_max':self.process.params.src_max, 
                        'dst_min':self.process.params.dst_min, 'dst_max':self.process.params.dst_max}
        if self.process.params.exponent:
            gdalcmd = '%(cmd)s -exponent %(exp)f' %{'cmd':gdalcmd, 'exp':self.exponent}
        gdalcmd = '%(cmd)s %(src)s %(dst)s' %{'cmd':gdalcmd, 
                        'src':srcLayer.FPN, 'dst':dstLayer.FPN}
        os.system(gdalcmd)
        
        
    
    
    def _GdalWarp(self):
        pass
           
    def _SetFGBGparams(self): 
        #TGTODO ASSEMBEL ALL PARAMETERS IN A SEPARATE CLASS
        from math import atan, sin, cos
        for locus in self.process.srcLayerD:
            for datum in self.process.srcLayerD[locus]:
                dstLayerD = {}
                for dstcomp in self.process.dstLayerD[locus][datum]:
                    layerId = self.process.dstLayerD[locus][datum][dstcomp].comp.id
                    dstLayerD[layerId] = []
                srcLayerD = {}
                for srccomp in self.process.srcLayerD[locus][datum]:
                    if not (self.process.srcLayerD[locus][datum][srccomp]):
                        continue
                    layerId = self.process.srcLayerD[locus][datum][srccomp].comp.id
                    srcLayerD[layerId] = []
                if (self.process.srcLayerD[locus][datum][srccomp]):
                    break
        
        srcKeys = list(srcLayerD.keys())
        if 'xband' in srcKeys:
            self.xb = 'xband'
        elif 'bandx' in srcKeys:
            self.xb = 'bandx'
        else:
            exitstr= 'the fgbg process must have either xband or bandx as srcband id'
            exit(exitstr)
        if 'yband' in srcKeys:
            self.yb = 'yband'
        elif 'bandy' in srcKeys:
            self.yb = 'bandy'
        else:
            exitstr= 'the fgbg process must have either yband or bandy as srcband id'
            exit(exitstr)
        dstKeys = list(dstLayerD.keys())
        if 'fg' or 'bg' in dstKeys:
            if 'fg' in dstKeys:
                self.fg = 'fg'
            else:
                self.fg = False
            if 'bg' in dstKeys:
                self.bg = 'bg'
            else:
                self.bg = False
        else:
            exitstr= 'the fgbg process must have either fg or bg (or both) as dstband id'
            exit(exitstr)
        
        #Do the rotation
        angrad = -atan(self.process.params.slope)
        rangdeg = 180 * angrad / 3.1415
        rangdeg += 45
        #Convert degrees to radians
        rangrad = 3.1415 * rangdeg / 180
        #Get the sin and cos angles
        self.sinrang = sin(rangrad) 
        self.cosrang = cos(rangrad)
                              
    def _ConditionalOverlay(self,locus, datum, srcCompL, dstComp):
        dstNull = self.process.dstLayerD[locus][datum][dstComp].comp.cellnull
        #Create single dst array
        if self.process.params.method == 'avg': 
            dstBAND = np.zeros_like( self.process.srcLayerD[locus][datum][srcCompL[0]].layer.NPBAND )
            #Loop over the scr compositions
            for c in srcCompL:
                #Set null to nan
                BAND = self.process.srcLayerD[locus][datum][c].layer.NPBAND
                cellnull = self.process.srcLayerD[locus][datum][c].comp.cellnull
                BAND[BAND == cellnull] = np.nan
                #Copy array
                srcBAND = BAND*1
                #Set nan to zero
                srcBAND[np.isnan(srcBAND)] = 0
                dstBAND += srcBAND   
            #Set to nan where both input layers are non
            dstBAND[ (np.isnan( self.process.srcLayerD[locus][datum][srcCompL[0]].layer.NPBAND ) ) & (np.isnan( self.process.srcLayerD[locus][datum][srcCompL[1]].layer.NPBAND ) ) ] = dstNull
            #Set to average where both input bands are valid
            dstBAND[ (~np.isnan( self.process.srcLayerD[locus][datum][srcCompL[0]].layer.NPBAND ) ) & (~np.isnan( self.process.srcLayerD[locus][datum][srcCompL[1]].layer.NPBAND ) ) ] /= 2
        
        elif self.process.params.method == 'min':
            dstBAND = np.full_like( self.process.srcLayerD[locus][datum][srcCompL[0]].layer.NPBAND, np.nan )
            #Loop over the scr compositions
            for c in srcCompL:
                #Set null to nan
                BAND = self.process.srcLayerD[locus][datum][c].layer.NPBAND
                cellnull = self.process.srcLayerD[locus][datum][c].comp.cellnull
                BAND[BAND == cellnull] = np.nan
                #fmin
                dstBAND = np.fmin(dstBAND, BAND )
            #Set to nan where both input layers are non
            dstBAND[ (np.isnan( self.process.srcLayerD[locus][datum][srcCompL[0]].layer.NPBAND ) ) & (np.isnan( self.process.srcLayerD[locus][datum][srcCompL[1]].layer.NPBAND ) ) ] = dstNull
             
        elif self.process.params.method == 'max':
            dstBAND = np.full_like( self.process.srcLayerD[locus][datum][srcCompL[0]].layer.NPBAND, np.nan )
            #Loop over the scr compositions
            for c in srcCompL:
                #Set null to nan
                BAND = self.process.srcLayerD[locus][datum][c].layer.NPBAND
                cellnull = self.process.srcLayerD[locus][datum][c].comp.cellnull
                BAND[BAND == cellnull] = np.nan
                #fmax
                dstBAND = np.fmax(dstBAND, BAND )
            #Set to nan where both input layers are non
            dstBAND[ (np.isnan( self.process.srcLayerD[locus][datum][srcCompL[0]].layer.NPBAND ) ) & (np.isnan( self.process.srcLayerD[locus][datum][srcCompL[1]].layer.NPBAND ) ) ] = dstNull

        else:   
            ERRORIGEN
       
        #Create the dst layer
        self.process.dstLayerD[locus][datum][dstComp].layer = lambda:None
        #Set the np array as the band
        self.process.dstLayerD[locus][datum][dstComp].layer.NPBAND = dstBAND
        #copy the geoformat from the src layer
        self.process.dstLayerD[locus][datum][dstComp].CopyGeoformatFromSrcLayer(self.process.srcLayerD[locus][datum][c].layer)
        #write the results
        self.process.dstLayerD[locus][datum][dstComp].CreateDSWriteRasterArray()
        #Register the layer
        self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstComp], self.process.overwrite, self.process.delete)

    def _LinearTransformMODISSingleTile(self, locus, datum, srcLayerD, dstLayerD):
        #create empty 2D array for the output
        #dstNull = self.process.dstLayerD[locus][datum][dstComp].comp.cellnull
        srcKeys = list(srcLayerD.keys())

        for dst in dstLayerD:
            #Create the dst layer
            dstLayerD[dst].layer = lambda:None
            #Set the np array as the band
            dstLayerD[dst].layer.NPBAND = np.zeros(srcLayerD[srcKeys[0]].layer.NPBAND.shape) 

        #then run
        for dst in dstLayerD:
            for src in srcLayerD:
                xid = '%s%s' %(srcLayerD[src].comp.id, 
                               dstLayerD[dst].comp.id) 
                scalefacD = getattr(self.process.proc.transformscale, xid)
                scalefac = scalefacD['scalefac']
                offsetD = getattr(self.process.proc.transformoffset, srcLayerD[src].comp.id)
                offset = offsetD['offset']
                #ImageTransform is a numba JIT function
                dstLayerD[dst].layer.NPBAND = ImageTransform(dstLayerD[dst].layer.NPBAND,srcLayerD[src].layer.NPBAND,offset,scalefac)
                #dstLayerD[dst].layer.NPBAND += (srcLayerD[src].layer.NPBAND + offset)*scalefac
        #self._SetMask(locus, datum, srcLayerD, dstLayerD)
        MultibandMasking(locus, datum, srcLayerD, dstLayerD)
        
        for dst in dstLayerD:
            #copy the geoformat from the src layer
            dstLayerD[dst].CopyGeoformatFromSrcLayer(srcLayerD[srcKeys[0]].layer)
            #write the results
            dstLayerD[dst].CreateDSWriteRasterArray()
            #Register the layer
            self.session._InsertLayer(dstLayerD[dst], self.process.overwrite, self.process.delete)
               
    def _fgbgmodisSingleTile(self, locus, datum, srcLayerD, dstLayerD):
        '''
        '''
        #srcKeys = list(srcLayerD.keys())
        X = srcLayerD[self.xb].layer.NPBAND
        Y = srcLayerD[self.yb].layer.NPBAND
        if self.fg:
            dstLayerD['fg'].layer = lambda:None
            #Set the np array as the band
            #self.process.dstLayerD[locus][datum][dstComp].layer.NPBAND = dstBAND
            #dstLayerD[dst].layer.NPBAND = np.zeros(srcLayerD[srcKeys[0]].layer.NPBAND.shape) 
            dstLayerD['fg'].layer.NPBAND = ImageFgBg(self.process.params.rescalefac, 
                    self.sinrang, self.cosrang, X, Y, self.process.params.intercept, self.process.params.calibfac)
            #dstLayerD['fg'].layer.NPBAND = self.process.params.rescalefac * ((self.sinrang*(x+y-self.process.params.intercept) + self.cosrang*(-x+y-self.process.params.intercept)) / 
            #                 (self.sinrang*(x-y+self.process.params.intercept) + self.cosrang*( x+y-self.process.params.intercept) + self.process.params.calibfac ))
        #FG =  5942*( ( self.sinrang*(x+y+2080) + self.cosrang*(-x+y+2080) ) / ( self.sinrang*(x - y - 2080)+self.cosrang*(x+y+2080) + 7000 ) )
        #twi = 5942*( ( _sinrang*(x+y+2080) + _cosrang*(-x+y+2080) ) / ( _sinrang*(x - y - 2080)+_cosrang*(x+y+2080) + 7000 ) )
        if self.bg:
            dstLayerD['bg'].layer = lambda:None
            
            #BG = self.process.rescalefac * ((self.sinrang*(x+y+self.process.intercept) + self.cosrang*(-x+y+self.process.intercept)) / 
            #                     (self.sinrang*(x-y-self.process.intercept) + self.cosrang*( x+y+self.process.intercept) + self.process.calibfac ))
        MultibandMasking(locus, datum, srcLayerD, dstLayerD)
        #self._SetMask(locus, datum, srcLayerD, dstLayerD)

        for dst in dstLayerD:
            #copy the geoformat from the src layer
            dstLayerD[dst].CopyGeoformatFromSrcLayer(srcLayerD[self.xb].layer)
            #write the results
            dstLayerD[dst].CreateDSWriteRasterArray()
            #Register the layer
            self.session._InsertLayer(dstLayerD[dst], self.process.overwrite, self.process.delete)
                            
class GDALstuff:
    def __init__(self, srcFPN, dstFPN, params):
        """The constructor is .""" 
        self.srcFPN = srcFPN
        self.dstFPN = dstFPN
        self.params = params
  
    def _GetRasterMeta(self):
        self.spatialRef, self.metadata = mj_gis.GetRasterMetaData(self.srcFPN)
        
    def _GetVectorProjection(self):
        self.spatialRef = mj_gis.GetVectorProjection(self.srcFPN)
        
    def MosaicRaster(self):
        '''
        '''
        gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdalbuildvrt ' 
        gdalcmd = '%(cmd)s -input_file_list %(lisfpn)s %(tar)s ' %{'cmd':gdalcmd, 'lisfpn':self.srcFPN, 'tar':self.dstFPN}
        os.system(gdalcmd)
        
    def SetClipBoxLLminmax(self,regExt):
        self.ulx = min(regExt[0:2])
        self.uly = max(regExt[2:4])
        self.lrx = max(regExt[4:6])
        self.lry = min(regExt[6:8])
             
    def SetClipBox(self,ulx,uly,lrx,lry):
        self.ulx = ulx
        self.uly = uly
        self.lrx = lrx
        self.lry = lry
        
    def SetTarEPSGCoordinates(self,srcEPSG,tarEPSG):
        srcProj = mj_gis.MjProj()
        srcProj.SetFromEPSG(srcEPSG)
        #if not roiLayer.Exists() or self.process.overwrite: #or overwrite
        #    mj_gis.CreateESRIPolygonPtL(roiLayer.FPN, roiLayer.fieldDefL, region.BoundsPtL, projection.proj_cs, region.regionid)
        #Get the bounds in the original projection
        #boundsD = mj_gis.GetFeatureBounds(roiLayer.FPN,'REGIONID')
        boundsPtL = [(self.ulx, self.uly),(self.lrx,self.uly), (self.lrx,self.lry), (self.ulx,self.lry)]
        #Set lonlat projection
        tarProj = mj_gis.MjProj()
        tarProj.SetFromEPSG(tarEPSG)
        #Get the corners in lonlat
        #print 'srcProj.proj_cs',srcProj.proj_cs
        #print 'tarProj.proj_cs',tarProj.proj_cs
        coords = mj_gis.ReprojectBoundsToAny(boundsPtL,srcProj.proj_cs, tarProj.proj_cs)
        #print coords
        self.ulx = min(coords['ulx'],coords['llx'])
        self.uly = max(coords['uly'],coords['ury'])
        self.lrx = max(coords['urx'],coords['lrx'])
        self.lry = min(coords['lly'],coords['lry'])
        #print self.ulx,self.uly,self.lrx,self.lry
            
    def ClipRasterOld(self,geoFormatD):
    
        layerIn_proj = mj_gis.MjProj()
        layerIn_proj.SetFromWKT(geoFormatD['projection'])
        layerIn_proj.ReadSpatialRef()
        #BALLE
        #if layerIn_proj.epsg != tar_epsg:
        #print  bounds_epsg, tar_epsg
        if self.params.t_epsg == 0 or self.params.t_epsg == layerIn_proj.epsg:
            print ('setting dst epsg to src epsg')
            dst_epsg = layerIn_proj.epsg
        elif self.params.t_epsg != self.params.bounds_epsg:
            self.SetTarEPSGCoordinates(self.params.bounds_epsg, self.params.t_epsg)

        if self.params.t_epsg == layerIn_proj.epsg:
            gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdal_translate'
            gdalcmd = '%(cmd)s -projwin %(ulx)f %(uly)f %(lrx)f %(lry)f' %{'cmd':gdalcmd, 'ulx':self.ulx, 'uly':self.uly, 'lrx':self.lrx, 'lry':self.lry}
            gdalcmd = '%(cmd)s %(src)s %(dst)s' %{'cmd':gdalcmd, 'src':self.layerInFPN, 'dst':self.layerOutFPN}
            print (gdalcmd)
            BALLE
            os.system(gdalcmd)                  
                #print '    Registering',self.layerOutD[key].comp.band,self.layerOutD[key].FPN
        else:
            layerOut_proj = mj_gis.MjProj()
            layerOut_proj.SetFromEPSG(self.params.t_epsg)
            layerOut_proj.SetProj4()
            gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdalwarp -t_srs "%(cs)s"' %{'cs':layerOut_proj.proj4}
            #if self.aD['celltype'] != 'auto':
            #    gdalcmd = ' %(cmd)s -ot %(ot)s ' %{'cmd':gdalcmd, 'ot': self.aD['celltype']}

            if self.params.xres == self.params.yres == 0:
                pass 
            else:
                gdalcmd = ' %(cmd)s -tr %(tr_x)1.9f %(tr_y)1.9f ' %{'cmd':gdalcmd, 'tr_x': self.params.tr_xres, 'tr_y': self.params.tr_yres}
            gdalcmd = ' %(cmd)s -te %(ulx)s %(lry)s %(lrx)s %(uly)s ' %{'cmd':gdalcmd, 'ulx': self.ulx,'uly': self.uly, 'lrx': self.lrx, 'lry': self.lry}
            if self.params.resample:
                gdalcmd = ' %(cmd)s -r %(resample)s ' %{'cmd':gdalcmd, 'resample': self.params.resample}
    
            #linestr = ' %(s1)s -dstnodata %(nodata)d -r %(resample)s -ot %(celltype)s -of ERS' %{'s1':linestr, 'nodata': self.nodata, 'resample': self.resample, 'celltype': self.cellType}
            gdalcmd = '%(cmd)s %(src)s %(dst)s' %{'cmd':gdalcmd, 'src':self.srcFPN, 'dst':self.dstFPN}
            print (gdalcmd)
            os.system(gdalcmd)
            
    def ClipRaster(self):
        '''
        '''
        #Get the projection and other metadata for the source layer (mosaic)
        self._GetRasterMeta()
        
        if self.params.t_epsg == 0 or self.params.t_epsg == self.spatialRef.epsg:
            print ('setting dst epsg to src epsg')
            dst_epsg = self.spatialRef.epsg
        elif self.params.t_epsg != self.spatialRef.epsg:
            if self.params.t_epsg != 4326:
                #not lonlat
                self.SetTarEPSGCoordinates(self.params.bounds_epsg, self.params.t_epsg)
            else:
                dstEPSG = 4326

        if self.params.t_epsg == self.spatialRef.epsg:
            #No reprojection, just cut
            gdalcmd = os.path.join(GDALpath,'gdal_translate')
            gdalcmd = '%(cmd)s -projwin %(ulx)f %(uly)f %(lrx)f %(lry)f' %{'cmd':gdalcmd, 'ulx':self.ulx, 'uly':self.uly, 'lrx':self.lrx, 'lry':self.lry}
            gdalcmd = '%(cmd)s %(src)s %(dst)s' %{'cmd':gdalcmd, 'src':self.layerInFPN, 'dst':self.layerOutFPN}
            print (gdalcmd)
            BALLE
            os.system(gdalcmd)                  
            #print '    Registering',self.layerOutD[key].comp.band,self.layerOutD[key].FPN
        else:
            dstProj = mj_gis.MjProj()
            dstProj.SetFromEPSG(dstEPSG)
            dstProj.SetProj4()
            gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdalwarp -t_srs "%(cs)s"' %{'cs':dstProj.proj4}    
            if self.params.xres == self.params.yres == 0:
                pass 
            else:
                gdalcmd = ' %(cmd)s -tr %(tr_x)1.9f %(tr_y)1.9f ' %{'cmd':gdalcmd, 'tr_x': self.params.xres, 'tr_y': self.params.yres}
            gdalcmd = ' %(cmd)s -te %(ulx)s %(lry)s %(lrx)s %(uly)s ' %{'cmd':gdalcmd, 'ulx': self.ulx,'uly': self.uly, 'lrx': self.lrx, 'lry': self.lry}
            if self.params.resample:
                gdalcmd = ' %(cmd)s -r %(resample)s ' %{'cmd':gdalcmd, 'resample': self.params.resample}

            gdalcmd = '%(cmd)s %(src)s %(dst)s' %{'cmd':gdalcmd, 'src':self.srcFPN, 'dst':self.dstFPN}
            os.system(gdalcmd)

    def ClipVector(self):
        '''
        '''
        #Get the projection and other metadata for the source layer (mosaic)
        self._GetVectorProjection()

        if self.params.t_epsg == 0 or self.params.t_epsg == self.spatialRef.epsg:
            dstEPSG = self.spatialRef.epsg
        elif self.params.t_epsg != self.spatialRef.epsg:
            if self.params.t_epsg != 4326:
                #not lonlat
                self.SetTarEPSGCoordinates(self.params.bounds_epsg, self.params.t_epsg)
            else:
                dstEPSG = 4326

        if self.params.t_epsg == self.spatialRef.epsg:
            #No reprojection, just cut
            gdalcmd = os.path.join(GDALpath,'ogr2ogr')
            gdalcmd = '%(cmd)s -overwrite -clipdst %(ulx)f %(lry)f %(lrx)f %(uly)f ' %{'cmd':gdalcmd, 'ulx':self.ulx, 'uly':self.uly, 'lrx':self.lrx, 'lry':self.lry}
            if self.params.tolerance: 
                gdalcmd += '-simplify %s ' %(self.params.tolerance)
            
            gdalcmd = '%(cmd)s %(dst)s %(src)s ' %{'cmd':gdalcmd, 'src':self.srcFPN, 'dst':self.dstFPN}

            os.system(gdalcmd)                  
            #print '    Registering',self.layerOutD[key].comp.band,self.layerOutD[key].FPN
        else:
            dstProj = mj_gis.MjProj()
            dstProj.SetFromEPSG(dstEPSG)
            dstProj.SetProj4()
            gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/ogr2ogr -t_srs "%(cs)s"' %{'cs':dstProj.proj4}    
            if self.params.xres == self.params.yres == 0:
                pass 
            else:
                gdalcmd = ' %(cmd)s -tr %(tr_x)1.9f %(tr_y)1.9f ' %{'cmd':gdalcmd, 'tr_x': self.params.xres, 'tr_y': self.params.yres}
            gdalcmd = ' %(cmd)s -te %(ulx)s %(lry)s %(lrx)s %(uly)s ' %{'cmd':gdalcmd, 'ulx': self.ulx,'uly': self.uly, 'lrx': self.lrx, 'lry': self.lry}
            if self.params.resample:
                gdalcmd = ' %(cmd)s -r %(resample)s ' %{'cmd':gdalcmd, 'resample': self.params.resample}

            gdalcmd = '%(cmd)s %(src)s %(dst)s' %{'cmd':gdalcmd, 'src':self.srcFPN, 'dst':self.dstFPN}
            print (gdalcmd)

            os.system(gdalcmd)
            
    def SetTargetProj(self,tarEPSG):
        self.tarProj = mj_gis.MjProj()
        self.tarProj.SetFromEPSG(tarEPSG)
        
    def WarpRaster(self, asscript=False):
        self.tarProj.SetProj4()
        gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdalwarp -t_srs "%(cs)s"' %{'cs':self.tarProj.proj4}
        if self.params.celltype != 'auto':
            gdalcmd = ' %(cmd)s -ot %(ot)s ' %{'cmd':gdalcmd, 'ot': self.params.celltype}
 
        if self.params.xres == self.params.yres == 0:
            pass 
        else:
            gdalcmd = '%(cmd)s -tr %(trx)1.9f %(try)1.9f ' %{'cmd':gdalcmd, 'trx': self.params.xres, 'try':self.params.yres}
        gdalcmd = '%(cmd)s -te %(ulx)s %(lry)s %(lrx)s %(uly)s ' %{'cmd':gdalcmd, 'ulx': self.ulx,'uly': self.uly, 'lrx': self.lrx, 'lry': self.lry}
        if self.params.resample:
            gdalcmd = '%(cmd)s -r %(resample)s ' %{'cmd':gdalcmd, 'resample': self.params.resample}
   
        #linestr = ' %(s1)s -dstnodata %(nodata)d -r %(resample)s -ot %(celltype)s -of ERS' %{'s1':linestr, 'nodata': self.nodata, 'resample': self.resample, 'celltype': self.cellType}
        gdalcmd = '%(cmd)s %(src)s %(dst)s;' %{'cmd':gdalcmd, 'src':self.srcFPN, 'dst':self.dstFPN}
        if asscript:
            return gdalcmd
        print (gdalcmd)
        os.system(gdalcmd)
        return False
    
    def TransformOutSize(self,xsize,ysize,resample,ot):
        gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_translate' 
        gdalcmd = '%(cmd)s -ot %(ot)s -outsize %(xsize)s %(ysize)d -r %(r)s' %{'cmd':gdalcmd, 'ot':ot, 'xsize':xsize, 'ysize':ysize, 'r':resample}
        gdalcmd = '%(cmd)s %(src)s %(dst)s;' %{'cmd':gdalcmd, 'src':self.srcFPN, 'dst':self.dstFPN}
        os.system(gdalcmd)
        
    def TransformOT(self,ot):
        gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_translate' 
        gdalcmd = '%(cmd)s -ot %(ot)s' %{'cmd':gdalcmd, 'ot':ot}
        gdalcmd = '%(cmd)s %(src)s %(dst)s;' %{'cmd':gdalcmd, 'src':self.srcFPN, 'dst':self.dstFPN}
        os.system(gdalcmd)
        
    def TransformOF(self, of):
        '''
        '''
 
        gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_translate'
        

        gdalcmd = '%(cmd)s -of %(of)s -B 1 -B 2 -B 3' %{'cmd':gdalcmd, 'of':of}

        gdalcmd = '%(cmd)s %(src)s %(dst)s;' %{'cmd':gdalcmd, 'src':self.srcFPN, 'dst':self.dstFPN}
        print (gdalcmd)
        os.system(gdalcmd)
        
    def ResampleRaster(self,asscript=False):
        gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_translate' 
        gdalcmd = '%(cmd)s -tr %(xres)s %(yres)s -r %(r)s' %{'cmd':gdalcmd,  'xres':self.params.xres, 'yres':self.params.yres, 'r':self.params.resample}
        gdalcmd = '%(cmd)s %(src)s %(dst)s;' %{'cmd':gdalcmd, 'src':self.srcFPN, 'dst':self.dstFPN}
        if asscript:
            return gdalcmd
        os.system(gdalcmd)
        return False
   
    def TPI(self):
        gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdaldem TPI %(src)s %(dst)s' %{'src':self.srcFPN, 'dst':self.dstFPN}
        if self.params.compute_edges:
            gdalcmd = '%(cmd)s -compute_edges;' %{'cmd':gdalcmd}
        os.system(gdalcmd)
        
    def TRI(self):
        gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdaldem TRI %(src)s %(dst)s' %{'src':self.srcFPN, 'dst':self.dstFPN}
        if self.params.compute_edges:
            gdalcmd = '%(cmd)s -compute_edges;' %{'cmd':gdalcmd}
        os.system(gdalcmd)
        
    def Roughness(self):
        gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdaldem roughness %(src)s %(dst)s' %{'src':self.srcFPN, 'dst':self.dstFPN}
        if self.params.compute_edges:
            gdalcmd = '%(cmd)s -compute_edges;' %{'cmd':gdalcmd}
        os.system(gdalcmd)

    def HillShade(self):
        gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdaldem hillshade %(src)s %(dst)s' %{'src':self.srcFPN, 'dst':self.dstFPN}
        if self.params.compute_edges:
            gdalcmd = '%(cmd)s -compute_edges' %{'cmd':gdalcmd}
        if self.params.zfactor:
            gdalcmd = '%(cmd)s -z %(zf)f' %{'cmd':gdalcmd, 'zf':self.params.zfactor}
        if self.params.scale:
            gdalcmd = '%(cmd)s -s %(sc)f' %{'cmd':gdalcmd, 'sc':self.params.scale}
        gdalcmd = '%(cmd)s -az %(az)f' %{'cmd':gdalcmd, 'az':self.params.azimuth}
        gdalcmd = '%(cmd)s -alt %(alt)f' %{'cmd':gdalcmd, 'alt':self.params.altitude}
        if self.params.combined:
            gdalcmd = '%(cmd)s -combined' %{'cmd':gdalcmd}
        if self.params.combined:
            gdalcmd = '%(cmd)s -multidirectional' %{'cmd':gdalcmd}
        os.system(gdalcmd)
        
    def Slope(self):
        gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdaldem slope %(src)s %(dst)s' %{'src':self.srcFPN, 'dst':self.dstFPN}
        if self.params.compute_edges:
            gdalcmd = '%(cmd)s -compute_edges' %{'cmd':gdalcmd}
        if self.params.ZevenbergenThorne:
            gdalcmd = '%(cmd)s -alg ZevenbergenThorne' %{'cmd':gdalcmd}   
        if self.params.percent:
            gdalcmd = '%(cmd)s -p' %{'cmd':gdalcmd}
        if self.params.scale:
            gdalcmd = '%(cmd)s -s %(sc)f' %{'cmd':gdalcmd, 'sc':self.params.scale}
        os.system(gdalcmd)
        
    def Aspect(self):
        gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdaldem aspect %(src)s %(dst)s' %{'src':self.srcFPN, 'dst':self.dstFPN}
        if self.params.compute_edges:
            gdalcmd = '%(cmd)s -compute_edges' %{'cmd':gdalcmd}
        if self.params.ZevenbergenThorne:
            gdalcmd = '%(cmd)s -alg ZevenbergenThorne' %{'cmd':gdalcmd}
        if self.params.trigonometric:
            gdalcmd = '%(cmd)s -trigonometric' %{'cmd':gdalcmd}
        if self.params.zero_for_flat:
            gdalcmd = '%(cmd)s -zero_for_flat' %{'cmd':gdalcmd}

        os.system(gdalcmd)
        
    def Delete(self):
        gdalcmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdalmanage delete %(dst)s' %{'dst':self.dstFPN}
        os.system(gdalcmd)

        