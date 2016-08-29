# -*- coding: utf-8 -*-

class MaterialManager ( list ):

  def __init__ ( self, matProps ):

    matType = matProps.type

    cmdStr = 'from pyfem.materials.' + matType + ' import ' + matType + ' as material'

    exec cmdStr
    
    self.mat = material( matProps )
    self.iIter = -1

  def reset( self ):

    self.iIter  = -1

  def getStress ( self, kinematic , iSam = -1 ):

    if iSam == -1:
      self.iIter += 1
      iSam = self.iIter

    self.mat.setIter( iSam )
     
    return self.mat.getStress( kinematic )
    
  def getHistoryParameter ( self , name , iSam = -1 ):
    
    if iSam == -1:
      self.iIter += 1
      iSam = self.iIter
      
    self.mat.setIter( iSam )
    
    return self.mat.getHistoryParameter(name)
    

  def commitHistory( self ):
    self.mat.commitHistory()
    
  def commitStress( self ):
    self.mat.commitStress()
    
  def clearHistory( self ):
    self.mat.clearHistory()

