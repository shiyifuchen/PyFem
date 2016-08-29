# -*- coding: utf-8 -*-
from numpy import outer, ones, zeros
from pyfem.materials.MaterialManager import MaterialManager

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class Element ( list ):

  dofTypes = []

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def __init__ ( self, elnodes , props ):
    list.__init__( self, elnodes )

    self.history = {}
    self.current = {}
    
    #props是与单元相关的所有属性
    for name,val in props:
      if name is "material":
        self.matProps = val
        self.mat = MaterialManager( self.matProps )
      else:
        setattr( self, name, val )

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def dofCount ( self ):

    return len( self ) * len( self.dofTypes )

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getNodes ( self ):
    return self

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getType ( self ):
    return self.elemType

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def appendNodalOutput ( self, outputName, globdat, outmat, outw = None ):

    if outw == None:
      outw = ones( outmat.shape[0] )

    if not hasattr( globdat, outputName ):

      globdat.outputNames.append( outputName )

      setattr( globdat, outputName, zeros( ( len(globdat.nodes), outmat.shape[1] ) ) )
      setattr( globdat, outputName + 'Weights', zeros( len(globdat.nodes) ) )

    outMat     = getattr( globdat, outputName )
    outWeights = getattr( globdat, outputName + 'Weights' )

    if outmat.shape[1] != outMat.shape[1] or outmat.shape[0] != len(self):
      raise RuntimeError("Appended output vector has incorrect size.")

    indi = globdat.nodes.getIndices( self )

    #outMat[ indi ]     += outer( outw, ones(outmat.shape[1]) ) * outmat
#    if self.name =="Elem1":
#      return

    outMat[ indi ]     += outmat
    outWeights[ indi ] += outw

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def setHistoryParameter ( self, name, val ):
    self.current[name] = val

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getHistoryParameter ( self, name ):
    return self.history[name]

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def commitHistory ( self ):
    
#    print "Element CommitHistory ......"
    self.history = self.current.copy()
    self.current = {}

    if hasattr( self , "mat" ):
      self.mat.commitHistory()
      
  def commitStress( self ):
    
    self.history = self.current.copy()
    self.current = {}

    if hasattr( self , "mat" ):
      self.mat.commitStress()    
