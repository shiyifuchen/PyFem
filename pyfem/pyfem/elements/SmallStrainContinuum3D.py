# -*- coding: utf-8 -*-

from .Element import Element
from pyfem.util.shapeFunctions  import getElemShapeData
from pyfem.util.kinematics      import Kinematics
from numpy import zeros, dot, outer, ones , eye , ndarray
from scipy.linalg import inv

class SmallStrainContinuum3D( Element ):

  #dofs per element
  dofTypes = [ 'u' , 'v' , 'w' ]
  
  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )

  def __type__ ( self ):
    return name

#------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):

    sData = getElemShapeData( elemdat.coords )

    kin = Kinematics(3,6)
    
#    elemdat.outlabel.append("stresses")
#    elemdat.outdata  = zeros( shape=(len(elemdat.nodes),6) )
#    elemdat.extrapolation = zeros(shape=(len(elemdat.nodes),len(elemdat.nodes)))

    for i,iData in enumerate(sData):
      
      b = self.getBmatrix( iData.dhdx )
#      elemdat.extrapolation[i,:]=iData.h

      kin.strain  = dot ( b , elemdat.state )
      kin.dstrain = dot ( b , elemdat.Dstate )
      
      sigma,tang = self.mat.getStress( kin )

      elemdat.stiff += dot ( b.transpose() , dot ( tang , b ) ) * iData.weight
      elemdat.fint  += dot ( b.transpose() , sigma ) * iData.weight

      #elemdat.outdata += outer( ones(len(self)), sigma )
#      elemdat.outdata += outer( eye(len(self),1,-i), sigma )
      
    #elemdat.outdata *= 1.0 / len(sData)
#    elemdat.outdata = dot(inv(elemdat.extrapolation),elemdat.outdata)
    #print elemdat.outdata

  
     
#-------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

#    elemdat.outlabel.append("stresses")
#    elemdat.outdata  = zeros( shape=(len(elemdat.nodes),6) )
#    elemdat.extrapolation = zeros(shape=(len(elemdat.nodes),len(elemdat.nodes)))

    kin = Kinematics(3,6)

    for i,iData in enumerate(sData):
      b = self.getBmatrix( iData.dhdx )
#      elemdat.extrapolation[i,:]=iData.h

      kin.strain  = dot ( b , elemdat.state )
      kin.dstrain = dot ( b , elemdat.Dstate )

      sigma,tang = self.mat.getStress( kin )

      elemdat.fint    += dot ( b.transpose() , sigma ) * iData.weight
      #elemdat.outdata += outer( ones(len(self)), sigma )
#      elemdat.outdata += outer( eye(len(self),1,-i), sigma )
      
    #elemdat.outdata *= 1.0 / len(sData)
#    elemdat.outdata = dot(inv(elemdat.extrapolation),elemdat.outdata)


#----------------------------------------------------------------------
    
  def getMassMatrix ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )
    
    if not hasattr(self,"rho"):
      setattr(self,"rho",0.0);

    rho = self.rho * eye(3)

    for iData in sData:
      N  = self.getNmatrix( iData.h )
      Nt = N.transpose()

      elemdat.mass += dot ( Nt , dot( rho , N ) ) * iData.weight
      
    elemdat.lumped = sum(elemdat.mass)
   
#--------------------------------------------------------------------------

  def getBmatrix( self , dphi ):

    b = zeros( shape=( 6 , self.dofCount() ) )

    for i,dp in enumerate(dphi):
      b[0,i*3  ] = dp[0]
      b[1,i*3+1] = dp[1]
      b[2,i*3+2] = dp[2]
      b[3,i*3  ] = dp[1]
      b[3,i*3+1] = dp[0]
      b[4,i*3+1] = dp[2]
      b[4,i*3+2] = dp[1]
      b[5,i*3  ] = dp[2]
      b[5,i*3+2] = dp[0]
   
    return b

#------------------------------------------------------------------------------

  def getNmatrix( self , h ):

    N = zeros( shape=( 3 , 3*len(h) ) )

    for i,a in enumerate( h ):
      N[0,3*i  ] = a
      N[1,3*i+1] = a
      N[1,3*i+2] = a
   
    return N
    
    
  def getOutputData( self , elemdat ):
    
    sData = getElemShapeData( elemdat.coords )
    extrapolation = zeros(shape=(len(elemdat.nodes),len(elemdat.nodes)))
    elemdat.outdata = []
    for label in elemdat.outlabel:
      column = 1
      if isinstance(self.mat.getHistoryParameter(label,0),ndarray):
        column = len(self.mat.getHistoryParameter(label,0))
      elemdat.outdata.append(zeros((len(elemdat.nodes),column)))
    
    for i,iData in enumerate(sData):
      extrapolation[i,:]=iData.h
      for j,label in enumerate(elemdat.outlabel):
        elemdat.outdata[j]+=outer(eye(len(self),1,-i),\
          self.mat.getHistoryParameter(label,i))
          
    for i,label in enumerate(elemdat.outlabel):
      if label == "equivalent_strain":
        elemdat.outdata[i]=elemdat.outdata[i][[0,4,6,2,1,5,7,3]]
        continue
#      elemdat.outdata[i] = elemdat.outdata[i][[0,4,6,2,1,5,7,3]]
      elemdat.outdata[i] = dot(inv(extrapolation),elemdat.outdata[i])
