# -*- coding: utf-8 -*-

from pyfem.materials.BaseMaterial import BaseMaterial
from numpy import zeros, dot

class Isotropic( BaseMaterial ):

  def __init__ ( self, props ):

    #Call the BaseMaterial constructor
    BaseMaterial.__init__( self, props )

    #Create the hookean matrix
    self.H = zeros( (6,6) )

    fac = 1.0 / (2.0 * self.nu * self.nu + self.nu - 1.0 );
  
    self.H[0,0] = fac * self.E * ( self.nu - 1.0 );
    self.H[0,1] = -1.0 * fac * self.E * self.nu;
    self.H[0,2] = self.H[0,1];
    self.H[1,0] = self.H[0,1];                                  
    self.H[1,1] = self.H[0,0];                                  
    self.H[1,2] = self.H[0,1];                                  
    self.H[2,0] = self.H[0,1];                                  
    self.H[2,1] = self.H[0,1];                                  
    self.H[2,2] = self.H[0,0];                                  
    self.H[3,3] = self.E / ( 2.0 + 2.0 * self.nu );          
    self.H[4,4] = self.H[3,3];                                  
    self.H[5,5] = self.H[3,3];
    
    self.setHistoryParameter("stresses",zeros(6))
    self.commitHistory()

  def getStress( self, deformation ):
#    if hasattr(self,'history'):
#          stresses = self.getHistoryParameter("stresses")
#          print "stresses = ",stresses  
    sigma = dot( self.H, deformation.strain )
    self.setHistoryParameter("stresses",sigma)
    return sigma, self.H

  def getTangent( self ):
  
    return self.H