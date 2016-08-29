from pyfem.materials.BaseMaterial import BaseMaterial
from numpy import zeros, dot

class PlaneStrain( BaseMaterial ):

  def __init__ ( self, props ):

    #Call the BaseMaterial constructor
    BaseMaterial.__init__( self, props )

    #Create the hookean matrix
    self.H = zeros( (3,3) )

    self.H[0,0] = self.E*(1.-self.nu)/((1+self.nu)*(1.-2.*self.nu))
    self.H[0,1] = self.H[0,0]*self.nu/(1-self.nu)
    self.H[1,0] = self.H[0,1]
    self.H[1,1] = self.H[0,0]
    self.H[2,2] = self.H[0,0]*0.5*(1.-2.*self.nu)/(1.-self.nu)
    
    self.setHistoryParameter("stresses",zeros(4))
    self.commitHistory()

  def getStress( self, deformation ):
    sigma0 = self.getHistoryParameter("stresses")
    sigma = zeros(4);
    sigma[[0,1,3]] = sigma0[[0,1,3]]+dot( self.H, deformation.dstrain )
    sigma[2] = self.nu*(sigma[0]+sigma[1])
    self.setHistoryParameter("stresses",sigma)

    return sigma[[0,1,3]], self.H

  def getTangent( self ):
  
    return self.H

