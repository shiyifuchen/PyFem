# -*- coding: utf-8 -*-

from pyfem.materials.BaseMaterial import BaseMaterial
from numpy import zeros,dot,outer
from numpy.linalg import matrix_rank,inv
from math import sqrt
from pyfem.util.tensor import decomp,norm,unitFlow,getII,getI_bar

class J2Plastic(BaseMaterial):
  def __init__(self,props):
    BaseMaterial.__init__(self,props)
    
    self.K = self.E/3.0/(1-2*self.nu)
    self.G = self.E/2.0/(1+self.nu)
    self.C = self.K*getII(3)+2*self.G*getI_bar(3)
    self.Cep = self.C.copy()
    self.setHistoryParameter("stresses",zeros(6))
    self.setHistoryParameter("equivalent_strain",0.0)
    self.commitHistory()
    
  def getStress( self, deformation ):
    
    sigma0=self.getHistoryParameter("stresses")
    eqStrain=self.getHistoryParameter("equivalent_strain")
      
    dstrain = deformation.dstrain.copy()
    dstrain[3:] *= 0.5
    sph_sigma0,dev_sigma0 = decomp(sigma0)
    sph_dstrain,dev_dstrain = decomp(dstrain)
    
    deltaGamma = 0.0
    alpha = eqStrain
    S_trial = dev_sigma0 + 2 * self.G * dev_dstrain
    
    g = norm(S_trial) - sqrt(2./3) * self.__getYieldStress(alpha)\
        - 2 * self.G * deltaGamma
    
    if g<=0:
      sigma = sigma0 + dot(self.C, deformation.dstrain )
      self.setHistoryParameter("stresses",sigma)
      self.setHistoryParameter("equivalent_strain",eqStrain)
      return sigma, self.C
    
    k = 0
    while (abs(g) > 1e-8):
      k = k + 1
      if k > 10:
        raise RuntimeError('Newton-Raphson iterations did not converge!')
      Dg = -2 * self.G - 2./3 * self.H
      
      deltaGamma -= g/Dg
      alpha = eqStrain + sqrt(2./3) * deltaGamma
      g = norm(S_trial) - sqrt(2./3) * self.__getYieldStress(alpha) \
        - 2 * self.G * deltaGamma
      
    n = unitFlow(S_trial)
    eqStrain = alpha
    sigma =  sph_sigma0 + 3*self.K*sph_dstrain + S_trial \
        -2*self.G*deltaGamma*n
    self.Cep = self.C - 6*self.G*self.G/(3*self.G+self.H)*outer(n,n)
    
    self.Cep = self.Cep + 4*self.G*self.G*deltaGamma/norm(S_trial)*\
        (outer(n,n)-getI_bar(3))
        
    print matrix_rank(self.Cep-self.C)
    print dot(inv(self.C-self.Cep),self.C-self.Cep)
    
    
    self.setHistoryParameter("stresses",sigma)
    self.setHistoryParameter("equivalent_strain",eqStrain)
    return sigma, self.Cep

  def getTangent( self ):
    return self.Cep
    
  def __getYieldStress(self,alpha):
    return self.H * alpha + self.sigmaY
