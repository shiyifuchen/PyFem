# -*- coding: utf-8 -*-

from pyfem.materials.BaseMaterial import BaseMaterial
from numpy import zeros, dot,outer,inv
from numpy.linalg import eig
from math import sqrt
from pyfem.util.tensor import decomp,norm,unitFlow,getII,getI_bar

class TrescaPlaneStrain( BaseMaterial ):
  def __init__ ( self, props ):
    BaseMaterial.__init__( self, props )
    
    self.K = self.E/3.0/(1-2*self.nu)
    self.G = self.E/2.0/(1+self.nu)
    self.C = self.K*getII(2)+2*self.G*getI_bar(2)
    self.Cep = self.C.copy()
    self.setHistoryParameter("stresses",zeros(4))
    self.setHistoryParameter("equivalent_strain",0.0)
    self.commitHistory()
  
  def getStress( self, deformation ):
    sigma0=self.getHistoryParameter("stresses")
    eqStrain=self.getHistoryParameter("equivalent_strain") 
    
    dstrain = zeros(4)
    dstrain[[0,1,3]] = deformation.dstrain.copy()
    dstrain[3] *= 0.5
    sph_sigma0,dev_sigma0 = decomp(sigma0)
    sph_dstrain,dev_dstrain = decomp(dstrain)
    p_trial = sph_sigma0[0] + self.K*sum(sph_dstrain)
    deltaGamma = 0.0
    alpha = eqStrain
    S_trial = dev_sigma0 + 2 * self.G * dev_dstrain
    eigVals,eigVects = eig(S_trial)
    #eig计算的特征值未排序，特征向量是按特征值的顺序按列向量存储
    sorted = eigVals.argsort()  #返回一个排序的索引
    S1_trial = eigVals[sorted[0]]
    S2_trial = eigVals[sorted[1]]
    S3_trial = eigVals[sorted[2]]
    
    e1 = eigVects[:,sorted[0]]
    e2 = eigVects[:,sorted[1]]
    e3 = eigVects[:,sorted[2]]
    E1 = outer(e1,e1)
    E2 = outer(e2,e2)
    E3 = outer(e3,e3)
    
    g = S1_trial-S3_trial-4*self.G*deltaGamma-self.__getYieldStress(alpha)
    
    if g<=0:
      sigma = zeros(4)
      sigma[[0,1,3]] = sigma0[[0,1,3]] + dot(self.C, deformation.dstrain )
      sigma[2] = sigma0[2]+(self.K-2./3*self.G)*(deformation.dstrain[0]+\
              deformation.dstrain[1])
      self.setHistoryParameter("stresses",sigma)
      self.setHistoryParameter("equivalent_strain",eqStrain)
      return sigma[[0,1,3]], self.C      
    
    S , alpha = self.__oneVector([S1_trial,S2_trial,S3_trial] , alpha)
    
    if S[0]>S[1] and S[1]>S[2]:
      eqStrain = alpha
      sigma = (S[0]+p_trial)*E1 + (S[1]+p_trial)*E2 + (S[2]+p_trial)*E3
      
      ################### self.Cep #############################
      self.Cep = self.C                                        #
      ##########################################################
      self.setHistoryParameter("stresses",sigma)
      self.setHistoryParameter("equivalent_strain",eqStrain)
      return sigma, self.Cep
    
    # Return to corner
    
    if (S1_trial+S3_trial-2*S2_trial)>0:
      # Return to right corner
      S , alpha = self.__twoVectorRight([S1_trial,S2_trial,S3_trial] , alpha)
      eqStrain = alpha
      sigma = (S[0]+p_trial)*E1 + (S[1]+p_trial)*E2 + (S[2]+p_trial)*E3
      
      ################### self.Cep #############################
      self.Cep = self.C                                        #
      ##########################################################
      self.setHistoryParameter("stresses",sigma)
      self.setHistoryParameter("equivalent_strain",eqStrain)
      return sigma, self.Cep
    else:
      # Return to left corner
      S , alpha = self.__twoVectorLeft([S1_trial,S2_trial,S3_trial] , alpha)
      eqStrain = alpha
      sigma = (S[0]+p_trial)*E1 + (S[1]+p_trial)*E2 + (S[2]+p_trial)*E3
      
      ################### self.Cep #############################
      self.Cep = self.C                                        #
      ##########################################################
      self.setHistoryParameter("stresses",sigma)
      self.setHistoryParameter("equivalent_strain",eqStrain)
      return sigma, self.Cep    
    
    
    
  
  def getTangent( self ):
    pass
  
  def __getYieldStress( self , alpha ):
    return self.H * alpha + self.sigmaY
  
  # One-vector return mapping to main plane
  def __oneVector( self , S_trial , eqStrain ):
    deltaGamma = 0.0
    alpha = eqStrain
    g = S_trial[0] - S_trial[2] - self.__getYieldStress(alpha)
    
    k = 0
    while (abs(g) > 1e-8):
      k = k + 1
      if k > 10:
        raise RuntimeError('Newton-Raphson iterations did not converge!')
      Dg = -4 * self.G - self.H
      
      deltaGamma -= g/Dg
      alpha = eqStrain + deltaGamma
      g = S_trial[0]-S_trial[2]-4*self.G*deltaGamma-self.__getYieldStress(alpha)
      
    S = zeros(3)
    S[0] = S_trial[0] - 2 * self.G * deltaGamma
    S[1] = S_trial[1]
    S[2] = S_trial[2] + 2 * self.G * deltaGamma
    
    return S , alpha
    
  # Two-vector return mapping to right corner
  def __twoVectorRight( self , S_trial , eqStrain ):
    deltaGamma_a = 0.0
    deltaGamma_b = 0.0
    deltaGamma_bar = deltaGamma_a + deltaGamma_b
    deltaGamma = zeros(2)
    alpha = eqStrain
    
    g_a = S_trial[0]-S_trial[2]-2*self.G*(2*deltaGamma_a+deltaGamma_b)-self.__getYieldStress(alpha)
    g_b = S_trial[0]-S_trial[1]-2*self.G*(deltaGamma_a+2*deltaGamma_b)-self.__getYieldStress(alpha)
    g = zeros(2)
    g[0] = g_a
    g[1] = g_b
    Dg = zeros((2,2))
    Dg[0,0] = Dg[1,1] = -4*self.G - self.H
    Dg[0,1] = Dg[1,0] = -2*self.G - self.H
    
    k=0
    while (abs(g_a)+abs(g_b)) > 1e-8:
      k = k + 1
      if k > 10:
        raise RuntimeError('Newton-Raphson iterations did not converge!')
      deltaGamma = deltaGamma - dot(inv(Dg),g)
      deltaGamma_a = deltaGamma[0]
      deltaGamma_b = deltaGamma[1]
      deltaGamma_bar = deltaGamma_a+deltaGamma_b
      alpha = eqStrain + deltaGamma_bar
      g_a = S_trial[0]-S_trial[2]-2*self.G*(2*deltaGamma_a+deltaGamma_b)-self.__getYieldStress(alpha)
      g_b = S_trial[0]-S_trial[1]-2*self.G*(deltaGamma_a+2*deltaGamma_b)-self.__getYieldStress(alpha)
      g[0] = g_a
      g[1] = g_b
      
    S = zeros(3)
    S[0] = S_trial[0] - 2 * self.G * ( deltaGamma_a + deltaGamma_b )
    S[1] = S_trial[1] + 2 * self.G * deltaGamma_b
    S[2] = S_trial[2] + 2 * self.G * deltaGamma_a
    
    return S , alpha
    
    
    
  # Two-vector return mapping to left corner
  def __twoVectorLeft( self , S_trial , eqStrain ):
    deltaGamma_a = 0.0
    deltaGamma_b = 0.0
    deltaGamma_bar = deltaGamma_a + deltaGamma_b
    deltaGamma = zeros(2)
    alpha = eqStrain
    
    g_a = S_trial[0]-S_trial[2]-2*self.G*(2*deltaGamma_a+deltaGamma_b)-self.__getYieldStress(alpha)
    g_b = S_trial[1]-S_trial[2]-2*self.G*(deltaGamma_a+2*deltaGamma_b)-self.__getYieldStress(alpha)
    g = zeros(2)
    g[0] = g_a
    g[1] = g_b
    Dg = zeros((2,2))
    Dg[0,0] = Dg[1,1] = -4*self.G - self.H
    Dg[0,1] = Dg[1,0] = -2*self.G - self.H
    
    k=0
    while (abs(g_a)+abs(g_b)) > 1e-8:
      k = k + 1
      if k > 10:
        raise RuntimeError('Newton-Raphson iterations did not converge!')
      deltaGamma = deltaGamma - dot(inv(Dg),g)
      deltaGamma_a = deltaGamma[0]
      deltaGamma_b = deltaGamma[1]
      deltaGamma_bar = deltaGamma_a+deltaGamma_b
      alpha = eqStrain + deltaGamma_bar
      g_a = S_trial[0]-S_trial[2]-2*self.G*(2*deltaGamma_a+deltaGamma_b)-self.__getYieldStress(alpha)
      g_b = S_trial[0]-S_trial[1]-2*self.G*(deltaGamma_a+2*deltaGamma_b)-self.__getYieldStress(alpha)
      g[0] = g_a
      g[1] = g_b
      
    S = zeros(3)
    S[0] = S_trial[0] - 2 * self.G * deltaGamma_a
    S[1] = S_trial[1] - 2 * self.G * deltaGamma_b
    S[2] = S_trial[2] + 2 * self.G * ( deltaGamma_a + deltaGamma_b )
    
    return S , alpha