from pyfem.materials.BaseMaterial import BaseMaterial
from numpy import zeros,dot,linalg
from math import log10,pi,sin,cos
from pyfem.util.tensor import getII,getI_bar

class DuncanChangEv(BaseMaterial):
  def __init__(self,props):
    BaseMaterial.__init__(self,props)
    
#    self.K = self.E/3.0/(1-2*self.nu)
#    self.G = self.E/2.0/(1+self.nu)
#    self.C = self.K*getII(3)+2*self.G*getI_bar(3)    
    
    self.setHistoryParameter("stresses",zeros(6))
    self.setHistoryParameter("confining_pressure",0.0)
    self.setHistoryParameter("max_stress_diff",0.0)
    self.setHistoryParameter("max_stress_level",0.0)
    self.commitHistory()
    
    #PROPS:KK,n,Rf,c,Phi0,G,D,F,Kur,Pa,DPhi
  def getStress( self, deformation ):
    sigma0=self.getHistoryParameter("stresses")
    CP=self.getHistoryParameter("confining_pressure")
    MSD=self.getHistoryParameter("max_stress_diff")
    MSL=self.getHistoryParameter("max_stress_level")
    PS = self.__calPrincipal(sigma0)
    PS1=-PS[2]
    PS3=-PS[0]
    Phi = self.Phi0-self.DPhi*log10(CP/self.Pa)
    Phi=Phi*pi/180
    S=self.__updateEnu(PS1,PS3,Phi)
    self.K = self.E/3.0/(1-2*self.nu)
    self.G = self.E/2.0/(1+self.nu)
    self.C = self.K*getII(3)+2*self.G*getI_bar(3)    
    dstrain = deformation.dstrain.copy()
    dstrain[3:] *= 0.5
    sigma = sigma0 + dot(self.C, deformation.dstrain )
    PS = self.__calPrincipal(sigma)
    PS1=-PS[2]
    PS3=-PS[0]    
    S=self.__updateEnu(PS1,PS3,Phi)
    self.K = self.E/3.0/(1-2*self.nu)
    self.G = self.E/2.0/(1+self.nu)
    self.C = self.K*getII(3)+2*self.G*getI_bar(3)  
    if PS3>CP:
      self.setHistoryParameter("confining_pressure",PS3)
    if (PS3-PS1)>MSD:
      self.setHistoryParameter("max_stress_diff",PS3-PS1)
    if S>MSL:
      self.setHistoryParameter("max_stress_level",S)
    self.setHistoryParameter("stresses",sigma)
    return sigma,self.C
  
  def getTangent( self ):
    return self.C

  def __calPrincipal(self,stress):
    A=zeros([3,3])
    A[0,0]=stress[0]
    A[1,1]=stress[1]
    A[2,2]=stress[2]
    A[0,1]=A[1,0]=stress[3]
    A[1,2]=A[2,1]=stress[4]
    A[2,0]=A[0,2]=stress[5]
    return linalg.eig(A)[0]
  def __updateEnu(self,PS1,PS3,Phi):
    CP=self.getHistoryParameter("confining_pressure")
    MSD=self.getHistoryParameter("max_stress_diff")
    MSL=self.getHistoryParameter("max_stress_level")
    S=(1-sin(Phi))*(PS1-PS3)
    PSFEI = PS3
    if PS3<0.1:
      PSFEI=0.1
    S=S/(2*self.c*cos(Phi)+2*PSFEI*sin(Phi))
    if S >= 0.99:
      S=0.99
    A=self.D*(PS1-PS3)
    A=A/(self.KK*self.Pa*((CP/self.Pa)**self.n))
    A=A/(1-self.Rf*S)
    self.nu=self.G-self.F*log10(CP/self.Pa)
    self.nu=self.nu/(1-A)/(1-A)
    if self.nu>=0.49:
      self.nu=0.49
    self.E=self.KK*self.Pa*((CP/self.Pa)**self.n)
    self.E=self.E*((1-self.Rf*S)**2)
    if (S<MSL) and ((PS1-PS3)<MSD):
      self.E=self.Kur*self.Pa*((CP/self.Pa)**self.n)
    return S
    
    