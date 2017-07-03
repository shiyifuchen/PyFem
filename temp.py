# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 13:48:48 2017

@author: LUX
"""

from numpy import zeros,array
from numpy import linalg

def calPrincipal(stress):
  A=zeros([3,3])
  A[0,0]=stress[0]
  A[1,1]=stress[1]
  A[2,2]=stress[2]
  A[0,1]=A[1,0]=stress[3]
  A[1,2]=A[2,1]=stress[4]
  A[2,0]=A[0,2]=stress[5]
  return linalg.eig(A)[0]
  
#  I1=stress[0]+stress[1]+stress[2]
#  I2=-stress[0]*stress[1]-stress[1]*stress[2]-stress[2]*stress[0]\
#    +stress[3]**2+stress[4]**2+stress[5]**2
#  I3=stress[0]*stress[1]*stress[2]+2*stress[3]*stress[4]*stress[5]\
#    -stress[0]*(stress[4]**2)-stress[1]*(stress[5]**2)-stress[2]*(stress[3]**2)
if __name__=="__main__":
  sigma=array([10,-10,10,0,0,-10])
  print calPrincipal(sigma)
  