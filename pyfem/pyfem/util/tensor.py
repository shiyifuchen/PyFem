# -*- coding: utf-8 -*-

from numpy import zeros,ones,eye
from math import sqrt

def decomp(T):
  if len(T)==6 or len(T)==4:
    spherical_part = zeros(len(T))
    deviatoric_part = T.copy()
    p = (T[0]+T[1]+T[2])/3.0
    for i in range(3):
      spherical_part[i]=p
      deviatoric_part[i]-=p
    return spherical_part,deviatoric_part
  elif len(T)==3:
    spherical_part = zeros(3)
    deviatoric_part = T.copy()
    p = (T[0]+T[1])/2.0
    for i in range(2):
      spherical_part[i]=p
      deviatoric_part[i]-=p
    return spherical_part,deviatoric_part
  else:
    raise RuntimeError("length of the vector should be 3 or 6 !")

  
def norm(T):
  if len(T)==6 or len(T)==4:
    dim = 3
  elif len(T)==3:
    dim = 2
  else:
    raise RuntimeError("length of the vector should be 3 or 6 !")
  sum = 0.0
  for i in range(len(T)):
    sum += T[i] * T[i]
    if i>=dim:
      sum += T[i] * T[i]
  return sqrt(sum)
  
def unitFlow(T):
  return T/norm(T)
  
def invariant1(T):
  if len(T)==6 or len(T)==4:
    dim = 3
  elif len(T)==3:
    dim = 2
  else:
    raise RuntimeError("length of the vector should be 3 or 6 !")
  return T[:dim].sum()
  
def invariant2(T):
  return 0.5*norm(T)  
  
def getI(dim):
  length = 0
  if dim==2:
    length = 3
  elif dim==3:
    length = 6
  else:
    raise RuntimeError("dimension of the tensor should be 2 or 3 !")
  T=eye(length)
  T[dim:,dim:] *= 0.5
  return T
  
def getII(dim):
  length = 0
  if dim==2:
    length = 3
  elif dim==3:
    length = 6
  else:
    raise RuntimeError("dimension of the tensor should be 2 or 3 !")
  T=zeros((length,length))
  T[:dim,:dim]=ones(dim)
  return T

def getI_bar(dim):
  return getI(dim)-1./3*getII(dim)