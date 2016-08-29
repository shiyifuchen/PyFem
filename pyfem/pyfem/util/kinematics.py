# -*- coding: utf-8 -*-
from numpy import zeros

class Kinematics:
  
  def __init__( self , nDim , nStr ):

    self.F      = zeros( shape=( nDim , nDim ) )
    self.E      = zeros( shape=( nDim , nDim ) )
    self.strain = zeros( nStr )
