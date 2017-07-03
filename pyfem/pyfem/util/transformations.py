# -*- coding: utf-8 -*-
from numpy import array
from scipy.linalg import norm


def getRotationMatrix ( el_coords ):

  #Check the dimension of physical space
  if el_coords.shape[1] != 2:
    raise NotImplementedError('Rotation matrix only implemented for 2D situation')

  #Compute the (undeformed) element length
  l0 = norm( el_coords[1]-el_coords[0] )

  #Set up the rotation matrix to rotate a globdal
  #coordinate to an element coordinate (see Ch 1.3)
  sinalpha = (el_coords[1,1]-el_coords[0,1])/l0
  cosalpha = (el_coords[1,0]-el_coords[0,0])/l0

  return array([[cosalpha,sinalpha],[-sinalpha,cosalpha]])