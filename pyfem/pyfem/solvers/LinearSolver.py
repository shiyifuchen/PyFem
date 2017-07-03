# -*- coding: utf-8 -*-

from pyfem.util.BaseModule import BaseModule

from numpy import zeros, array
from pyfem.fem.Assembly import assembleInternalForce, assembleTangentStiffness,assembleMassMatrix

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class LinearSolver ( BaseModule ):

  def __init__( self , props , globdat ):

    self.tol     = 1.0e-3
    self.iterMax = 10 

    BaseModule.__init__( self , props )

#    self.fext  = zeros( len(globdat.dofs) )  
 
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
   
  def run( self , props , globdat ):

    globdat.cycle += 1
    
    if globdat.hasGravity:
      M,lumped = assembleMassMatrix(props,globdat)
      fhat = globdat.fhat + lumped * globdat.gravity
    else:
      fhat = globdat.fhat
      
    K,fint = assembleTangentStiffness( props, globdat )
    
         
    globdat.state = globdat.dofs.solve( K, fhat )

    globdat.fint = assembleInternalForce( props, globdat )
    
    globdat.elements.commitHistory()

    globdat.active = False 
