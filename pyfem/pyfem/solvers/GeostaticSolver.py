# -*- coding: utf-8 -*-
from pyfem.util.BaseModule import BaseModule

from numpy import zeros
from pyfem.fem.Assembly import assembleMassMatrix, assembleTangentStiffness

class GeostaticSolver ( BaseModule ):
  
  def __init__  (self , props , globdat ):
    self.tol = 1.0e-6
    self.iterMax = 10
    
    BaseModule.__init__( self , props )
    
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
  def run( self , props , globdat ):
    
    globdat.cycle += 1
    
    dofCount = len(globdat.dofs)
    
    a = globdat.state
    Da = globdat.Dstate
    
    M,lumped = assembleMassMatrix( props , globdat )
    fext = globdat.fhat + lumped * globdat.gravity
    norm = globdat.dofs.norm(fext)
    
    Da[:] = zeros(dofCount)
    fint = zeros(dofCount)
    
    globdat.iiter = 0
    
    K,fint = assembleTangentStiffness(props,globdat)
    
    error = 1.
    
    while error > self.tol:
      globdat.iiter += 1
      da = globdat.dofs.solve(K,fext-fint)
      
      Da[:] += da[:]
      a[:] += da[:]
      
      K,fint = assembleTangentStiffness( props , globdat )
      
      if norm < 1.0e-16:
        error = globdat.dofs.norm(fext-fint)
      else:
        error = globdat.dofs.norm(fext-fint) / norm
        
      print ' Iter',globdat.iiter,' :',error
      
      if globdat.iiter == self.iterMax:
        raise RuntimeError('Newton-Raphson iterations did not converge!')
        
    globdat.elements.commitStress()
    a[:] -= Da[:]
    Da[:] = zeros(dofCount)
    fint = zeros(dofCount)
    
    globdat.iiter = 0
    
    K,fint = assembleTangentStiffness( props, globdat )
    error = 1.
    
    while error > self.tol:
      globdat.iiter += 1
      da = globdat.dofs.solve( K, fext-fint)
      
      Da[:] += da[:]
      a[:] += da[:]
      
      K,fint = assembleTangentStiffness( props, globdat )
      
      if norm < 1.0e-16:
        error = globdat.dofs.norm(fext-fint)
      else:
        error = globdat.dofs.norm(fext-fint) / norm
        
      print ' Iter',globdat.iiter,' :', error
      if globdat.iiter == self.iterMax:
        raise RuntimeError('Newton-Raphson iterations did not converge!')
      
    if da.max() > 1.0e-3:
      raise RuntimeError('GeoStatic Analysis did not converge!')
    
    globdat.elements.commitHistory()

    Da[:] = zeros(len(globdat.dofs))
    
    globdat.fint = fint
    
    globdat.active = False
    
    
