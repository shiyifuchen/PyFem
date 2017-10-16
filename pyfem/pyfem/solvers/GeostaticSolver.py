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
    print "#######################   run"
    
    globdat.cycle += 1
    
    dofCount = len(globdat.dofs)
    
    a = globdat.state
    Da = globdat.Dstate
    
    M,lumped = assembleMassMatrix( props , globdat )
    fext = globdat.fhat + lumped * globdat.gravity
    norm = globdat.dofs.norm(fext,globdat.stepName)
    
    Da[:] = zeros(dofCount)
    fint = zeros(dofCount)
    
    globdat.iiter = 0
    
    K,fint = assembleTangentStiffness(props,globdat)
    K0=K.todense()
    
    error = 1.
    
    while error > self.tol:
      globdat.iiter += 1
      print "AAAA"
      da = globdat.dofs.solve(K,fext-fint,globdat.stepName)
      
      Da[:] += da[:]
      a[:] += da[:]
      
      K,fint = assembleTangentStiffness( props , globdat )
      print K0-K.todense()
      
      if norm < 1.0e-16:
        error = globdat.dofs.norm(fext-fint,globdat.stepName)
      else:
        error = globdat.dofs.norm(fext-fint,globdat.stepName) / norm
        
      print ' Iter',globdat.iiter,' :',error
      
      if globdat.iiter == self.iterMax:
        raise RuntimeError('Newton-Raphson iterations did not converge!')
#        
#    globdat.elements.commitStress()
#    a[:] -= Da[:]
#    Da[:] = zeros(dofCount)
##    fint = zeros(dofCount)
#    
#    globdat.iiter = 0
#    
##    K,fint = assembleTangentStiffness( props, globdat )
##    print fint
##    print fext
#
#    error = 1.
#
#    while error > self.tol:
#      globdat.iiter += 1
#      da = globdat.dofs.solve( K, fext-fint,globdat.stepName)
#      
#      Da[:] += da[:]
#      a[:] += da[:]
#      
#      K,fint = assembleTangentStiffness( props, globdat )
#      
#      if norm < 1.0e-16:
#        error = globdat.dofs.norm(fext-fint,globdat.stepName)
#      else:
#        error = globdat.dofs.norm(fext-fint,globdat.stepName) / norm
#        
#      print ' Iter',globdat.iiter,' :', error
#      if globdat.iiter == self.iterMax:
#        raise RuntimeError('Newton-Raphson iterations did not converge!')
#      
#    if da.max() > 1.0e-3:
#      raise RuntimeError('GeoStatic Analysis did not converge!')
#      
    globdat.elements.commitHistory()
    a[:] = 0

    Da[:] = zeros(len(globdat.dofs))
    
    globdat.fint = fint
    
    globdat.active = False
    
    
