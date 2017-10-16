# -*- coding: utf-8 -*-

class Solver:

  def __init__( self , props , globdat ):

    solverProps = getattr( props, "solver" )

    solverType = solverProps.type

    exec "from pyfem.solvers."+solverType+" import "+solverType

    props.currentModule = "solver"

    self.solver = eval(solverType+"( props , globdat )")
    if globdat.fhatMap.has_key(globdat.stepName):
      globdat.fhat = globdat.fhatMap[globdat.stepName]
    

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
        
  def run( self , props , globdat ):

    self.solver.run( props , globdat )
