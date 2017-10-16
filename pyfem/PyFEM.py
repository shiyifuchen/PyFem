# -*- coding: utf-8 -*-

from pyfem.io.InputReader   import InputReader
from pyfem.io.OutputManager import OutputManager
from pyfem.solvers.Solver   import Solver

import sys
import time

t0=time.time()
print "******  InputReader Start ..."

props,globdat = InputReader( sys.argv )
print props.Elem.material.E
props.kill = []

t1=time.time()
print "******  InputReader Time Elapse: ",t1-t0

output = OutputManager(props,globdat)

for step in props.steps:
  print "Step '",step,"' Start ..."
  step_props = getattr(props,step)
  globdat.stepName = step
  
  if hasattr(step_props,"kill"):
    props.kill = list(set(props.kill + step_props.kill))
  if hasattr(step_props,"active"):
    for name in step_props.active:
      if name in props.kill:
        props.kill.remove(name)

  solver = Solver(step_props,globdat)
  
  while globdat.active:
    solver.run(props,globdat)
    print globdat.getDisp(100,'w')
    output.run(props,globdat)
  globdat.active = True
  globdat.cycle = 0
  print "Step '",step,"' Finished ..."
  
print "PyFem analysis terminated successfully."