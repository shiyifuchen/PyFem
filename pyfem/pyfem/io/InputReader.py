# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 20:35:12 2014

@author: LUX
"""

from pyfem.util.dataStructures import GlobalData

from pyfem.fem.NodeSet    import NodeSet
from pyfem.fem.ElementSet import ElementSet
from pyfem.fem.DofSpace   import DofSpace

from pyfem.util.fileParser import fileParser

import getopt,os.path

def InputReader( argv ):

#  options, remainder = getopt.getopt( argv , 'a:k:v', ['all','author='])

  proFileName  = argv[1]
  
  print "Input file name is : ",proFileName 

  props        = fileParser( proFileName )

  dataFileName = props.input
  
  print "Data file name is : ",dataFileName

  nodes = NodeSet()
  nodes.readFromFile( dataFileName )
  
  elems = ElementSet( nodes , props )
  elems.readFromFile( dataFileName )
  
  dofs = DofSpace( elems )
  dofs.readFromFile( dataFileName )

  globdat = GlobalData( nodes, elems, dofs ) 

  globdat.readFromFile( dataFileName )

  globdat.active = True
  globdat.cycle  = 0
  globdat.prefix = os.path.splitext(proFileName)[0]
	
  return props,globdat