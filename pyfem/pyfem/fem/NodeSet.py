# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 21:03:59 2014

@author: LUX
"""

from numpy import array
from pyfem.util.itemList import itemList
import re

class NodeSet( itemList ):

  def getNodeCoords( self, nodeIDs ):
    return array( self.get( nodeIDs ) )
    
  def readFromFile( self, fname ):
    
    fin = open( fname )
  
    while True:
    
      line = fin.readline()  
  
      if line.startswith('<Nodes>') == True:
      
        while True:
          line = fin.readline()  

          if line.startswith('</Nodes>') == True:
            return

          #对于空格数大于等于2的情况，将被一个空格代替  
          line = re.sub('\s{2,}',' ',line)
          a = line.split(';')
     
          for a in a[:-1]:
            b = a.strip().split(' ')
            
            if b[0].startswith("//") or b[0].startswith("#"):
              break
            if len(b) > 1 and type(eval(b[0])) == int:
              self.add( eval(b[0]), [eval(crd) for crd in b[1:]] )          
