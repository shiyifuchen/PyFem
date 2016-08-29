# -*- coding: utf-8 -*-

from numpy import zeros, ones, ix_
from pyfem.util.dataStructures import Properties
from pyfem.util.dataStructures import elementData
from scipy.sparse import lil_matrix
import time


#######################################
# General array assembly routine for: # 
# * assembleInternalForce             #
# * assembleTangentStiffness          #
#######################################

def assembleArray ( props, globdat, rank, action ):
  t0=time.time()
  #Initialize the global array A with rank 2
  A = lil_matrix((len(globdat.dofs),len(globdat.dofs)))
  B = zeros( len(globdat.dofs) * ones(1,dtype=int) )

  globdat.resetNodalOutput()
  outlabel=[]
  if hasattr(props,'outlabel'):
    outlabel = getattr(props,'outlabel')

  #Loop over the element groups
  for elementGroup in globdat.elements.iterGroupNames():

    #Get the properties corresponding to the elementGroup
    el_props = getattr( props, elementGroup )
    
    #Loop over the elements in the elementGroup
    for element in globdat.elements.iterElementGroup( elementGroup ):

      #Get the element nodes
      el_nodes = element.getNodes()

      #Get the element coordinates
      el_coords = globdat.nodes.getNodeCoords( el_nodes )

      #Get the element degrees of freedom
      el_dofs = globdat.dofs.get( el_nodes )

      #Get the element state
      el_a  = globdat.state [el_dofs].copy()
      el_Da = globdat.Dstate[el_dofs].copy()
      
      factor1 = 1.0
      factor2 = 1.0
      
      if elementGroup in props.kill:
        el_a = zeros(el_a.shape)
        el_Da = zeros(el_Da.shape)
        factor1 = 0.0
        factor2 = 1e-6
        if hasattr(element,"mat"):
          element.mat.clearHistory()
#      if elementGroup == 'Elem1':
#        el_a = zeros(el_a.shape)
#        el_Da = zeros(el_Da.shape)

      #Create the an element state to pass through to the element
      #el_state = Properties( { 'state' : el_a, 'Dstate' : el_Da } )
      elemdat = elementData( el_a , el_Da )

      elemdat.coords   = el_coords
      elemdat.nodes    = el_nodes
      elemdat.props    = el_props
      elemdat.outlabel = outlabel
      
      if hasattr( element , "matProps" ):
        elemdat.matprops = element.matProps

      if hasattr( element , "mat" ):
        element.mat.reset()

      #Get the element contribution by calling the specified action
      getattr( element, action )( elemdat )
      
#      for label in elemdat.outlabel:	
#        element.appendNodalOutput( label , globdat , elemdat.outdata )
      

      if rank == 0:
        if elementGroup in props.kill:
          continue
        for i,label in enumerate(elemdat.outlabel):
          element.appendNodalOutput( label , globdat , elemdat.outdata[i] )
      elif rank == 1:
        B[el_dofs] += elemdat.fint*factor1
      elif rank == 2 and action is "getTangentStiffness":  
        A[ix_(el_dofs,el_dofs)] += elemdat.stiff*factor2
        B[el_dofs] += elemdat.fint*factor1
      elif rank == 2 and action is "getMassMatrix":  
        A[ix_(el_dofs,el_dofs)] += elemdat.mass*factor1
        B[el_dofs] += elemdat.lumped*factor1
      else:
        raise NotImplementedError('assemleArray is only implemented for vectors and matrices.')
#  A=A.tocsr()
  t1=time.time()
  print "Time Elapse for Assembly: ",t1-t0
  if rank == 1:
    return B
  elif rank == 2:
    return A.tocsr(),B


##########################################
# Internal force vector assembly routine # 
##########################################

def assembleInternalForce ( props, globdat ):
  return assembleArray( props, globdat, rank = 1, action = 'getInternalForce' )


#############################################
# Tangent stiffness matrix assembly routine # 
#############################################

def assembleTangentStiffness ( props, globdat ):
  return assembleArray( props, globdat, rank = 2, action = 'getTangentStiffness' )

#############################################
# Mass matrix assembly routine              # 
#############################################

def assembleMassMatrix ( props, globdat ):
  return assembleArray( props, globdat, rank = 2, action = 'getMassMatrix' )
  
def assembleOutputData ( props, globdat ):
  return assembleArray( props, globdat, rank = 0, action = 'getOutputData' )