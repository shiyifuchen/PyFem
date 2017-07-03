# -*- coding: utf-8 -*-
from .Element import Element
from pyfem.util.transformations import getRotationMatrix

from numpy import zeros, dot, eye, outer, empty
from scipy.linalg import norm

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class KirchhoffBeam ( Element ):

  #dofs per element
  dofTypes = [ 'u' , 'v' , 'rz' ]

  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )

    self.EA = self.E * self.A
    self.EI = self.E * self.I

    self.intpoints = zeros(2)
    self.weights   = zeros(2)
    
    self.weights[0] = 1.0
    self.weights[1] = 1.0
        
    self.intpoints[0] = -0.57735
    self.intpoints[1] =  0.57735

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def __type__ ( self ):
    return self.name

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):

    l0  = norm( elemdat.coords[1]-elemdat.coords[0] )
    jac = 0.5 * l0
    
    a_bar = self.glob2Elem( elemdat.state , elemdat.coords )
            
    fint  = zeros(6);
    stiff = zeros( elemdat.stiff.shape ); 
           
    for xi,alpha in zip( self.intpoints , self.weights ):
      
      bu = self.getBu( l0 , xi )
      bw = self.getBw( l0 , xi )
      c  = self.getC ( l0 , xi )
      
      epsl = dot( bu , a_bar ) + 0.5*(dot( bw , a_bar ) )**2
      chi  = dot( c  , a_bar )
             
      N = self.EA * epsl
      M = self.EI * chi
              
      wght = jac * alpha
    
      fint  += N * bu * wght
      fint  += ( N * dot( bw , a_bar ) * bw + M * c ) * wght
    
      stiff += self.EA * outer( bu , bu ) * wght
      stiff += self.EA * dot( bw , a_bar ) * outer( bu , bw ) * wght
      stiff += self.EA * dot( bw , a_bar ) * outer( bw , bu ) * wght
      stiff += ( self.EI * outer( c , c ) + \
                 self.EA * (dot( bw , a_bar ))**2 * outer( bw , bw ) + \
                 N  * outer( bw , bw ) ) * wght
    
    elemdat.fint  = self.elem2Glob( fint  , elemdat.coords )
    elemdat.stiff = self.elem2Glob( stiff , elemdat.coords )    

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
    
  def getInternalForce ( self, elemdat ):

    EA = elemdat.props.E * elemdat.props.A
    EI = elemdat.props.E * elemdat.props.I
    l0 = norm( elemdat.coords[1]-elemdat.coords[0] )

    intpoints = zeros(2)
    weights   = zeros(2)
    
    weights[0] = 1.0
    weights[1] = 1.0
        
    intpoints[0] = -0.57735
    intpoints[1] =  0.57735

    a_bar = elemdat.state
            
    for xi,alpha in zip(intpoints,weights):
      
      bu = self.getBu( l0 , xi )
      bw = self.getBw( l0 , xi )
      c  = self.getC ( l0 , xi )
      
      epsl = dot( bu , a_bar ) + 0.5*(dot( bw , a_bar ) )**2
      chi  = dot( c  , a_bar )
             
      N    = EA * epsl
      M    = EI * chi
              
      wght = 0.5 * l0 * alpha
    
      elemdat.fint  += N * bu * wght
      elemdat.fint  += ( N * dot( bw , a_bar ) * bw + M * c ) * wght
  
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getBu( self , l0 , xi ):

    Bu = zeros(6 )

    Bu[0] = -1.0/l0
    Bu[3] =  1.0/l0
  
    return Bu

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getBw( self , l0 , xi ):

    Bw = zeros(6)

    Bw[1] = 1.5/l0*(xi*xi-1.0)
    Bw[2] = 0.25*(3*xi*xi-2.0*xi-1.0)
    Bw[4] = -1.5/l0*(xi*xi-1.0)
    Bw[5] = 0.25*(3*xi*xi+2.0*xi-1.0)
        
    return Bw

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getC( self , l0 , xi ):

    C = zeros( 6 )

    l02 = l0*l0
    
    C[1] = 6.0*xi/l02
    C[2] = (3.0*xi-1.0)/l0
    C[4] = -6.0*xi/l02
    C[5] = (3.0*xi+1.0)/l0
       
    return C

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def glob2Elem( self , a , coords ):

    a_bar = empty( a.shape )

    R     = eye( 6)
    crd   = zeros( shape=(2,2) )
        
    crd[0,:] = coords[0,:]
    crd[1,:] = coords[1,:]

    R[:2,:2]   = getRotationMatrix( crd )
    R[3:5,3:5] = R[:2,:2]

    a_bar      = dot( R , a )

    return a_bar

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def elem2Glob( self , a_bar , coords ):

    a = empty( a_bar.shape )
    R     = eye( 6)
    crd   = zeros( shape=(2,2) )
    
    crd[0,:] = coords[0,:]
    crd[1,:] = coords[1,:]
    
    R[:2,:2]   = getRotationMatrix( crd )
    R[3:5,3:5] = R[:2,:2]
    
    if len(a_bar.shape) == 1:
      a        = dot( R.transpose() , a_bar )
    elif len(a_bar.shape) == 2:
      a        = dot( R.transpose(), dot ( a_bar , R ) )
    return a
    
  def getOutputData( self , elemdat ):
    pass