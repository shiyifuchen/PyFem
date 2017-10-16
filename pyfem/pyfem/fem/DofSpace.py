# -*- coding: utf-8 -*-

from numpy import array, dot, zeros
import scipy.linalg
from pyfem.util.itemList import itemList
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import time

class DofSpace:

  def __init__ ( self, elements ):

    self.dofTypes = elements.getDofTypes()
    self.dofs     = array( range( len(elements.nodes) * len(self.dofTypes) ) ).reshape( ( len(elements.nodes), len(self.dofTypes) ) )
    self.nodes    = elements.nodes

    #Create the ID map
    self.IDmap = itemList()
    for ind,ID in enumerate(elements.nodes):
      self.IDmap.add( ID, ind )

    self.constrainedDofs = []
    self.constrainedVals = []
    self.constrainedFac  = 1.
    #########2017.07.08
    # For Multi-Step Analysis
    self.constrainedDofsMap={}
    self.constrainedValsMap={}
    ###################

  def __str__ ( self ):
    return str(self.dofs)

  def __len__ ( self ):
    return len(self.dofs.flatten())

  def setConstrainFactor( self , fac ):
    self.constrainedFac = fac
    
  def readFromFile( self, fname ):
    
    fin = open( fname )

    while True:
      line = fin.readline()  
  
      if line.startswith('<NodeConstraints>') == True:
        while True:
          line = fin.readline()  

          if line.startswith('</NodeConstraints>') == True:
            ##### 如何存在默认的边界条件，在解析完所有边界条件后先设置为默认的
            if self.constrainedDofsMap.has_key("Default"):
              self.constrainedDofs=self.constrainedDofsMap["Default"]
              self.constrainedVals=self.constrainedValsMap["Default"]
            return
            
          ######## 2017.07.08
          if line.startswith('<step') ==True:
            ####有初始的边界条件，需要保存
            if len(self.constrainedDofs)>0 and len(self.constrainedDofsMap)== 0:
              self.constrainedDofsMap["Default"]=self.constrainedDofs
              self.constrainedValsMap["Default"]=self.constrainedVals
            a=(line[1:-1]).strip().split(',')
            stepName=''
            opType='MOD'
            if len(a)==1:
              stepName=(a[0].split('=')[1]).strip()
            elif len(a)==2:
              stepName=(a[0].split('=')[1]).strip()
              opType=(a[1].split('=')[1]).strip()
              opType=(opType.split('>')[0]).strip()              
              
              if opType!="MOD":
                opType="NEW"
                #####如果为NEW模式，则清除之前的边界条件
                self.constrainedDofs = []
                self.constrainedVals = []
            while True:
              line = fin.readline()
              if line.startswith('</step>') == True:
                #####分析步字段结束时以分析步名称保存当前的边界条件
                self.constrainedDofsMap[stepName]=self.constrainedDofs
                self.constrainedValsMap[stepName]=self.constrainedVals
                break
              a = line.strip().split(';')
              if len(a) == 2:
                b = a[0].split('=')
                if len(b) == 2:
                  c = b[0].split('[')
                  dofType = c[0]
                  nodeID  = eval(c[1].split(']')[0])
                  self.constrain( nodeID , dofType , eval(b[1]))
                
        
          a = line.strip().split(';')
      
          if len(a) == 2:
            b = a[0].split('=')
        
            if len(b) == 2:
              c = b[0].split('[')
              
              dofType = c[0]
              nodeID  = eval(c[1].split(']')[0])
              
              self.constrain( nodeID , dofType , eval(b[1]))
              
  def constrain ( self, nodeID, dofTypes , val = 0. ):

    if not nodeID in self.nodes:
      raise RuntimeError('Node ID ' + str(nodeID) + ' does not exist')

    ind = self.IDmap.get( nodeID )

    if isinstance( dofTypes, str ):
      dofTypes = [dofTypes]

    #Check if the dofTypes exist
    for dofType in dofTypes:
      if dofType not in self.dofTypes:
        raise RuntimeError('DOF type "' + dofType + '" does not exist')
      
    for dofType in dofTypes:
      ######### 2017.07.09
      ##为了提供对多分析步的支持，允许对指定自由度的约束值进行修改
      dofIndex = self.dofs[ind,self.dofTypes.index(dofType)]
      if dofIndex in self.constrainedDofs:
        self.constrainedVals[self.constrainedDofs.index(dofIndex)] = val
      else:
        self.constrainedDofs.append( dofIndex )
        self.constrainedVals.append( val )
      ############################################################

  #获得一组节点的某个自由度的索引列表
  def getForType ( self, nodeIDs, dofType ):
    return self.dofs[self.IDmap.get( nodeIDs ),self.dofTypes.index(dofType)]
    
  #获得一组节点所有自由度的索引列表
  def get ( self, nodeIDs ):
    return self.dofs[self.IDmap.get(nodeIDs)].flatten()

  def getConstraintsMatrix ( self ,name):
    ##### 如果没有保存的名字则沿用上一个分析步的约束
    if self.constrainedDofsMap.has_key(name):
      self.constrainedDofs=self.constrainedDofsMap[name]
      self.constrainedVals=self.constrainedValsMap[name]
    n_constrained = len( self.constrainedDofs )
    n             = len( self )

    C = lil_matrix( (n,n-n_constrained),dtype=int )
    j = 0
    
    for i in range(n):
      
      if i in self.constrainedDofs:
        continue

      C[i,j] = 1.
      j+=1
      
    Ct=C.transpose()
    
    return C.tocsr(),Ct.tocsr()

  def solve ( self, A, b ,name):
    t0=time.time()
    if len(A.shape) == 2:
      C,Ct = self.getConstraintsMatrix(name)

      a = zeros(len(self))
      a[self.constrainedDofs] = self.constrainedFac * array(self.constrainedVals)
      
      A_constrained = Ct*A*C
      b_constrained = Ct*( b + A* (-a ) )
      x_constrained = spsolve( A_constrained, b_constrained )
      x =  C.dot( x_constrained )

      x[self.constrainedDofs] = self.constrainedFac * array(self.constrainedVals)
    
    elif len(A.shape) == 1:
      x = b / A

      x[self.constrainedDofs] = self.constrainedFac * array(self.constrainedVals)
    t1=time.time()
    print "Time Elapse for Solve: ",t1-t0
    return x
    
  def eigensolve( self, A , B ,name):
    #To be implemented!
    return A,B

  def norm ( self, r,name ):

    C,Ct = self.getConstraintsMatrix(name)
    return scipy.linalg.norm( Ct*r )
