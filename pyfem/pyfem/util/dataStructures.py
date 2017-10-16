# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 17:30:21 2014

@author: LUX
"""

from numpy import zeros,linalg

class Properties:
  def __init__ ( self, dictionary = {} ):

    for key in dictionary.iterkeys():
      setattr( self, key, dictionary[key] )

  def __str__ ( self ):

    myStr  = ''
    for att in dir( self ):
      
      #Ignore private members and standard routines
      if att.startswith('__'):
        continue
      
      myStr += 'Attribute: ' + att + '\n'
      myStr += str( getattr(self,att) ) + '\n'

    return myStr

  def __iter__ ( self ):

    propsList = []
    for att in dir( self ):
      
      #Ignore private members and standard routines
      if att.startswith('__'):
        continue
      
      propsList.append( ( att, getattr(self,att) ) )

    return iter(propsList)

  def store ( self , key , val ):
    setattr( self , key , val )

class GlobalData ( Properties ):
  
  def __init__ ( self, nodes, elements, dofs ):

    Properties.__init__( self, { 'nodes' : nodes, 'elements' : elements, 'dofs' : dofs } )

    self.state  = zeros( len( self.dofs ) )
    self.Dstate = zeros( len( self.dofs ) )
    self.fint   = zeros( len( self.dofs ) )
    self.fhat   = zeros( len( self.dofs ) )
    self.gravity = zeros( len( self.dofs ) )

    self.velo   = zeros( len( self.dofs ) )
    self.acce   = zeros( len( self.dofs ) )

    self.cycle  = 0
    self.iiter  = 0
    self.time   = 0.0
   
    self.outputNames = []
    self.hasGravity = False
    self.stepName = ''
    self.fhatMap = {}

  def readFromFile( self , fname ):

    fin = open( fname )
  
    while True:
      line = fin.readline()  
  
      if line.startswith('<ExternalForces>') == True:
        while True:
          line = fin.readline()  

          if line.startswith('</ExternalForces>') == True:
          ##### 如果存在默认的载荷，在解析完所有载荷条件后先设置为默认的
            if self.fhatMap.has_key("Default"):
              self.fhat = self.fhatMap["Default"]
            return
          
          ########## 2017.07.09
          if line.startswith('<Gravity>') == True:
            while True:
              line = fin.readline()
              
              if line.startswith('</Gravity>') == True:
                break
              a = line.strip().split(';')
              if len(a) == 2:
                b = a[0].split('=')
                dofType = b[0].strip()
                self.gravity[self.dofs.getForType(self.nodes.keys(),dofType)] = eval(b[1])
                self.hasGravity = True                
          
          if line.startswith('<step') ==True:
            #####又初始的载荷条件，需要保存
            if linalg.norm(self.fhat) > 0 and len(self.fhatMap) == 0:
              self.fhatMap['Default'] = self.fhat
            a = a=(line[1:-1]).strip().split(',')
            stepName=''
            opType='MOD'
            if len(a) == 1:
              stepName=(a[0].split('=')[1]).strip()
            elif len(a) == 2:
              stepName=(a[0].split('=')[1]).strip()
              opType=(a[1].split('=')[1]).strip()
              opType=(opType.split('>')[0]).strip()
              if opType!='MOD':
                opType='NEW'
                #####如果为NEW模式，则清除之前的载荷条件
                self.fhat = zeros( len( self.dofs ) )
            while True:
              line = fin.readline()
              if line.startswith('</step>') == True:
                #####分析步字段结束时以分析步名称保存当前的边界条件
                self.fhatMap[stepName]=self.fhat
                break
              a = line.strip().split(';')
              if len(a) == 2:
                b = a[0].split('=')
                if len(b) == 2:
                  c = b[0].split('[')
                  dofType = c[0]
                  nodeID  = eval(c[1].split(']')[0])
                  self.fhat[self.dofs.getForType(nodeID,dofType)] = eval(b[1])
        
          a = line.strip().split(';')
          
          if len(a) == 2:
            b = a[0].split('=')
        
            if len(b) == 2:
              c = b[0].split('[')
              
              dofType = c[0]
              nodeID  = eval(c[1].split(']')[0])
              
              self.fhat[self.dofs.getForType(nodeID,dofType)] = eval(b[1])
              
#      if line.startswith('<Gravity>') == True:
#        while True:
#          line = fin.readline()  
#
#          if line.startswith('</Gravity>') == True:
#            break
#        
#          a = line.strip().split(';')
#          
#          if len(a) == 2:
#            b = a[0].split('=')
#            dofType = b[0].strip()
#            self.gravity[self.dofs.getForType(self.nodes.keys(),dofType)] = eval(b[1])
#            self.hasGravity = True

#---------------------------------------------------------------------------------
#
#---------------------------------------------------------------------------------

  def printNodes( self , inodes=None ):
    

    if inodes is None:
      inodes = self.nodes.keys() 
	
    print '   Node | ',
    
    for dofType in self.dofs.dofTypes:
      print "  '%s'        " % dofType,

    if hasattr( self , 'fint' ):
      for dofType in self.dofs.dofTypes:
        print " fint'%s'   " % dofType,

    for name in self.outputNames:
      print "       '%s'             " % name,

    print 
    print '-------------------------------------------------------',
    print '-------------------------------------------------------'

    for nodeID in inodes:
      print '  %4i  | ' % (nodeID),
      for dofType in self.dofs.dofTypes:
        print ' %10.3e ' % self.state[self.dofs.getForType(nodeID,dofType)],
      for dofType in self.dofs.dofTypes:
        print ' %10.3e ' % self.fint[self.dofs.getForType(nodeID,dofType)],

      for name in self.outputNames:
        data = self.getData( name , nodeID )
      
        for a in data:
          print ' %10.3e ' % a,
      print
    print

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
  def getDisp( self , nodeID , dofType ):
    return self.state[self.dofs.getForType(nodeID,dofType)]
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
      
  def getData( self , outputName , inodes ):

    data    = getattr( self, outputName )
    weights = getattr( self, outputName + 'Weights' )

    if type(inodes) is int:
      i = self.nodes.keys().index(inodes)
      if weights[i] != 0:
        return data[i,:] / weights[i]
      else:
        return data[i,:]
    else:
      outdata=[]

      for row,w in zip(data[inodes,:],weights[inodes]):
        if w != 0:
          outdata.append(row/w)
        else:
          outdata.append(row)

      return outdata

  def resetNodalOutput ( self ):

    for outputName in self.outputNames:
      delattr( self , outputName )
      delattr( self , outputName + 'Weights' )

    self.outputNames=[]

class elementData():

  def __init__( self , elstate , elDstate ):

    nDof        = len(elstate)

    self.state  = elstate
    self.Dstate = elDstate
    self.stiff  = zeros( shape=( nDof,nDof ) )
    self.fint   = zeros( shape=(nDof) )
    self.mass   = zeros( shape=( nDof,nDof ) )
    self.lumped = zeros( shape=(nDof) )

    self.outlabel = []
   
  def __str__( self ):

    return self.state