# -*- coding: utf-8 -*-
class BaseMaterial:

  def __init__ ( self, props ):
    for name,val in props:
      setattr( self, name, val )

    self.initHistory={}
    self.current = []
    self.iIter = -1

  def setIter( self, iIter ):
    
    self.iIter = iIter

  def setHistoryParameter( self , name , val ):

    if self.iIter == -1:
      self.initHistory[name]=val
      return

    if len(self.current) == self.iIter:
      self.current.append(self.initHistory.copy())
 
    self.current[self.iIter][name] = val
    
    return 
   
  def getHistoryParameter( self , name ):

    if len(self.history) == 0:
      return self.initHistory[name]
    else:
      return self.history[self.iIter][name]
    
  def commitHistory( self ):
#    print "Material CommitHistory ......"
    self.history = []

    #对积分点循环
    for h in self.current:
      #此处需要使用copy()！
      self.history.append(h.copy())
      
  def commitStress( self ):
    
    self.history = []

    for i,h in enumerate(self.current):
      self.history.append(self.initHistory.copy())
      self.history[i]["stresses"] = self.current[i]["stresses"]
      
      
  def clearHistory( self ):
    if hasattr(self,"history"):
      self.history = []
