# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 20:47:10 2014

@author: LUX
"""

#可以扩展成一个完善的数据库管理系统

class itemList ( dict ):

  def add ( self, ID, item ):
    
    if ID in self:
      raise RuntimeError( 'ID ' + str(ID) + ' already exists in ' + type(self).__name__ )

    self[ID] = item

  def get ( self, IDs ):

    if isinstance(IDs,int):
      return self[IDs]
    elif isinstance(IDs,list):
      return [ self[ID] for ID in IDs ]
      
    raise RuntimeError('illegal argument for itemList.get')  

  def getIndices ( self, IDs ):
    
    if isinstance(IDs,int):
      return self.keys().index( IDs )
    elif isinstance(IDs,list):
      return [ self.keys().index( ID ) for ID in IDs ]
      
    raise RuntimeError('illegal argument for itemList.getIndices')  