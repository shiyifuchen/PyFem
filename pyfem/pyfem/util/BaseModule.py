# -*- coding: utf-8 -*-

class BaseModule:

  def __init__ ( self, props ):

    if hasattr(props,'currentModule') and hasattr(props,props.currentModule):
      currentModule = props.currentModule
    else:
      currentModule = self.__class__.__name__
      
      if currentModule.endswith("olver") is "olver":
        currentModule = "solver"
      
    if hasattr(props,currentModule):
      myProps = getattr(props,currentModule)

      for name,val in myProps:
        setattr( self, name, val )
