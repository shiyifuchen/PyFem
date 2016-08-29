from pyfem.fem.Assembly import assembleOutputData
import time

class OutputManager:

  def __init__( self , props , globdat ):

    self.outman = []

    outputModules = props.outputModules

    for name in outputModules:
   
      props.currentModule = name

      ioType = name

      if hasattr( props , name):
       moduleProps = getattr( props, name )
       if hasattr( moduleProps , "type" ):
         ioType = moduleProps.type

      exec "from pyfem.io."+ioType+" import "+ioType

      self.outman.append(eval(ioType+"( props , globdat )"))

  def run( self , props , globdat ):
    assembleOutputData(props,globdat)

    for i,output in enumerate(self.outman):
      print "******  ",props.outputModules[i]," Start ..."
      t0 = time.time()
      output.run( props , globdat )
      t1 = time.time()
      print "Time Elapse for ",props.outputModules[i],": ",t1-t0