from pyfem.util.BaseModule import BaseModule

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class Mesh3DWriter ( BaseModule ):
 
  def __init__( self , props , globdat ):
	
    self.prefix       = globdat.prefix
    self.elementGroup = "All"
    self.k            = 0
    self.interval     = 1

    BaseModule.__init__( self , props )

  def run( self , props , globdat ):
    
    if not globdat.cycle%self.interval == 0:
      return

    vtkfile = open( self.prefix + '-' + str(self.k) + '.vtu' ,'w' )

    vtkfile.write('<?xml version="1.0"?>\n')
    vtkfile.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n')
    vtkfile.write('<UnstructuredGrid>\n')
    vtkfile.write('<Piece NumberOfPoints="'+str(len(globdat.nodes))+'" NumberOfCells="')
    vtkfile.write(str(globdat.elements.elementGroupCount( self.elementGroup))+'">\n')
    vtkfile.write('<PointData>\n')
    vtkfile.write('<DataArray type="Float64" Name="displacement" NumberOfComponents="3" format="ascii" >\n')

    for nodeID in globdat.nodes.keys():
      for dofType in globdat.dofs.dofTypes:
        vtkfile.write(str(globdat.state[globdat.dofs.getForType(nodeID,dofType)])+' ')
   
      vtkfile.write('\n')
 
    vtkfile.write('</DataArray>\n')
  
#    stress = globdat.getData( "stresses" , range(len(globdat.nodes)) )

#    vtkfile.write('<DataArray type="Float64" Name="sigma_xx" NumberOfComponents="1" format="ascii" >\n')
#    for i in range(len(globdat.nodes)):
#      vtkfile.write( str(stress[i][0]) + " \n" )
#    vtkfile.write('</DataArray>\n')
#
#    vtkfile.write('<DataArray type="Float64" Name="sigma_yy" NumberOfComponents="1" format="ascii" >\n')
#    for i in range(len(globdat.nodes)):
#      vtkfile.write( str(stress[i][1]) + " \n" )
#    vtkfile.write('</DataArray>\n')
#
#    vtkfile.write('<DataArray type="Float64" Name="sigma_zz" NumberOfComponents="1" format="ascii" >\n')
#    for i in range(len(globdat.nodes)):
#      vtkfile.write( str(stress[i][2]) + " \n" )
#    vtkfile.write('</DataArray>\n')
#    
#    vtkfile.write('<DataArray type="Float64" Name="sigma_xy" NumberOfComponents="1" format="ascii" >\n')
#    for i in range(len(globdat.nodes)):
#      vtkfile.write( str(stress[i][3]) + " \n" )
#    vtkfile.write('</DataArray>\n')
#
#    vtkfile.write('<DataArray type="Float64" Name="sigma_yz" NumberOfComponents="1" format="ascii" >\n')
#    for i in range(len(globdat.nodes)):
#      vtkfile.write( str(stress[i][4]) + " \n" )
#    vtkfile.write('</DataArray>\n')
#
#    vtkfile.write('<DataArray type="Float64" Name="sigma_zx" NumberOfComponents="1" format="ascii" >\n')
#    for i in range(len(globdat.nodes)):
#      vtkfile.write( str(stress[i][5]) + " \n" )
#    vtkfile.write('</DataArray>\n')    

    print "Output label list : ",props.outlabel
    for label in props.outlabel:
      outdata = globdat.getData(label , range(len(globdat.nodes)) )
      column = len(outdata[0])
      vtkfile.write('<DataArray type="Float64" Name="'+label+'" NumberOfComponents="'+str(column)+'" format="ascii" >\n')
      for i in range(len(globdat.nodes)):
        for j in range(column):
          vtkfile.write( str(outdata[i][j]) + ' ' )
        vtkfile.write(" \n")
      vtkfile.write('</DataArray>\n')
	
    vtkfile.write('</PointData>\n')
    vtkfile.write('<CellData>\n')
    vtkfile.write('</CellData>\n')
    vtkfile.write('<Points>\n')
    vtkfile.write('<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">\n')
  
    for nodeID in globdat.nodes.keys():
      crd = globdat.nodes.getNodeCoords(nodeID)
      vtkfile.write( str(crd[0]) + ' ' + str(crd[1]) + ' ' + str(crd[2]) + " \n" )
    
    vtkfile.write('</DataArray>\n')
    vtkfile.write('</Points>\n')
    vtkfile.write('<Cells>\n')
    vtkfile.write('<DataArray type="Int64" Name="connectivity" format="ascii">\n')

    #--Store elements-----------------------------

    for element in globdat.elements.iterElementGroup( self.elementGroup ):
      el_nodes = globdat.nodes.getIndices(element.getNodes())

      if len(el_nodes) == 3 or len(el_nodes) == 4:
        for node in el_nodes:
          vtkfile.write(str(node)+' ')
  
      elif len(el_nodes) == 6 or len(el_nodes) == 8:
        for node in el_nodes:
          vtkfile.write(str(node)+' ')
  
      vtkfile.write('\n')
  
    vtkfile.write('</DataArray>\n')
    vtkfile.write('<DataArray type="Int64" Name="offsets" format="ascii">\n')

    for i,element in enumerate(globdat.elements.iterElementGroup( self.elementGroup )):
      el_nodes = globdat.nodes.getIndices(element.getNodes())
      vtkfile.write(str(len(el_nodes)*(i+1))+'\n')

    vtkfile.write('</DataArray>\n')
    vtkfile.write('<DataArray type="UInt8" Name="types" format="ascii" RangeMin="12" RangeMax="12">\n')

    for i in range(globdat.elements.elementGroupCount( self.elementGroup)):
      vtkfile.write('12\n')

    vtkfile.write('</DataArray>\n')
    vtkfile.write('</Cells>\n')
    vtkfile.write('</Piece>\n')
    vtkfile.write('</UnstructuredGrid>\n')
    vtkfile.write('</VTKFile>\n')
 
    #--Write pvd file

    f = open( self.prefix + '.pvd' ,'w' )

    f.write("<VTKFile byte_order='LittleEndian' type='Collection' version='0.1'>\n")
    f.write("<Collection>\n")
  
    for i in range(self.k+1):
      f.write("<DataSet file='"+self.prefix+'-'+str(i)+".vtu' groups='' part='0' timestep='"+str(i)+"'/>\n")
   
    f.write("</Collection>\n")
    f.write("</VTKFile>\n")

    f.close()
  
    self.k = self.k+1
