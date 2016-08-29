input = "Test.dat";

Elem =
{
  type = "SmallStrainContinuum3D";

  material = 
  {
    type = "Isotropic";

    E = 2.1E11;
    nu   = 0.3;
  };
};

solver =
{
  #type = "LinearSolver";
  type = "NonlinearSolver";
  maxLam = 5.0;
};

#outputModules = ["vtk","output"];
outputModules = ["vtk"];

vtk =
{
  type = "Mesh3DWriter";
};

#output =
#{
#  type = "OutputWriter";
#
#  onScreen = true;
#};