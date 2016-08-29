input = "OneElem.dat";

Elem =
{
  type = "SmallStrainContinuum3D";
  rho = 0;
  material = 
  {
    type = "Isotropic";
#    type = "J2Plastic";

    E = 2.E8;
    nu   = 0.3;
#    sigmaY = 2.E4;
#    H = 1.E8;
  };
};

steps=["Step-1"];

Step-1=
{
	type = "Static";
	solver = 
	{
#		type = "NonlinearSolver";
		type = "LinearSolver";
#		maxCycle = 1;
#		iterMax = 100;
#		tol = 1.0e-6;
	};
};

outlabel = ["stresses"];

outputModules = ["vtk","output"];

vtk =
{
  type = "Mesh3DWriter";
};

output =
{
  type = "OutputWriter";

  onScreen = true;
};