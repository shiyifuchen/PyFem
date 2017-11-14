input = "OneElem.dat";
Elem =
{
  type = "SmallStrainContinuum3D";

  material = 
  {
    type = "J2Plastic"; 
    E = 2.0E8;
    nu   = 0.3;
    sigmaY = 1.0E5;
    H = 2.0E8;
  };
};

steps=["Step-1"];

Step-1=
{
	type = "Static";
	solver = 
	{
		type = "NonlinearSolver";
		maxCycle = 1;
		iterMax = 100;
		tol = 1.0e-6;
	};
};

outlabel = ["stresses","equivalent_strain"];

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