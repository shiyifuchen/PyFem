input = "2DSlope.dat";

Elem =
{
	type = "SmallStrainContinuum";
	
	rho = 2000;
	
	material = 
	{
		#type = "PlaneStrain";
		type = "J2PlaneStrain";

		E = 2.0E7;
		nu = 0.49;
          sigmaY = 5.E4;
          H = 1.E7;
	};
};

#solver =
#{
	#type = "LinearSolver";
#};

solver =
{
  type = "NonlinearSolver";
  maxCycle = 10;
  iterMax = 100;
  tol = 1.0e-5;
};

outlabel=['stresses',"equivalent_strain"];

outputModules = ["output","vtk"];

output =
{
	type = "OutputWriter";
	
	onScreen = true;
};

vtk =
{
  type = "MeshWriter";
};