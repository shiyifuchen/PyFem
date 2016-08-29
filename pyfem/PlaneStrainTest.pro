input = "PlaneStrainTest.dat";

Elem =
{
	type = "SmallStrainContinuum";
	rho = 2000;
	material = 
	{
		type = "PlaneStrain";
#         type = "J2PlaneStrain";
		
		E = 2.0E8;
		nu = 0.3;
#         sigmaY = 1.E4;
#         H = 1.E8;
	};
};

steps=["Step-1"];

Step-1=
{
	type = "Static";
	solver = 
	{
		type = "GeostaticSolver";
	};
};



outlabel = ["stresses"];

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
