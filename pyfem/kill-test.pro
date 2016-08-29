input = "kill-test.dat";

Elem1 =
{
	type = "SmallStrainContinuum";
	
	rho = 2000;
	
	material = 
	{
		type = "J2PlaneStrain";
#		type = "PlaneStrain";

		E = 2.0E8;
		nu = 0.3;
         sigmaY = 9.0E4;
         H = 2.0E8;
	};
};

Elem2 =
{
	type = "SmallStrainContinuum";
	
	rho = 2000;
	
	material = 
	{
		type = "J2PlaneStrain";
#		type = "PlaneStrain";

		E = 2.0E8;
		nu = 0.3;
         sigmaY = 9.0E4;
         H = 2.0E8;
	};
};

#steps=["Step-1","Step-2"];
steps=["Geo-Step"];
Geo-Step=
{
	type = "Static";
	solver = 
	{
		type = "GeostaticSolver";
	};
};

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

Step-2=
{
	type = "Static";
	solver = 
	{
		type = "NonlinearSolver";
		maxCycle = 1;
		iterMax = 100;
		tol = 1.0e-6;
	};
	kill = ["Elem1"];
};

#outlabel=['stresses','equivalent_strain'];
outlabel=['stresses'];

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