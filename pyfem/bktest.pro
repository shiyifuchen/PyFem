input = "OneElem.dat";

Elem =
{
  type = "SmallStrainContinuum3D";
#  rho = 0;
  material = 
  {
    type = "DuncanChangEv";
    KK=1000;
    n=0.5;
    Rf=0.8;
    c=10;
    Phi0=30;
    GG=0.3;
    D=0.0;
    F=0.0;
    Kur=1500;
    Pa=100;
    DPhi=0.0;  
  };
};

steps=["Step-1","Step-2"];
Step-1=
{
	type = "Static";
	solver = 
	{
		type = "GeostaticSolver";
	};
};

Step-2=
{
	type = "Static";
	solver = 
	{
		type = "NonlinearSolver";
		maxCycle = 1;
		iterMax = 10;
		tol = 1.0e-6;
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