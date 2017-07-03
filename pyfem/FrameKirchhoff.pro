input = "Frame2.dat";

BeamElem =
{
  type = "TimoshenkoBeam";
#  type = "KirchhoffBeam";
  E    = 7.2e6;
#  A    = 6.0;
#  I    = 2.0;
#  G    = 5e6;
  A    = 1.0;
  I    = 0.08333;
  G    = 2.77e6;
};

steps=["Step-1"];

Step-1=
{
	type = "Static";
	solver = 
	{
		type = "LinearSolver";
	};
};
outputModules = ["output"];

output =
{
  type = "OutputWriter";

  onScreen = true;
};