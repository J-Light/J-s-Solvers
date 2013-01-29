dictionary solnDict = mesh.solutionDict().subDict("RHOCONTROL");
const dictionary residualDict(solnDict.subOrEmptyDict("residualControl"));

dimensionedScalar HeatFluxResidualRequired
(
	"HeatFluxResidual", 
	dimless, 
	0	
);

HeatFluxResidualRequired.value() = 
	readScalar(residualDict.lookup("wallHeatFlux"));

//Info<<"\nThe Residual value is: "<<HeatFluxResidual.value();
