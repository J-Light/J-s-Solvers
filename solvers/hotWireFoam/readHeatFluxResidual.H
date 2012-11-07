dictionary solnDict = mesh.solutionDict().subDict("RHOCONTROL");
const dictionary residualDict(solnDict.subOrEmptyDict("residualControl"));

dimensionedScalar HeatFluxResidual
(
	"HeatFluxResidual", 
	dimless, 
	0	
);

HeatFluxResidual.value() = readScalar(residualDict.lookup("wallHeatFlux"));

//Info<<"\nThe Residual value is: "<<HeatFluxResidual.value();
