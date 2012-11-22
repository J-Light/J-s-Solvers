volScalarField wallHeatFlux
(
    IOobject
    (
        "wallHeatFlux",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


surfaceScalarField heatFlux
(
    IOobject
    (
        "heatFlux",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("heatFlux", wallHeatFlux.dimensions(), 0.0)
);

