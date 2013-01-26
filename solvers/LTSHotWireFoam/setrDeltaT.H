{

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.dimensionedInternalField() = max
    (
      1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
      fvc::surfaceSum(amaxSf)().dimensionedInternalField()
              /(maxCo*mesh.V())
    );

    Info<< "    Flow        = "
            << gMin(1/rDeltaT.internalField()) << ", "
            << gMax(1/rDeltaT.internalField()) << endl;

    // Limit the largest time scale
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rDeltaT.max(1/maxDeltaT);


    // Spatially smooth the time scale field
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar rDeltaTSmoothingCoeff = 0.1;
    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

    Info<< "    Overall     = " << min(1/rDeltaT).value()
        << ", " << max(1/rDeltaT).value() << nl << endl;

}
