/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    wallHeatFlux

Description
    Calculates and writes the heat flux for all patches as the boundary field
    of a volScalarField and also prints the integrated flux for all wall
    patches.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallFvPatch.H"
#include "OFstream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
#   include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
    
    fileName outputFile(runTime.path()+"/heatFlux.txt");
    OFstream heatflxFile(outputFile);
    
    

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();
        Info<<SMALL;
		
		if (runTime.value() > SMALL)
		{
			#include "createFields.H"
			
			Info<< "\nWall heat fluxes [W]" << endl;
			heatflxFile	<< "Wall heat fluxes [W] at time " 
						<< runTime.timeName()
						<< endl;
			forAll(wallHeatFlux.boundaryField(), patchi)
			{
				dimensionedScalar tHeatFlux= gSum
						   (
							   mesh.magSf().boundaryField()[patchi]
							  *wallHeatFlux.boundaryField()[patchi]
						   );
						 
				if (isA<wallFvPatch>(mesh.boundary()[patchi]))
				{
					Info<< mesh.boundary()[patchi].name()
						<< ":"
						<< tHeatFlux.value()
						<< endl;
					heatflxFile<< mesh.boundary()[patchi].name()
						<< ":"
						<< tHeatFlux.value()
						<< endl;
				}
			}
			Info<< endl;
		}
    }

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
