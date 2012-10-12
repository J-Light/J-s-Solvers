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
    rhoCentralFoam

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"
#include "wallFvPatch.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedRhoFvPatchScalarField.H"
#include "boundMinMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "readThermophysicalProperties.H"
    #include "readTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"
  
    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
    dimensionedScalar smallU("smallU", dimVelocity, ROOTVSMALL);
    dimensionedScalar gamma("gamma",dimless, 1.4);
    dimensionedScalar oldQdot("oldQdot",dimless, 0);
	dimensionedScalar NormFactHeatFlux("NormFactHeatFlux",dimless, 0);
    volScalarField eOld = e;

    Info<< "\nStarting time loop\n" << endl;
    

    while (runTime.run())
    {
        // --- upwind interpolation of primitive fields on faces

        #include "readFieldBounds.H"
        
        surfaceScalarField rho_pos
        (
            fvc::interpolate(rho, pos, "reconstruct(rho)")
        );
        surfaceScalarField rho_neg
        (
            fvc::interpolate(rho, neg, "reconstruct(rho)")
        );

        surfaceVectorField rhoU_pos
        (
            fvc::interpolate(rhoU, pos, "reconstruct(U)")
        );
        surfaceVectorField rhoU_neg
        (
            fvc::interpolate(rhoU, neg, "reconstruct(U)")
        );

        volScalarField rPsi(1.0/psi);
        surfaceScalarField rPsi_pos
        (
            fvc::interpolate(rPsi, pos, "reconstruct(T)")
        );
        surfaceScalarField rPsi_neg
        (
            fvc::interpolate(rPsi, neg, "reconstruct(T)")
        );

        surfaceScalarField e_pos
        (
            fvc::interpolate(e, pos, "reconstruct(T)")
        );
        surfaceScalarField e_neg
        (
            fvc::interpolate(e, neg, "reconstruct(T)")
        );

        surfaceVectorField U_pos(rhoU_pos/rho_pos);
        surfaceVectorField U_neg(rhoU_neg/rho_neg);

        surfaceScalarField p_pos(rho_pos*rPsi_pos);
        surfaceScalarField p_neg(rho_neg*rPsi_neg);

        surfaceScalarField phiv_pos(U_pos & mesh.Sf());
        surfaceScalarField phiv_neg(U_neg & mesh.Sf());

        volScalarField c(sqrt(thermo.Cp()/thermo.Cv()*rPsi));
        surfaceScalarField cSf_pos
        (
            fvc::interpolate(c, pos, "reconstruct(T)")*mesh.magSf()
        );
        surfaceScalarField cSf_neg
        (
            fvc::interpolate(c, neg, "reconstruct(T)")*mesh.magSf()
        );

        surfaceScalarField ap
        (
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );
        surfaceScalarField am
        (
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );

        surfaceScalarField a_pos(ap/(ap - am));

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        surfaceScalarField aSf(am*a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg(1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos(phiv_pos - aSf);
        surfaceScalarField aphiv_neg(phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        #include "compressibleCourantNo.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

        surfaceVectorField phiUp
        (
            (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
          + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf()
        );

        surfaceScalarField phiEp
        (
            aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
          + aSf*p_pos - aSf*p_neg
        );

        volScalarField muEff(turbulence->muEff());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));
        
        //Info<<"\nMin Value of rho inviscid:"<<min(rho).value()<<"  "<<max(rho).value()<<'\n';
        boundMinMax(rho,rhoMin, rhoMax);

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));
        
        //Info<<"\nMin Value of RhoUmag inviscid:"<<min(sqrt(magSqr(rhoU))).value()<<"  "<<max(sqrt(magSqr(rhoU))).value()<<'\n';

        U.dimensionedInternalField() =
            rhoU.dimensionedInternalField()
           /rho.dimensionedInternalField();
        U.correctBoundaryConditions();
        
        volScalarField magU = mag(U);
        
        if (max(magU) > UMax)
        {
        	Info<< "bounding " << U.name()
            	<< " max: " << max(magU).value()
                << endl;

			volScalarField UDiff = magU - UMax;
			volScalarField Ulimiter = ::pos(UDiff) * UMax / (magU + smallU) + ::neg(UDiff);
            Ulimiter.max(scalar(0));
            Ulimiter.min(scalar(1));

            U *= Ulimiter;
            U.correctBoundaryConditions();
            rhoU = rho*U;
		}
        
        rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();
        magU = mag(U);
        //Info<<"\nMin Value of Umag inviscid:"<<min(magU).value()<<"  "<<max(magU).value()<<'\n';
                
        volScalarField rhoBydt(rho/runTime.deltaT());

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tauMC)
            );
            rhoU = rho*U;
        }

        // --- Solve energy
        surfaceScalarField sigmaDotU
        (
            (
                fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
              + (mesh.Sf() & fvc::interpolate(tauMC))
            )
            & (a_pos*U_pos + a_neg*U_neg)
        );

        solve
        (
            fvm::ddt(rhoE)
          + fvc::div(phiEp)
          - fvc::div(sigmaDotU)
        );
        
      //  Info<<"\nMin Value of rhoE inviscid:"<<min(rhoE).value()<<'\n';
        boundMinMax(rhoE, rhoEMin, rhoEMax);
        
        e = rhoE/rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();
        
      //  Info<<"\nMin Value of e inviscid:"<<min(e).value()<<'\n';
        boundMinMax(e, eMin, eMax);
        
      //  Info<<"\nMin Value of T inviscid:"<<min(T).value()<<'\n';
        thermo.correct();
        
        rhoE.boundaryField() =
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            );

        if (!inviscid)
        {
            volScalarField k("k", thermo.Cp()*muEff/Pr);
            solve
            (
                fvm::ddt(rho, e) - fvc::ddt(rho, e)
              - fvm::laplacian(turbulence->alphaEff(), e)
              + fvc::laplacian(turbulence->alpha(), e)
              - fvc::laplacian(k, T)
            );
           // Info<<"\nMin Value of E viscous:"<<min(e).value()<<'\n';
            boundMinMax(e, eMin, eMax);
            
            thermo.correct();
            rhoE = rho*(e + 0.5*magSqr(U));
            boundMinMax(rhoE, rhoEMin, rhoEMax);
        }

        p.dimensionedInternalField() =
            rho.dimensionedInternalField()
           /psi.dimensionedInternalField();
        p.correctBoundaryConditions();
        rho.boundaryField() = psi.boundaryField()*p.boundaryField();

        turbulence->correct();
        
        surfaceScalarField heatFlux
        (
            fvc::interpolate(turbulence->alphaEff())*gamma*fvc::snGrad(e)
        );

        const surfaceScalarField::GeometricBoundaryField& patchHeatFlux =
            heatFlux.boundaryField();
		
		forAll(patchHeatFlux, patchi)
		{
			if (isA<wallFvPatch>(mesh.boundary()[patchi]))
			{
				dimensionedScalar tHeatFlux= gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchHeatFlux[patchi]
                       );
				Info<<"\nWall heat fluxes [W]: "
					<< tHeatFlux.value() <<endl;
				
				dimensionedScalar 
					resHeatFlux 
					= 	
					mag(tHeatFlux-oldQdot)/runTime.deltaTValue();
					
				Info<<"Heat Flux residual: "
					<<resHeatFlux.value()<<endl;
				oldQdot=tHeatFlux;
			}
			
		}
        //label TimeIdx=runTime.timeIndex();
        //if (TimeIdx< 100)
		//{
			//Info<<"Setting old values";
			//eOld = e;
		//}
		//else
		//{
			//volScalarField eDiff = e-eOld;
			//volScalarField SquardRes = sqr(eDiff/runTime.deltaTValue());
			//dimensionedScalar ResidualE = sqrt(sum(SquardRes)); 
			//if (TimeIdx < 110 )
			//{
				//if (ResidualE > NormResE)
				//{
					//NormResE=ResidualE;
				//}
			//}
			//dimensionedScalar NResE=ResidualE/NormResE;
			//Info<<"The residual is: "<<NResE.value()<< endl;
		//}

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
