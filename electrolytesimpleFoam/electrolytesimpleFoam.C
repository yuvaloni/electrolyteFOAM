/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    icoFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

    \heading Solver details
    The solver uses the PISO algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}}
          + \div \left( \vec{U} \vec{U} \right)
          - \div \left(\nu \grad \vec{U} \right)
          = - \grad p
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, laminar flow"
        " of Newtonian fluids."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    dimensionedScalar T
    (
        dimensionSet(0, 0, -1, 0, 0, 0, 0),
        1
    );

    dimensionedScalar J
    (
        dimensionSet(0, 2, -1, 0, 0, 0, 0),
        1
    );

    dimensionedScalar Tinv
    (
        dimensionSet(0, 0, 1, 0, 0, 0, 0),
        1
    );
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        while(pimple.loop())
        {
            v1.storePrevIter();
            v2.storePrevIter();
            n1.storePrevIter();
            n2.storePrevIter();
            psi.storePrevIter();
            U.storePrevIter();
            p.storePrevIter();

            #include "CourantNo.H"
            #include "vCourantNo.H"

            // Momentum predictor

            fvVectorMatrix UEqn
            (
                 (1/J) * rho * fvm::ddt(U)
                + rho * fvm::div(phi, U)
                - fvm::laplacian(eta, U)
                - fvm::Sp(-(n1*zeta1+n2*zeta2),U)
                - n1*v1 - n2*v2
                == Tinv*fvOptions(U)
            );
            fvOptions.constrain(UEqn);

            if (pimple.momentumPredictor())
            {
                solve(UEqn == -fvc::grad(p));
                fvOptions.correct(U);
            }

            // --- PISO loop
            while (pimple.correct())
            {
                volScalarField rAU(1.0/UEqn.A());
                volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
                surfaceScalarField phiHbyA
                (
                    "phiHbyA",
                    fvc::flux(HbyA)*T
                + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
                );

                adjustPhi(phiHbyA, U, p);

                // Update the pressure BCs to ensure flux consistency
                constrainPressure(p, U, phiHbyA, rAU);

                // Non-orthogonal pressure corrector loop
                while (pimple.correctNonOrthogonal())
                {
                    // Pressure corrector

                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(T*rAU, p) == fvc::div(phiHbyA)
                    );

                    pEqn.setReference(pRefCell, pRefValue);

                    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

                    if (pimple.finalNonOrthogonalIter())
                    {
                        phi = Tinv*phiHbyA - Tinv*pEqn.flux();
                    }
                }

                #include "continuityErrs.H"

                U = HbyA - rAU*fvc::grad(p);
                U.correctBoundaryConditions();
                fvOptions.correct(U);
            }

                U.relax();
                p.relax();
                solve(fvm::laplacian(psi) == -I*(z1*n1 + z2*n2));
                psi.relax();
                v1 = zeta1*U - fvc::grad(Foam::log(Foam::mag(n1)) + z1*psi);
                v1.correctBoundaryConditions();
                v1.relax();
                v2 = zeta2*U - fvc::grad(Foam::log(Foam::mag(n2)) + z2*psi);
                v2.correctBoundaryConditions();
                v2.relax();
                solve(fvm::ddt(n1) + J*(1/zeta1)*fvm::div(fvc::flux(v1), n1) - J*fvm::laplacian(D1, n1));
                n1.relax();
                solve(fvm::ddt(n2) + J*(1/zeta2)*fvm::div(fvc::flux(v2), n2) - J*fvm::laplacian(D2, n2));
                n2.relax(); 
        }

        runTime.write();

        runTime.printExecutionTime(Info);

    }
    

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
