/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Global
    CourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar v1CoNum = 0.0;
scalar v1meanCoNum = 0.0;
scalar v2CoNum = 0.0;
scalar v2meanCoNum = 0.0;

{
    scalarField v1sumPhi
    (
        fvc::surfaceSum(mag(fvc::flux(v1)))().primitiveField()
    );

    scalarField v2sumPhi
    (
        fvc::surfaceSum(mag(fvc::flux(v2)))().primitiveField()
    );

    v1CoNum = 0.5*gMax(v1sumPhi/mesh.V().field())*runTime.deltaTValue();
    v2CoNum = 0.5*gMax(v2sumPhi/mesh.V().field())*runTime.deltaTValue();

    v1meanCoNum =
        0.5*(gSum(v1sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
            v2meanCoNum =
        0.5*(gSum(v2sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "v1: Courant Number mean: " << v1meanCoNum
    << " max: " << v1CoNum << endl;
Info<< "v2: Courant Number mean: " << v2meanCoNum
    << " max: " << v2CoNum << endl;

// ************************************************************************* //
