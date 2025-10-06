/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 24 "/home/yufan/xd/solvers/CFD/10-phaseField/testCase/0/U/#codeStream"
#include "fvCFD.H"
        #include "constants.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    void codeStream_2f852221fe1165d00a9aebc900bd46669dc58e91
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 37 "/home/yufan/xd/solvers/CFD/10-phaseField/testCase/0/U/#codeStream"
const IOdictionary& d = static_cast<const IOdictionary&>(dict);
        const fvMesh& mesh = refCast<const fvMesh>(d.db());

        const vectorField& points = mesh.C();
        vectorField u(mesh.nCells());
        const double t = mesh.time().value();

        const scalar twoPi = 2.0*Foam::constant::mathematical::pi;
        const scalar pi = Foam::constant::mathematical::pi;
        forAll (points, celli)
        {
            scalar x = points[celli].x();
            scalar y = points[celli].y();
            u[celli] =
                vector
                (
                  - Foam::sqr(Foam::sin(pi*x))*Foam::sin(twoPi*y)*Foam::cos(pi*t/4.0),
                    Foam::sqr(Foam::sin(pi*y))*Foam::sin(twoPi*x)*Foam::cos(pi*t/4.0),
                    0.0
                );
        }
        os  << "nonuniform " << u;
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

