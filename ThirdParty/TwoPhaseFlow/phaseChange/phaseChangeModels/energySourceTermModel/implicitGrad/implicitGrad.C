/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
	unofficial extension to OpenFOAM.
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

\*---------------------------------------------------------------------------*/


#include "implicitGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "tmp.H"
#include "volFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace TwoPhaseFlow
{
    defineTypeNameAndDebug(implicitGrad, 0);
    addToRunTimeSelectionTable(energySourceTermModel,implicitGrad, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TwoPhaseFlow::implicitGrad::implicitGrad
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const compressibleInterPhaseTransportModel& turbModel,
    const volScalarField& p,
    singleComponentSatProp& satModel,
    geometricVofExt::SimPLIC::reconstruction& reconstructor,
    const dictionary& dict
)
:
    energySourceTermModel
    (
        typeName,
        phase1,
        phase2,
        turbModel,
        p,
        satModel,
        reconstructor,
        dict
    ),
    impDiffFlux_(phase1.mesh())
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix> Foam::TwoPhaseFlow::implicitGrad::TSource1()
{
    reconstructor_.reconstruct(false);
    const volScalarField& TSat = satModel_.TSat();

    tmp<fvScalarMatrix> TSource1
    (
        impDiffFlux_.diffusiveFlux
        (
            reconstructor_.interfaceCell(),
            reconstructor_.interfaceS(),
            reconstructor_.interfaceC(),
            phase1_.thermo().T(),
            phase1_.thermo().kappaEff(turbModel_.alphat()),
            TSat,
            false,
            dimPower
        )
    );

    return TSource1;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::TwoPhaseFlow::implicitGrad::TSource2()
{
    reconstructor_.reconstruct(false);
    const volScalarField& TSat = satModel_.TSat();

    tmp<fvScalarMatrix> TSource2
    (
        impDiffFlux_.diffusiveFlux
        (
            reconstructor_.interfaceCell(),
            reconstructor_.interfaceS(),
            reconstructor_.interfaceC(),
            phase2_.thermo().T(),
            phase2_.thermo().kappaEff(turbModel_.alphat()),
            TSat,
            true,
            dimPower
        )
    );

    return TSource2;
}


Foam::tmp<Foam::volScalarField> Foam::TwoPhaseFlow::implicitGrad::energySource()
{
    reconstructor_.reconstruct(false);
    tmp<volScalarField> energySource
    (
        energyFlux1() + energyFlux2()
    );

    volScalarField& energySourceRef = energySource.ref();
    energySourceRef.ref() *= mag(reconstructor_.interfaceS().internalField())/phase1_.mesh().V();
    energySourceRef.correctBoundaryConditions();

    return energySource;

}


Foam::tmp<Foam::volScalarField> Foam::TwoPhaseFlow::implicitGrad::energyFlux1()
{
    const volScalarField& TSat = satModel_.TSat();
    reconstructor_.reconstruct(false);

    tmp<volScalarField> energyFlux1
    (
        impDiffFlux_.diffusiveFlux
        (
            reconstructor_.interfaceCell(),
            reconstructor_.interfaceS(),
            reconstructor_.interfaceC(),
            phase1_.thermo().T(),
            phase1_.thermo().kappaEff(turbModel_.alphat()),
            TSat,
            false
        )
    );

    return energyFlux1;

}


Foam::tmp<Foam::volScalarField> Foam::TwoPhaseFlow::implicitGrad::energyFlux2()
{
    const volScalarField& TSat = satModel_.TSat();
    reconstructor_.reconstruct(false);

    Foam::tmp<Foam::volScalarField> energyFlux2
    (
        impDiffFlux_.diffusiveFlux
        (
            reconstructor_.interfaceCell(),
            reconstructor_.interfaceS(),
            reconstructor_.interfaceC(),
            phase2_.thermo().T(),
            phase2_.thermo().kappaEff(turbModel_.alphat()),
            TSat,
            true
        )
    );

    return energyFlux2;
}
