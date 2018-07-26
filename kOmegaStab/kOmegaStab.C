/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "kOmegaStab.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaStab<BasicTurbulenceModel>::correctNut()
{
  this->nut_ = k_/max(omega_,lambda2_*beta_/(Cmu_*gamma_)*2*magSqr(symm(fvc::grad(this->U_)))/(2*magSqr(skew(fvc::grad(this->U_)))+pOmegaSmall_)*omega_); //nut as suggested by Larsen and (Fuhrman 2018)
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaStab<BasicTurbulenceModel>::kOmegaStab
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta",
            this->coeffDict_,
            0.072
        )
    ),
    gamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            0.52
        )
    ),
    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            0.5
        )
    ),
    alphaBS_ //Coefficient for the buoyancy production
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaBS",
            this->coeffDict_,
            1.36
        )
    ),
    lambda2_ //New limiter suggested by Larsen and Fuhrman (2018)
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "lambda2",
            this->coeffDict_,
            0.05
        )
    ),
    pOmegaSmall_("pOmegaSmall", dimless/(dimTime*dimTime), SMALL),
pOmega_
(
            IOobject
            (
                    "pOmega",
                   this-> runTime_.timeName(),
                   this-> mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
            ),
            this->mesh_,
	    dimensionedScalar("dime",dimensionSet(0, 0, 0, 0, 0), 1.0) 
    ),


    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
     ),
    gField_ //Needed for the buoyancy production term
(
IOobject
(
"g",
this->db().time().constant(),
this->db(),
IOobject::MUST_READ,
IOobject::NO_WRITE
)
 )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaStab<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        beta_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
	alphaBS_.readIfPresent(this->coeffDict());
	lambda2_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kOmegaStab<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }


    // Local references
    const alphaField& alpha = this->alpha_;
  const volVectorField& U = this->U_;
    const volScalarField& rho1 = U.db().objectRegistry::lookupObject<volScalarField>("rho");
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
  

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Include density and buoyancy production term. Implemented by Bjarke Eltard Larsen 18-07-2018

    // Calculate the Brunt-Vaisala frequency
 volScalarField N2 = gField_&fvc::grad(rho1)/rho1;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Calculate eddy viscosity using the stabilizing approach described in Larsen and Fuhrman (2018)

 volScalarField p0 = 2.0*magSqr(symm(fvc::grad(U))); // 2S_ij S_ij
 volScalarField pOmega= 2.0*magSqr(skew(fvc::grad(U))); // 2 Omega_ij Omega_ij

 volScalarField omegaTilde2 = max (omega_+this->omegaMin_, lambda2_*(beta_/(Cmu_*gamma_))*p0*omega_/(pOmega+pOmegaSmall_));
 pOmega_=lambda2_;
  volScalarField nut = k_/omegaTilde2;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G
    (
        this->GName(),
        nut*(tgradU() && dev(twoSymm(tgradU())))
    );
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        gamma_*alpha*rho*p0
	   - fvm::SuSp(((2.0/3.0)*gamma_)*alpha*rho*divU, omega_)
      - fvm::Sp(beta_*alpha*rho*omega_, omega_)
      + fvOptions(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G
	   - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_) 
      - fvm::Sp(Cmu_*alpha*rho*omega_, k_)
 - fvm::Sp(nut*alphaBS_*N2/max(k_,this->kMin_),k_) // Buoyancy production term
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);


     correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
