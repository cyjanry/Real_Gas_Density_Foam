    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();

    volScalarField& h = thermo.h();
//     volScalarField& e = thermo.e();

    volScalarField& p = thermo.p();
     const volScalarField& T = thermo.T();
//     const volScalarField& psi = thermo.psi();





    Info<< "Allocating field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );


// Add Acoustic speed, Jianhui qi
Info<< "Allocating field A\n" << endl;
    volScalarField A
    (
        IOobject
        (
            "A",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
       ),
       sqrt(thermo.Cp()/thermo.Cv()*p/rho)
    );


    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    );

#   include "compressibleCreatePhi.H"

    Info<< "Creating MRF model\n" << endl;
    MRFZones mrfZones(mesh);

    volVectorField URel
    (
        IOobject
        (
            "URel",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U
    );
    URel.correctBoundaryConditions();
    mrfZones.relativeVelocity(U,URel);

    Info<< "Creating turbulence model\n" << endl;
    // all turbulence models which use the skew of this velocity gradient will
    // fail for rotating meshs, like SpalartAllmaras or kOmegaSST_LowRe
    // Otherwise if using URel, then you will become problems at Rotor-Stator
    // Interfaces, as the grad(Urel) will be computed wrong there
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Create Riemann solver\n" << endl;
//   godunovFlux_Real<AUSMplusALERealFlux> Godunov(p, U, rho, thermo, turbulence());
      godunovFlux_Real<hllcALERealFlux> Godunov(p, U, rho, thermo, turbulence());
//    godunovFlux<hllcALEFlux> Godunov(p, U, rho, thermo, turbulence());
//     godunovFlux<AUSMplusALEFlux> Godunov(p, U, rho, thermo, turbulence());
//     godunovFlux<AUSMplusUpALEFlux> Godunov(p, U, rho, thermo, turbulence());
//     godunovFlux<roeALEFlux> Godunov(p, U, rho, thermo, turbulence());
//     godunovFlux<rusanovALEFlux> Godunov(p, U, rho, thermo, turbulence());

    Info<< "Allocating field rhoU\n" << endl;
    volVectorField rhoU
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*U
    );

    Info<< "Allocating field rhoE\n" << endl;
    volScalarField rhoE
    (
        IOobject
        (
            "rhoE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*(h + 0.5*magSqr(U) + turbulence->k()) - p
//         rho*(e + 0.5*magSqr(U) + turbulence->k())
    );

    Info<< "Allocating field e\n" << endl;
    volScalarField e
    (
        IOobject
        (
            "e",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
       rhoE/rho - 0.5*magSqr(U) - turbulence->k()
    );






    localTimeStep_Real localTimeStep_Real(URel, thermo, turbulence());

    Info<< "Allocating physDeltaT list for RK and Dual-Time Stepping\n" << endl;
    IOField<scalar> physDeltaT
    (
        IOobject
        (
            "physDeltaT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        3
    );

    label numberSubCycles = 1;
    scalarList beta(1,1.0);
    bool secondOrder = true;

#   include "readMultiStage.H"
#   include "createDualTimeSteppingFields.H"
