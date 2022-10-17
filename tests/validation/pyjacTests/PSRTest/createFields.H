    if (mesh.nCells() != 1)
    {
        FatalErrorInFunction
            << "Solver only applicable to single cell cases"
            << exit(FatalError);
    }

    #include "createBaseFields.H"

    autoPtr<fluidMulticomponentThermo> pThermo(fluidMulticomponentThermo::New(mesh));
    fluidMulticomponentThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h");

    autoPtr<basicChemistryModel> pChemistry(basicChemistryModel::New(thermo));

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

    volScalarField Rspecific
    (
        IOobject
        (
            "Rspecific",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            dimensionSet(dimEnergy/dimMass/dimTemperature),
            0
        )
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity, Zero)
    );

    #include "createPhi.H"