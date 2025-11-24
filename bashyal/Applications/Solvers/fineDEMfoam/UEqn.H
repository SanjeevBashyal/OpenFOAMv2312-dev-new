MRF.correctBoundaryVelocity(U);

    // Penalty Method for IBM
    // Enforce U = Usolid where epsilon > 0
    // Penalty coefficient should be large enough to dominate other terms
    // Dimensions: rho [1 -3 0] * rate [0 0 -1] = [1 -3 -1]
    // Resulting term Sp(C, U) -> [1 -3 -1] * [0 1 -1] = [1 -2 -2] (Force/Vol)
    volScalarField penaltyCoeff = epsilon * rho * dimensionedScalar("penalty", dimless/dimTime, 1.0e5);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
      + fvm::Sp(penaltyCoeff, U)
     ==
        fvOptions(rho, U)
      + penaltyCoeff * Usolid
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
