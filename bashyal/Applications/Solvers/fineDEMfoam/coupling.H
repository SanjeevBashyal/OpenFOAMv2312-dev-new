// coupling.H

// Construct gradP from p_rgh and rho*g
// p = p_rgh + rho*g.h
// grad(p) approx grad(p_rgh) + rho*g
volVectorField gradP = fvc::grad(p_rgh) + rho*g;

forAll(particles, pI)
{
    Bashyal::particle& p = particles[pI];
    vector force = vector::zero;
    vector torque = vector::zero;
    
    // Clear previous forces
    p.clearForceAndTorque();

    // Pre-calculate viscous stress term (laplacian of U)
    volVectorField lapUField = fvc::laplacian(U);

    vector fBuoy = vector::zero;
    vector fVisc = vector::zero;
    scalar totVol = 0.0;

    // Volume-based force integration
    // F = sum( -gradP * epsilon * V_cell )
    forAll(mesh.C(), celli)
    {
        if (epsilon[celli] > 0.001)
        {
            scalar dist = mag(mesh.C()[celli] - p.position());
            scalar radius = p.getDC() * 0.5;
            // Integrate over the entire support region (epsilon > 0)
            if (dist < radius * 2.5)
            {
                scalar vol = mesh.V()[celli] * epsilon[celli];
                totVol += vol;
                
                // Pressure force (Buoyancy)
                vector fb = -gradP[celli] * vol;
                force += fb;
                fBuoy += fb;
                
                // Viscous force (Drag)
                scalar mu = 0.001; // Water viscosity
                vector lapU = lapUField[celli];
                vector fv = mu * lapU * vol;
                force += fv;
                fVisc += fv;
            }
        }
    }

    // Apply force (Gravity is handled by particle class usually, or add here)
    // Adding gravity explicitly if particle class doesn't
    vector fGrav = p.mass() * g.value();
    force += fGrav;

    if (pI == 0)
    {
        Info << "Particle 0 Debug:" << endl;
        Info << "  Vol Integrated: " << totVol << " (Expected: " << pow3(p.getDC()) << ")" << endl;
        Info << "  Buoyancy: " << fBuoy << endl;
        Info << "  Viscous: " << fVisc << endl;
        Info << "  Gravity: " << fGrav << endl;
        Info << "  Total: " << force << endl;
    }

    p.applyForce(force, p.position());
}

// Update particles
forAll(particles, pI)
{
    particles[pI].update(runTime.deltaT().value());
}
