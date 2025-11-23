// coupling.H

volVectorField gradP = fvc::grad(p);

forAll(particles, pI)
{
    Bashyal::particle& p = particles[pI];
    vector force = vector::zero;
    vector torque = vector::zero;
    
    // Clear previous forces
    p.clearForceAndTorque();

    // Volume-based force integration
    // F = sum( -gradP * epsilon * V_cell )
    forAll(mesh.C(), celli)
    {
        if (epsilon[celli] > 0.001)
        {
            // Check if cell belongs to this particle (approximate association)
            if (mag(mesh.C()[celli] - p.position()) < p.getDC() * 0.5)
            {
                scalar vol = mesh.V()[celli] * epsilon[celli];
                force += -gradP[celli] * vol;
                
                // TODO: Add viscous forces
            }
        }
    }
    
    // Apply force (Gravity is handled by particle class usually, or add here)
    // Adding gravity explicitly if particle class doesn't
    force += p.mass() * g.value();

    p.applyForce(force, p.position());
}

// Update particles
forAll(particles, pI)
{
    particles[pI].update(runTime.deltaT().value());
}
