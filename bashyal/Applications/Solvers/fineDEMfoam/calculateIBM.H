// calculateIBM.H

// Reset fields
epsilon = dimensionedScalar("epsilon", dimless, 0.0);
Usolid = dimensionedVector("Usolid", dimVelocity, vector::zero);

forAll(particles, pI)
{
    const Bashyal::particle& p = particles[pI];
    const point& center = p.position();
    scalar radius = p.getDC() * 0.5; // Approximation using characteristic diameter

    // Iterate over all cells
    // TODO: Optimize using mesh search or bounding box
    forAll(mesh.C(), celli)
    {
        scalar dist = mag(mesh.C()[celli] - center);
        
        // Simple sphere check
        if (dist < radius)
        {
            epsilon[celli] = 1.0;
            
            // Calculate solid velocity at cell center
            // u_solid = v + omega x r
            vector r = mesh.C()[celli] - center;
            Usolid[celli] = p.velocity() + (p.angularVelocity() ^ r);
        }
    }
}
