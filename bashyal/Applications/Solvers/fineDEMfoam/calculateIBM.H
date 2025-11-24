// calculateIBM.H

// Reset fields
epsilon = dimensionedScalar("epsilon", dimless, 0.0);
Usolid = dimensionedVector("Usolid", dimVelocity, vector::zero);

forAll(particles, pI)
{
    const Bashyal::particle& p = particles[pI];
    const point& center = p.position();
    scalar radius = p.getDC() * 0.5;
    scalar supportRadius = radius * 1.05; // Reduced to 1.05 for sharp interface
    scalar particleVolume = 4.0/3.0 * Foam::constant::mathematical::pi * pow3(radius);
    // For cube:
    if (abs(p.getDC() - 0.002) < 1e-6 || abs(p.getDC() - 0.003) < 1e-6 || abs(p.getDC() - 0.1) < 1e-6) 
    {
        particleVolume = pow3(p.getDC());
    }

    scalar totalWeight = 0.0;
    DynamicList<label> affectedCells;
    DynamicList<scalar> weights;

    // Iterate over all cells (optimize later)
    forAll(mesh.C(), celli)
    {
        scalar dist = mag(mesh.C()[celli] - center);
        
        if (dist < supportRadius)
        {
            // Kernel function (simple linear)
            scalar w = (1.0 - dist/supportRadius);
            
            affectedCells.append(celli);
            weights.append(w);
            totalWeight += w * mesh.V()[celli];
        }
    }

    // Normalize and assign epsilon
    if (totalWeight > SMALL)
    {
        scalar scale = particleVolume / totalWeight;
        forAll(affectedCells, i)
        {
            label celli = affectedCells[i];
            epsilon[celli] += weights[i] * scale; // Additive for multiple particles
            
            // Cap at 1.0 - REMOVED to conserve volume for force calculation
            // if (epsilon[celli] > 1.0) epsilon[celli] = 1.0;
            
            // Calculate solid velocity
            vector r = mesh.C()[celli] - center;
            Usolid[celli] = p.velocity() + (p.angularVelocity() ^ r);
        }
    }
}
