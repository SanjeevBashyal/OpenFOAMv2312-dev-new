#ifndef particlesGenerator_H
#define particlesGenerator_H

#include "quickInclude.H"
#include "PtrList.H"
#include "point.H"
#include "vector.H"
#include "tensor.H"
#include "Random.H"

#include "boundary.H"
#include "faceList.H"
#include "particle.H"

namespace Bashyal
{
    class particlesGenerator
    {
    private:
        Foam::Random rng_;

        Foam::tensor makeDiagonalTensor(const Foam::scalar Ixx, const Foam::scalar Iyy, const Foam::scalar Izz) const;

        // Geometry helpers (standalone; no dependency on aggregate classes)
        Foam::pointField makeCubePoints(const Foam::scalar side) const;
        Foam::faceList makeCubeFaces() const;

        Foam::pointField makeCuboidPoints(const Foam::scalar ax, const Foam::scalar ay, const Foam::scalar az) const;
        Foam::faceList makeCuboidFaces() const;

        Foam::pointField makeSpherePoints(const Foam::scalar radius, const int resolution) const;
        Foam::faceList makeSphereLikeFaces(const int resolution) const; // for sphere and well-rounded

    public:
        particlesGenerator();

        // Create N cubical particles with side length uniformly sampled in [minSize, maxSize]
        Foam::PtrList<Bashyal::particle> createCubical(
            const Foam::label n,
            const Foam::scalar minSize,
            const Foam::scalar maxSize);

        // Create N cuboid particles with dimensions sampled per-axis in [minDims, maxDims]
        Foam::PtrList<Bashyal::particle> createCuboid(
            const Foam::label n,
            const Foam::vector& minDims,
            const Foam::vector& maxDims);

        // Create N well-rounded polyhedral particles
        // roundnessFactor in [0,1], resolution >= 3
        Foam::PtrList<Bashyal::particle> createWellRounded(
            const Foam::label n,
            const Foam::scalar s1,
            const Foam::scalar s2,
            const Foam::scalar roundnessFactor,
            const int resolution);

        // Create N spherical particles
        Foam::PtrList<Bashyal::particle> createSpherical(
            const Foam::label n,
            const Foam::scalar radius,
            const int resolution);
    };
}

#endif


