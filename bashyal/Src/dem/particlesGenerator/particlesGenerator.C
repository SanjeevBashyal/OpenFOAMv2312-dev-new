#include "particlesGenerator.H"

#include "face.H"
#include "boundBox.H"
#include "mathematicalConstants.H"

namespace Bashyal
{
    particlesGenerator::particlesGenerator()
    : rng_(12345)
    {
    }

    Foam::tensor particlesGenerator::makeDiagonalTensor(const Foam::scalar Ixx, const Foam::scalar Iyy, const Foam::scalar Izz) const
    {
        return Foam::tensor(
            Ixx, 0,   0,
            0,   Iyy, 0,
            0,   0,   Izz
        );
    }

    // Create a particle from raw points/faces by building a boundary directly
    static Bashyal::particle* makeParticleFromGeom(
        const Foam::pointField& localPoints,
        const Foam::faceList& faces,
        const Foam::scalar mass,
        const Foam::tensor& moi,
        const Foam::point& initialPosition)
    {
        Bashyal::boundary b(localPoints, faces);
        return new Bashyal::particle(b, mass, moi, initialPosition);
    }

    // Geometry helpers
    Foam::pointField particlesGenerator::makeCubePoints(const Foam::scalar side) const
    {
        Foam::pointField nodes(8);
        Foam::scalar s = side;
        nodes[0] = Foam::point(-s/2, -s/2, -s/2);
        nodes[1] = Foam::point(-s/2,  s/2, -s/2);
        nodes[2] = Foam::point( s/2,  s/2, -s/2);
        nodes[3] = Foam::point( s/2, -s/2, -s/2);
        nodes[4] = Foam::point(-s/2, -s/2,  s/2);
        nodes[5] = Foam::point(-s/2,  s/2,  s/2);
        nodes[6] = Foam::point( s/2,  s/2,  s/2);
        nodes[7] = Foam::point( s/2, -s/2,  s/2);
        return nodes;
    }

    Foam::faceList particlesGenerator::makeCubeFaces() const
    {
        Foam::faceList faces(6);
        faces[0] = Foam::face({0, 1, 2, 3});
        faces[1] = Foam::face({4, 7, 6, 5});
        faces[2] = Foam::face({0, 4, 5, 1});
        faces[3] = Foam::face({3, 2, 6, 7});
        faces[4] = Foam::face({0, 3, 7, 4});
        faces[5] = Foam::face({1, 5, 6, 2});
        return faces;
    }

    Foam::pointField particlesGenerator::makeCuboidPoints(const Foam::scalar ax, const Foam::scalar ay, const Foam::scalar az) const
    {
        Foam::pointField nodes(8);
        nodes[0] = Foam::point(-ax/2, -ay/2, -az/2);
        nodes[1] = Foam::point(-ax/2,  ay/2, -az/2);
        nodes[2] = Foam::point( ax/2,  ay/2, -az/2);
        nodes[3] = Foam::point( ax/2, -ay/2, -az/2);
        nodes[4] = Foam::point(-ax/2, -ay/2,  az/2);
        nodes[5] = Foam::point(-ax/2,  ay/2,  az/2);
        nodes[6] = Foam::point( ax/2,  ay/2,  az/2);
        nodes[7] = Foam::point( ax/2, -ay/2,  az/2);
        return nodes;
    }

    Foam::faceList particlesGenerator::makeCuboidFaces() const
    {
        return makeCubeFaces();
    }

    Foam::pointField particlesGenerator::makeSpherePoints(const Foam::scalar radius, const int resolution) const
    {
        const int pointsPerRing = resolution;
        const int totalPoints = (resolution - 1) * pointsPerRing + 2;
        Foam::pointField points(totalPoints);

        // North pole
        points[0] = Foam::point(0, 0, radius);

        // Intermediate rings
        for (int i = 1; i < resolution; ++i)
        {
            Foam::scalar theta = Foam::constant::mathematical::pi * i / resolution;
            Foam::scalar sinTheta = Foam::sin(theta);
            Foam::scalar cosTheta = Foam::cos(theta);

            for (int j = 0; j < resolution; ++j)
            {
                Foam::scalar phi = 2 * Foam::constant::mathematical::pi * j / resolution;
                Foam::scalar x = radius * sinTheta * Foam::cos(phi);
                Foam::scalar y = radius * sinTheta * Foam::sin(phi);
                Foam::scalar z = radius * cosTheta;
                int index = (i - 1) * pointsPerRing + j + 1;
                points[index] = Foam::point(x, y, z);
            }
        }

        // South pole
        points[totalPoints - 1] = Foam::point(0, 0, -radius);

        return points;
    }

    Foam::faceList particlesGenerator::makeSphereLikeFaces(const int resolution) const
    {
        Foam::faceList faces;
        int pointsPerRing = resolution;

        // Top triangular faces (north pole)
        for (int j = 0; j < resolution; ++j)
        {
            int p1 = 0;                         // North pole
            int p2 = j + 1;                     // Current point on first ring
            int p3 = (j + 1) % resolution + 1;  // Next point on first ring
            faces.append(Foam::face({p1, p2, p3}));
        }

        // Intermediate quadrilateral faces
        int offset = 1;
        for (int i = 1; i < resolution - 1; ++i)
        {
            for (int j = 0; j < resolution; ++j)
            {
                int p1 = offset + j;
                int p2 = offset + (j + 1) % resolution;
                int p3 = offset + pointsPerRing + (j + 1) % resolution;
                int p4 = offset + pointsPerRing + j;
                faces.append(Foam::face({p1, p2, p3, p4}));
            }
            offset += pointsPerRing;
        }

        // Bottom triangular faces (south pole)
        int bottomOffset = (resolution - 2) * pointsPerRing + 1;
        int southPole = bottomOffset + resolution;
        for (int j = 0; j < resolution; ++j)
        {
            int p1 = southPole;                 // South pole
            int p2 = bottomOffset + j;          // Current point on last ring
            int p3 = bottomOffset + (j + 1) % resolution; // Next point
            faces.append(Foam::face({p1, p2, p3}));
        }

        return faces;
    }

    Foam::PtrList<Bashyal::particle> particlesGenerator::createCubical(
        const Foam::label n,
        const Foam::scalar minSize,
        const Foam::scalar maxSize)
    {
        Foam::PtrList<Bashyal::particle> particles(n);

        for (Foam::label i = 0; i < n; ++i)
        {
            Foam::scalar r = rng_.sample01<Foam::scalar>();
            Foam::scalar side = minSize + (maxSize - minSize) * r;
            Foam::pointField nodes = makeCubePoints(side);
            Foam::faceList faces = makeCubeFaces();

            // Mass and inertia for a solid cube of side length 'side' and density=1
            Foam::scalar density = 1.0;
            Foam::scalar mass = density * Foam::pow3(side);
            Foam::scalar I = (mass * side * side) / 6.0; // Ixx = Iyy = Izz for cube about center
            Foam::tensor moi = makeDiagonalTensor(I, I, I);

            Foam::point initialPosition(0, 0, 0);
            particles.set(i, makeParticleFromGeom(nodes, faces, mass, moi, initialPosition));
        }

        return particles;
    }

    Foam::PtrList<Bashyal::particle> particlesGenerator::createCuboid(
        const Foam::label n,
        const Foam::vector& minDims,
        const Foam::vector& maxDims)
    {
        Foam::PtrList<Bashyal::particle> particles(n);

        for (Foam::label i = 0; i < n; ++i)
        {
            // Sample dimensions ax, ay, az
            Foam::scalar rx = rng_.sample01<Foam::scalar>();
            Foam::scalar ry = rng_.sample01<Foam::scalar>();
            Foam::scalar rz = rng_.sample01<Foam::scalar>();
            Foam::scalar ax = minDims.x() + (maxDims.x() - minDims.x()) * rx;
            Foam::scalar ay = minDims.y() + (maxDims.y() - minDims.y()) * ry;
            Foam::scalar az = minDims.z() + (maxDims.z() - minDims.z()) * rz;

            Foam::pointField nodes = makeCuboidPoints(ax, ay, az);
            Foam::faceList faces = makeCuboidFaces();

            // Mass and inertia for solid rectangular box ax, ay, az, density=1
            Foam::scalar density = 1.0;
            Foam::scalar mass = density * ax * ay * az;
            Foam::scalar Ixx = (mass / 12.0) * (ay*ay + az*az);
            Foam::scalar Iyy = (mass / 12.0) * (ax*ax + az*az);
            Foam::scalar Izz = (mass / 12.0) * (ax*ax + ay*ay);
            Foam::tensor moi = makeDiagonalTensor(Ixx, Iyy, Izz);

            Foam::point initialPosition(0, 0, 0);
            particles.set(i, makeParticleFromGeom(nodes, faces, mass, moi, initialPosition));
        }

        return particles;
    }

    Foam::PtrList<Bashyal::particle> particlesGenerator::createWellRounded(
        const Foam::label n,
        const Foam::scalar s1,
        const Foam::scalar s2,
        const Foam::scalar roundnessFactor,
        const int resolution)
    {
        Foam::PtrList<Bashyal::particle> particles(n);

        for (Foam::label i = 0; i < n; ++i)
        {
            // sample size in [s1, s2]
            Foam::scalar r = rng_.sample01<Foam::scalar>();
            Foam::scalar s = s1 + (s2 - s1) * r;

            // roundnessFactor in [0,1]; create sphere points scaled randomly per-axis
            Foam::scalar anisotropyScale = 1.0 - roundnessFactor; // 0 -> sphere, 1 -> max deviation
            Foam::scalar sx = 1.0 + anisotropyScale * (2.0 * rng_.sample01<Foam::scalar>() - 1.0);
            Foam::scalar sy = 1.0 + anisotropyScale * (2.0 * rng_.sample01<Foam::scalar>() - 1.0);
            Foam::scalar sz = 1.0 + anisotropyScale * (2.0 * rng_.sample01<Foam::scalar>() - 1.0);

            Foam::pointField points = makeSpherePoints(1.0, resolution);
            forAll(points, pi) { points[pi].x() *= sx; points[pi].y() *= sy; points[pi].z() *= sz; }

            // Scale to target size using bounding box span
            Foam::boundBox bb(points);
            Foam::vector span = bb.span();
            Foam::scalar maxExtent = Foam::max(Foam::max(span.x(), span.y()), span.z());
            Foam::scalar scaleFactor = s / maxExtent;
            forAll(points, pi) { points[pi] *= scaleFactor; }

            Foam::faceList faces = makeSphereLikeFaces(resolution);

            // Estimate ellipsoidal inertia from extents
            Foam::scalar ax = span.x() * scaleFactor;
            Foam::scalar ay = span.y() * scaleFactor;
            Foam::scalar az = span.z() * scaleFactor;

            // Approximate as equivalent ellipsoid with semi-axes aSemi, bSemi, cSemi
            Foam::scalar aSemi = ax/2.0, bSemi = ay/2.0, cSemi = az/2.0;
            Foam::scalar density = 1.0;
            Foam::scalar mass = density * (4.0/3.0) * Foam::constant::mathematical::pi * aSemi * bSemi * cSemi;
            Foam::scalar Ixx = (mass/5.0) * (bSemi*bSemi + cSemi*cSemi);
            Foam::scalar Iyy = (mass/5.0) * (aSemi*aSemi + cSemi*cSemi);
            Foam::scalar Izz = (mass/5.0) * (aSemi*aSemi + bSemi*bSemi);
            Foam::tensor moi = makeDiagonalTensor(Ixx, Iyy, Izz);

            Foam::point initialPosition(0, 0, 0);
            particles.set(i, makeParticleFromGeom(points, faces, mass, moi, initialPosition));
        }

        return particles;
    }

    Foam::PtrList<Bashyal::particle> particlesGenerator::createSpherical(
        const Foam::label n,
        const Foam::scalar radius,
        const int resolution)
    {
        Foam::PtrList<Bashyal::particle> particles(n);

        for (Foam::label i = 0; i < n; ++i)
        {
            Foam::pointField points = makeSpherePoints(radius, resolution);
            Foam::faceList faces = makeSphereLikeFaces(resolution);

            // Mass and inertia for a solid sphere of radius r, density=1
            Foam::scalar r = radius;
            Foam::scalar density = 1.0;
            Foam::scalar mass = density * (4.0/3.0) * Foam::constant::mathematical::pi * r*r*r;
            Foam::scalar I = (2.0/5.0) * mass * r*r;
            Foam::tensor moi = makeDiagonalTensor(I, I, I);

            Foam::point initialPosition(0, 0, 0);
            particles.set(i, makeParticleFromGeom(points, faces, mass, moi, initialPosition));
        }

        return particles;
    }
}


