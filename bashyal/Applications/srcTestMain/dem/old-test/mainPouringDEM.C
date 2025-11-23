/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    DEM test that pours convex polyhedral particles (cubes) into a
    five-face cuboid container (floor + 4 side walls, open at the top).

    - Generates particles at a prescribed interval above the container
    - Solves particle–particle and particle–wall contacts with the GJK model
    - Uses the existing timeRegistry utilities to write per-particle VTP files
      so the motion can be inspected in ParaView
\*---------------------------------------------------------------------------*/

#include "particle.H"
#include "wall.H"
#include "gjkContact.H"
#include "contactForce.H"
#include "timeRegistry.H"

#include "faceList.H"
#include "pointField.H"
#include "quaternion.H"
#include "List.H"
#include "OFstream.H"
#include "constants.H"
#include "mathematicalConstants.H"

#include <memory>
#include <random>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>

using namespace Foam;
using namespace Bashyal;

namespace
{

//- Helper to create a cube-shaped particle of a given size, mass and centre
particle makeCubeParticle
(
    const point& centre,
    const scalar size,
    const scalar mass,
    const quaternion& orientation = quaternion::I
)
{
    pointField vertices(8);
    scalar halfSize = 0.5*size;

    vertices[0] = point(-halfSize, -halfSize, -halfSize);
    vertices[1] = point( halfSize, -halfSize, -halfSize);
    vertices[2] = point( halfSize,  halfSize, -halfSize);
    vertices[3] = point(-halfSize,  halfSize, -halfSize);
    vertices[4] = point(-halfSize, -halfSize,  halfSize);
    vertices[5] = point( halfSize, -halfSize,  halfSize);
    vertices[6] = point( halfSize,  halfSize,  halfSize);
    vertices[7] = point(-halfSize,  halfSize,  halfSize);

    faceList faces(6);
    faces[0] = face({0, 3, 2, 1}); // Bottom
    faces[1] = face({4, 5, 6, 7}); // Top
    faces[2] = face({0, 1, 5, 4}); // Front
    faces[3] = face({2, 3, 7, 6}); // Back
    faces[4] = face({1, 2, 6, 5}); // Right
    faces[5] = face({3, 0, 4, 7}); // Left

    boundary cubeBoundary(vertices, faces);

    scalar I = (1.0/6.0)*mass*size*size;
    tensor moi = tensor::zero;
    moi.xx() = I;
    moi.yy() = I;
    moi.zz() = I;

    particle p(cubeBoundary, mass, moi, centre);
    p.setOrientation(orientation);
    return p;
}

}

int main(int argc, char* argv[])
{
    Info<< nl
        << "==============================================" << nl
        << "  Pouring DEM test (open-top cuboid container)" << nl
        << "==============================================" << nl
        << endl;

    // --- Simulation controls
    const scalar dt          = 5e-4;
    const scalar endTime     = 2.0;
    const scalar gravityMag  = 9.81;

    // Particle properties
    const scalar particleSize  = 0.008;   // 8 mm cube
    const scalar particleMass  = 1.5e-3;  // 1.5 grams
    const scalar crustWidth    = 5e-5;

    // Pouring controls
    const label  maxParticles     = 40;
    const scalar spawnInterval    = 0.025;        // seconds between particles
    scalar       timeSinceLastSpawn = 0.0;
    label        spawnedCount     = 0;

    // Container geometry (centre at origin, open top)
    const scalar halfWidth  = 0.05;  // X extent (±)
    const scalar halfDepth  = 0.05;  // Y extent (±)
    const scalar floorZ     = 0.0;
    const scalar wallHeight = 0.12;
    const scalar pourHeight = wallHeight + 0.02;

    // Contact model coefficients (Hookean)
    const scalar kn     = 4e5;
    const scalar gammaN = 2.6e3;
    const scalar gammaT = 8e4;
    const scalar muC    = 0.5;

    // --- Utilities
    gjkContact contactDetector;
    contactForce forceComputer(kn, gammaN, gammaT, muC, 0.0);

    timeRegistry registry(dt, endTime, "VTK_pouringDEM");
    registry.setWriteInterval(0.01);
    registry.setWriteAtStart(true);
    registry.setWriteAtEnd(true);

    // Random generator for spawn locations and slight orientation jitter
    std::mt19937 rng(2025);
    std::uniform_real_distribution<scalar> randX(-halfWidth + particleSize, halfWidth - particleSize);
    std::uniform_real_distribution<scalar> randY(-halfDepth + particleSize, halfDepth - particleSize);
    const scalar angSpread = Foam::constant::mathematical::pi/12.0;
    std::uniform_real_distribution<scalar> randAngle(-angSpread, angSpread);

    // Create five plane walls (floor + four sides)
    auto makePlane = [&](const point& pt, const vector& normal, const word& name) -> wall
    {
        wall w(pt, normal, name);
        w.setContactParameters(kn, gammaN, gammaT, muC, crustWidth);
        return w;
    };

    List<wall> walls(5);
    walls[0] = makePlane(point(0, 0, floorZ), vector(0, 0, 1), "floor");
    walls[1] = makePlane(point(-halfWidth, 0, 0),  vector(1, 0, 0), "wall_xMin");
    walls[2] = makePlane(point( halfWidth, 0, 0),  vector(-1, 0, 0), "wall_xMax");
    walls[3] = makePlane(point(0, -halfDepth, 0),  vector(0, 1, 0), "wall_yMin");
    walls[4] = makePlane(point(0,  halfDepth, 0),  vector(0,-1, 0), "wall_yMax");

    // Particle storage
    std::vector<std::unique_ptr<particle>> particles;
    particles.reserve(maxParticles);

    auto spawnParticle = [&](label id)
    {
        point c(randX(rng), randY(rng), pourHeight);
        quaternion q = quaternion(vector(0,1,0), randAngle(rng));

        auto pPtr = std::make_unique<particle>(makeCubeParticle(c, particleSize, particleMass, q));
        pPtr->setVelocity(vector(0, 0, 0));
        pPtr->setAngularVelocity(vector::zero);
        pPtr->setIsActive(true);

        std::ostringstream name;
        name << "particle_" << std::setw(3) << std::setfill('0') << id;
        registry.addParticle(*pPtr, name.str());

        particles.emplace_back(std::move(pPtr));
        Info<< "Spawned particle " << name.str() << " at " << c << endl;
    };

    // Spawn a couple of seed particles so container is not empty
    for (label i = 0; i < 3 && spawnedCount < maxParticles; ++i)
    {
        spawnParticle(spawnedCount++);
    }

    const vector gravity(0, 0, -gravityMag);
    const label nSteps = label(std::ceil(endTime/dt));

    for (label step = 0; step <= nSteps; ++step)
    {
        const scalar time = step*dt;
        timeSinceLastSpawn += dt;

        if (spawnedCount < maxParticles && timeSinceLastSpawn >= spawnInterval)
        {
            spawnParticle(spawnedCount++);
            timeSinceLastSpawn = 0.0;
        }

        // Clear forces and apply gravity
        for (auto& pPtr : particles)
        {
            particle& p = *pPtr;
            p.clearForceAndTorque();
            p.applyForce(p.mass()*gravity, p.position());
        }

        // Particle-particle contacts
        const label nParticles = label(particles.size());
        for (label i = 0; i < nParticles; ++i)
        {
            particle& pA = *particles[i];
            for (label j = i+1; j < nParticles; ++j)
            {
                particle& pB = *particles[j];
                ContactInfo cAB = contactDetector.detectContact(pA, pB, crustWidth, crustWidth);
                if (cAB.isContact)
                {
                    forces F_ab = forceComputer.computeForce(pA, pB, cAB);
                    pA.updateContactForces(F_ab);

                    forces F_ba;
                    F_ba.F = -F_ab.F;
                    F_ba.T = -F_ab.T;
                    pB.updateContactForces(F_ba);
                }
            }
        }

        // Particle-wall contacts
        for (auto& pPtr : particles)
        {
            particle& p = *pPtr;
            for (const wall& w : walls)
            {
                ContactInfo cw = contactDetector.detectContactParticleWall(p, w, crustWidth);
                if (cw.isContact)
                {
                    forces Fw = forceComputer.computeForceParticleWall(p, w, cw);
                    p.updateContactForces(Fw);
                }
            }
        }

        // Integrate motion
        for (auto& pPtr : particles)
        {
            pPtr->updateMovement(dt);
        }

        // Write VTP outputs according to registry schedule
        registry.advanceTime();

        if (step % 200 == 0 || step == nSteps)
        {
            Info<< "Time " << time << " s, particles: " << particles.size() << endl;
        }
    }

    registry.writeTimeSeries();
    Info<< "Finished pouring DEM test. VTP files in " << "VTK_pouringDEM" << endl;

    return 0;
}

