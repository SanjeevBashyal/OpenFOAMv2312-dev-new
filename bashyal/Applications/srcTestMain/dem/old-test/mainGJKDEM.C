/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Test main function for GJK-based 3D convex polyhedral DEM
    
    Tests the implementation of:
    - GJK contact detection algorithm
    - Wall classes (plane, cylinder, sphere)
    - Contact force computation (Hysteretic)
    - Particle-particle and particle-wall interactions
    - VTP Output
    
    Based on the Grains3D paper implementation.

\*---------------------------------------------------------------------------*/

#include "particle.H"
#include "wall.H"
#include "gjkContact.H"
#include "contactForce.H"
#include "timeRegistry.H"
#include "OFstream.H"
#include "quaternion.H"
#include "point.H"
#include "vector.H"
#include "faceList.H"
#include "pointField.H"
#include "constants.H"
#include "List.H"
#include "OSspecific.H" // For mkDir
#include <cmath>
#include <iostream>

using namespace Bashyal;
using namespace Foam;

// Helper function to create a cube particle
particle createCubeParticle(
    const point& center,
    const scalar size,
    const scalar mass,
    const quaternion& orientation = quaternion::I
)
{
    pointField vertices(8);
    scalar halfSize = 0.5 * size;
    vertices[0] = point(-halfSize, -halfSize, -halfSize);
    vertices[1] = point(halfSize, -halfSize, -halfSize);
    vertices[2] = point(halfSize, halfSize, -halfSize);
    vertices[3] = point(-halfSize, halfSize, -halfSize);
    vertices[4] = point(-halfSize, -halfSize, halfSize);
    vertices[5] = point(halfSize, -halfSize, halfSize);
    vertices[6] = point(halfSize, halfSize, halfSize);
    vertices[7] = point(-halfSize, halfSize, halfSize);

    faceList faces(6);
    faces[0] = face({0, 3, 2, 1});  // Bottom
    faces[1] = face({4, 5, 6, 7});  // Top
    faces[2] = face({0, 1, 5, 4});  // Front
    faces[3] = face({2, 3, 7, 6});  // Back
    faces[4] = face({1, 2, 6, 5});  // Right
    faces[5] = face({3, 0, 4, 7});  // Left

    boundary cubeBoundary(vertices, faces);
    
    // Compute moment of inertia for a cube
    scalar I = (1.0/6.0) * mass * size * size;
    tensor moi = tensor::zero;
    moi.xx() = I;
    moi.yy() = I;
    moi.zz() = I;
    
    particle p(cubeBoundary, mass, moi, center);
    p.setOrientation(orientation);
    
    return p;
}

int main(int argc, char *argv[])
{
    Info << "\n========================================\n"
         << "GJK-based 3D Convex Polyhedral DEM Packing Test\n"
         << "========================================\n" << endl;

    // Simulation parameters
    scalar endTime = 1.0;
    scalar dt = 1e-4;
    scalar g = 9.81;  // Gravity
    
    // Contact parameters
    scalar kn = 4e5;           // Normal stiffness (loading)
    scalar knUnloading = 4e5;  // Normal stiffness (unloading) - set higher for hysteresis if desired
    scalar kt = 3e5;           // Tangential stiffness
    scalar gammaN = 200.0;     // Normal damping
    scalar gammaT = 100.0;     // Tangential damping
    scalar muC = 0.5;          // Coulomb friction
    scalar crustWidth = 1e-4;  // Crust width for GJK
    
    Info << "Simulation parameters:" << endl;
    Info << "  Time step: " << dt << " s" << endl;
    Info << "  End time: " << endTime << " s" << endl;
    Info << "  Contact stiffness kn: " << kn << " N/m" << endl;
    Info << "  Friction coefficient: " << muC << endl;
    
    // Create output directory
    Foam::mkDir("VTK");

    // Create GJK contact detector and force computer
    gjkContact contactDetector;
    contactForce forceComputer(kn, knUnloading, kt, gammaN, gammaT, muC);

    // Create particles (Grid)
    Info << "Creating particles..." << endl;
    List<particle> particles;
    
    scalar size = 0.02;
    scalar spacing = 0.04;
    scalar mass = 0.01;
    
    int nx = 3;
    int ny = 3;
    int nz = 3;
    
    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                point center(
                    (i - nx/2.0) * spacing + spacing/2.0,
                    (j - ny/2.0) * spacing + spacing/2.0,
                    0.1 + k * spacing
                );
                
                // Random orientation
                // For simplicity, just some rotation
                quaternion q(vector(1, 1, 1), 0.1 * (i+j+k));
                q.normalise();
                
                particles.append(createCubeParticle(center, size, mass, q));
            }
        }
    }
    
    Info << "Created " << particles.size() << " particles." << endl;

    // Create walls (Box)
    Info << "\nCreating walls..." << endl;
    List<wall> walls;
    
    // Bottom
    walls.append(wall(point(0, 0, 0), vector(0, 0, 1), "bottom"));
    // Top (open)
    
    // Sides
    walls.append(wall(point(-0.1, 0, 0), vector(1, 0, 0), "left"));
    walls.append(wall(point(0.1, 0, 0), vector(-1, 0, 0), "right"));
    walls.append(wall(point(0, -0.1, 0), vector(0, 1, 0), "front"));
    walls.append(wall(point(0, 0.1, 0), vector(0, -1, 0), "back"));
    
    forAll(walls, i)
    {
        walls[i].setContactParameters(kn, gammaN, gammaT, muC, crustWidth);
    }

    // Time integration
    Info << "\nStarting time integration..." << endl;
    
    scalar t = 0.0;
    int step = 0;
    int writeInterval = 100;  // Write every 100 steps (0.01s)
    
    while (t < endTime + SMALL)
    {
        // 1. Clear forces
        forAll(particles, i)
        {
            particles[i].clearForceAndTorque();
            // Apply gravity
            particles[i].applyForce(particles[i].mass() * vector(0, 0, -g), particles[i].position());
        }
        
        // 2. Detect and compute contacts
        
        // Particle-Particle (N^2 / 2)
        // Ideally use broadphase, but N is small (~27)
        for (int i = 0; i < particles.size(); ++i)
        {
            for (int j = i + 1; j < particles.size(); ++j)
            {
                // Simple bounding box check (optional optimization)
                if (mag(particles[i].position() - particles[j].position()) > (particles[i].getDC() + particles[j].getDC()))
                {
                    continue;
                }
                
                ContactInfo contact = contactDetector.detectContact(particles[i], particles[j], crustWidth, crustWidth);
                if (contact.isContact)
                {
                    forceComputer.computeForce(particles[i], particles[j], contact, dt);
                }
            }
        }
        
        // Particle-Wall
        for (int i = 0; i < particles.size(); ++i)
        {
            forAll(walls, w)
            {
                ContactInfo contact = contactDetector.detectContactParticleWall(particles[i], walls[w], crustWidth);
                if (contact.isContact)
                {
                    forceComputer.computeForceParticleWall(particles[i], walls[w], contact, dt);
                }
            }
        }
        
        // 3. Update particles
        forAll(particles, i)
        {
            particles[i].update(dt);
        }
        
        // 4. Write Output
        if (step % writeInterval == 0)
        {
            Info << "Time: " << t << " s" << endl;
            
            // Write VTP for each particle
            // Better: Write one VTP with all particles (MultiBlock or just append)
            // But particle::writeVtp writes a single file.
            // Let's write individual files for now: VTK/particle_ID_step.vtp
            // Or better: VTK/particles_step.vtp containing all geometry.
            // Since writeVtp is a member function, it writes that particle.
            // We can modify it or just loop.
            // Looping is easiest for now.
            
            forAll(particles, i)
            {
                fileName fName = "VTK/particle_" + name(i) + "_" + name(step) + ".vtp";
                particles[i].writeVtp(fName);
            }
        }
        
        t += dt;
        step++;
    }

    Info << "\nSimulation completed!" << endl;
    return 0;
}
