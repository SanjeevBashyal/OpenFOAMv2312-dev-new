/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Test main function for head-on collision of two particles
    
\*---------------------------------------------------------------------------*/

#include "particle.H"
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
         << "Head-on Collision Test\n"
         << "========================================\n" << endl;

    // Simulation parameters
    scalar endTime = 0.1; // Short simulation
    scalar dt = 1e-5;
    scalar g = 0.0;  // No gravity for collision test
    
    // Contact parameters (same as stable packing)
    scalar kn = 2e5;
    scalar knUnloading = 2e5;
    scalar kt = 1.5e5;
    scalar gammaN = 50.0;
    scalar gammaT = 25.0;
    scalar muC = 0.5;
    scalar crustWidth = 1e-4;
    
    Info << "Simulation parameters:" << endl;
    Info << "  Time step: " << dt << " s" << endl;
    Info << "  End time: " << endTime << " s" << endl;
    
    // Create output directory
    Foam::mkDir("VTK_Collision");

    // Create GJK contact detector and force computer
    gjkContact contactDetector;
    contactForce forceComputer(kn, knUnloading, kt, gammaN, gammaT, muC);

    // Create Time Registry
    timeRegistry timeReg(dt, endTime, "VTK_Collision");
    timeReg.setWriteInterval(0.001); // Write frequently

    // Particle List
    List<particle> particles;

    scalar size = 0.02;
    scalar mass = 0.01;
    
    // Particle 1: Moving Right
    particle p1 = createCubeParticle(point(-0.05, 0, 0), size, mass);
    p1.setVelocity(vector(1.0, 0, 0)); // 1 m/s
    particles.append(p1);
    
    // Particle 2: Moving Left
    particle p2 = createCubeParticle(point(0.05, 0, 0), size, mass);
    p2.setVelocity(vector(-1.0, 0, 0)); // -1 m/s
    particles.append(p2);
    
    // Register particles
    forAll(particles, i)
    {
        timeReg.addParticle(particles[i], "particle_" + name(i));
    }
    
    Info << "Created " << particles.size() << " particles." << endl;

    // Time integration
    Info << "\nStarting time integration..." << endl;
    
    scalar t = 0.0;
    int step = 0;
    
    while (t < endTime + SMALL)
    {
        // 1. Clear forces
        forAll(particles, i)
        {
            particles[i].clearForceAndTorque();
            // No gravity
        }
        
        // 2. Detect and compute contacts
        for (int i = 0; i < particles.size(); ++i)
        {
            for (int j = i + 1; j < particles.size(); ++j)
            {
                ContactInfo contact = contactDetector.detectContact(particles[i], particles[j], crustWidth, crustWidth);
                if (contact.isContact)
                {
                    forceComputer.computeForce(particles[i], particles[j], contact, dt);
                    if (step % 100 == 0)
                    {
                         // Info << "Contact detected between " << i << " and " << j << endl;
                    }
                }
            }
        }
        
        // 3. Update particles
        forAll(particles, i)
        {
            particles[i].update(dt);
        }
        
        // 4. Write Output
        timeReg.setCurrentTime(t);
        timeReg.advanceTime();
        
        // Progress
        if (step % 100 == 0)
        {
             Info << "Time = " << t << " s (" << (t/endTime)*100.0 << "%)" 
                  << " P0 x: " << particles[0].position().x() 
                  << " P1 x: " << particles[1].position().x() << endl;
        }

        t += dt;
        step++;
    }
    
    timeReg.writeTimeSeries();

    Info << "\nSimulation completed!" << endl;
    return 0;
}
