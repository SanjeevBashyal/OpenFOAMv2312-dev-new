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

// Helper function to create a wall particle (flattened box)
particle createWallParticle(
    const point& center,
    const vector& dimensions,
    const quaternion& orientation = quaternion::I
)
{
    pointField vertices(8);
    scalar dx = 0.5 * dimensions.x();
    scalar dy = 0.5 * dimensions.y();
    scalar dz = 0.5 * dimensions.z();

    vertices[0] = point(-dx, -dy, -dz);
    vertices[1] = point(dx, -dy, -dz);
    vertices[2] = point(dx, dy, -dz);
    vertices[3] = point(-dx, dy, -dz);
    vertices[4] = point(-dx, -dy, dz);
    vertices[5] = point(dx, -dy, dz);
    vertices[6] = point(dx, dy, dz);
    vertices[7] = point(-dx, dy, dz);

    faceList faces(6);
    faces[0] = face({0, 3, 2, 1});  // Bottom
    faces[1] = face({4, 5, 6, 7});  // Top
    faces[2] = face({0, 1, 5, 4});  // Front
    faces[3] = face({2, 3, 7, 6});  // Back
    faces[4] = face({1, 2, 6, 5});  // Right
    faces[5] = face({3, 0, 4, 7});  // Left

    boundary wallBoundary(vertices, faces);
    
    // Infinite mass/inertia for static body (handled by bodyOperation=0, but good to set high)
    scalar mass = 1e9;
    tensor moi = tensor::zero;
    moi.xx() = 1e9;
    moi.yy() = 1e9;
    moi.zz() = 1e9;
    
    particle p(wallBoundary, mass, moi, center);
    p.setOrientation(orientation);
    p.setBodyOperation(0); // 0: static body
    
    return p;
}

int main(int argc, char *argv[])
{
    Info << "\n========================================\n"
         << "GJK-based 3D Convex Polyhedral DEM Packing Test\n"
         << "========================================\n" << endl;

    // Simulation parameters
    scalar endTime = 1.0;
    scalar dt = 1e-4; // Reduced time step for stability with higher stiffness
    scalar g = 9.81;  // Gravity
    
    // Contact parameters
    scalar kn = 2e5;           // Normal stiffness (increased to reduce overlap)
    scalar knUnloading = 2e5;  // Normal stiffness (unloading)
    scalar kt = 1.5e5;         // Tangential stiffness
    scalar gammaN = 50.0;      // Normal damping (adjusted for higher stiffness)
    scalar gammaT = 25.0;      // Tangential damping
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

    // Create Time Registry
    timeRegistry timeReg(dt, endTime, "VTK");
    timeReg.setWriteInterval(0.01); // Write every 0.01s

    // Particle List
    List<particle> particles;

    // 1. Create Container (Static Particles)
    Info << "Creating container..." << endl;
    
    scalar containerWidth = 0.2;
    scalar containerDepth = 0.2;
    scalar containerHeight = 0.3; // Make it tall enough
    scalar wallThickness = 0.05; // Increased thickness to prevent penetration
    
    // Bottom
    particles.append(createWallParticle(
        point(0, 0, -wallThickness/2.0), 
        vector(containerWidth + 2*wallThickness, containerDepth + 2*wallThickness, wallThickness)
    ));
    
    // Left
    particles.append(createWallParticle(
        point(-containerWidth/2.0 - wallThickness/2.0, 0, containerHeight/2.0),
        vector(wallThickness, containerDepth, containerHeight)
    ));
    
    // Right
    particles.append(createWallParticle(
        point(containerWidth/2.0 + wallThickness/2.0, 0, containerHeight/2.0),
        vector(wallThickness, containerDepth, containerHeight)
    ));
    
    // Front
    particles.append(createWallParticle(
        point(0, -containerDepth/2.0 - wallThickness/2.0, containerHeight/2.0),
        vector(containerWidth + 2*wallThickness, wallThickness, containerHeight)
    ));
    
    // Back
    particles.append(createWallParticle(
        point(0, containerDepth/2.0 + wallThickness/2.0, containerHeight/2.0),
        vector(containerWidth + 2*wallThickness, wallThickness, containerHeight)
    ));

    Info << "Created " << particles.size() << " static wall particles." << endl;

    // 2. Create Dynamic Particles (Grid)
    Info << "Creating dynamic particles..." << endl;
    
    scalar size = 0.02;
    scalar spacing = 0.04;
    scalar mass = 0.01;
    
    int nx = 3;
    int ny = 3;
    int nz = 3;
    
    int startIdx = particles.size();

    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                point center(
                    (i - nx/2.0) * spacing + spacing/2.0,
                    (j - ny/2.0) * spacing + spacing/2.0,
                    0.1 + k * spacing // Start above the bottom
                );
                
                // Random orientation
                quaternion q(vector(1, 1, 1), 0.1 * (i+j+k));
                q.normalise();
                
                particle p = createCubeParticle(center, size, mass, q);
                particles.append(p);
            }
        }
    }
    
    Info << "Total particles: " << particles.size() << " (" << particles.size() - startIdx << " dynamic)" << endl;

    // Register particles with timeRegistry after list is stable
    forAll(particles, i)
    {
        timeReg.addParticle(particles[i], "particle_" + name(i));
    }
    

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
            
            // Apply gravity only to dynamic particles
            if (particles[i].getBodyOperation() != 0)
            {
                particles[i].applyForce(particles[i].mass() * vector(0, 0, -g), particles[i].position());
            }
        }
        
        // 2. Detect and compute contacts
        // Particle-Particle (N^2 / 2)
        for (int i = 0; i < particles.size(); ++i)
        {
            for (int j = i + 1; j < particles.size(); ++j)
            {
                // Skip if both are static
                if (particles[i].getBodyOperation() == 0 && particles[j].getBodyOperation() == 0)
                {
                    continue;
                }

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
        
        // 3. Update particles
        forAll(particles, i)
        {
            particles[i].update(dt);
        }
        
        // 4. Write Output (handled by timeRegistry)
        timeReg.setCurrentTime(t);
        timeReg.advanceTime();
        
        // Progress reporting
        if (step % 100 == 0)
        {
            Info << "Time = " << t << " s (" << (t/endTime)*100.0 << "%)" << endl;
        }

        t += dt;
        step++;
    }
    
    timeReg.writeTimeSeries();

    Info << "\nSimulation completed!" << endl;
    return 0;
}
