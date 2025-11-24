/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Test Case 1: DEM Validation - The Hysteretic Bounce
    Validates GJK collision detection, Time Integration, and Walton-Braun model.
    
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
#include "OSspecific.H" 
#include <cmath>
#include <iostream>
#include <fstream>

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

// Helper function to write wall VTP
void writeWallVtp(const Foam::fileName& filename, const wall& w, scalar size = 2.0)
{
    Foam::OFstream vtpFile(filename);
    if (!vtpFile.good())
    {
        FatalErrorIn("writeWallVtp") << "Cannot open file " << filename << Foam::exit(Foam::FatalError);
    }

    // Create a quad for the floor (centered at 0,0,0, normal 0,0,1)
    // Vertices: (-1,-1,0), (1,-1,0), (1,1,0), (-1,1,0) scaled by size/2
    scalar s = size / 2.0;
    point p0(-s, -s, 0);
    point p1(s, -s, 0);
    point p2(s, s, 0);
    point p3(-s, s, 0);

    vtpFile << "<?xml version=\"1.0\"?>" << Foam::endl;
    vtpFile << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">" << Foam::endl;
    vtpFile << "  <PolyData>" << Foam::endl;
    vtpFile << "    <Piece NumberOfPoints=\"4\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"1\">" << Foam::endl;
    
    vtpFile << "      <Points>" << Foam::endl;
    vtpFile << "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << Foam::endl;
    vtpFile << "          " << p0.x() << " " << p0.y() << " " << p0.z() << Foam::endl;
    vtpFile << "          " << p1.x() << " " << p1.y() << " " << p1.z() << Foam::endl;
    vtpFile << "          " << p2.x() << " " << p2.y() << " " << p2.z() << Foam::endl;
    vtpFile << "          " << p3.x() << " " << p3.y() << " " << p3.z() << Foam::endl;
    vtpFile << "        </DataArray>" << Foam::endl;
    vtpFile << "      </Points>" << Foam::endl;

    vtpFile << "      <Polys>" << Foam::endl;
    vtpFile << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << Foam::endl;
    vtpFile << "          0 1 2 3" << Foam::endl;
    vtpFile << "        </DataArray>" << Foam::endl;
    vtpFile << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << Foam::endl;
    vtpFile << "          4" << Foam::endl;
    vtpFile << "        </DataArray>" << Foam::endl;
    vtpFile << "      </Polys>" << Foam::endl;

    vtpFile << "    </Piece>" << Foam::endl;
    vtpFile << "  </PolyData>" << Foam::endl;
    vtpFile << "</VTKFile>" << Foam::endl;
}

int main(int argc, char *argv[])
{
    Info << "\n========================================\n"
         << "Test Case 1: The Hysteretic Bounce\n"
         << "========================================\n" << endl;

    // 1. Setup Parameters
    scalar dt = 1e-5;           // Time step (s)
    scalar endTime = 1.5;       // Sufficient time for bounce (s)
    scalar g = 9.81;            // Gravity (m/s^2)
    
    // Particle Properties
    scalar L = 0.1;             // Side length (m)
    scalar m = 1.0;             // Mass (kg)
    scalar H0 = 1.0;            // Drop height (bottom face) (m)
    
    // Contact Parameters (Walton-Braun)
    scalar k1 = 1e5;            // Loading stiffness (N/m)
    scalar k2 = 4e5;            // Unloading stiffness (N/m)
    scalar mu = 0.0;            // Friction coefficient
    scalar gammaN = 0.0;        // No viscous damping for this validation
    scalar gammaT = 0.0;
    scalar kt = 0.0;            // No tangential stiffness needed for vertical drop
    scalar crustWidth = 1e-4;

    Info << "Parameters:" << endl;
    Info << "  Mass (m): " << m << " kg" << endl;
    Info << "  Size (L): " << L << " m" << endl;
    Info << "  Drop Height (H0): " << H0 << " m" << endl;
    Info << "  Loading Stiffness (k1): " << k1 << " N/m" << endl;
    Info << "  Unloading Stiffness (k2): " << k2 << " N/m" << endl;
    Info << "  Expected Coeff. Restitution (e): " << Foam::sqrt(k1/k2) << endl;
    Info << "  Expected Rebound Height (H1): " << H0 * (k1/k2) << " m" << endl;

    // 2. Create Objects
    Foam::mkDir("VTK_Bounce");
    
    // Wall (Floor at z=0)
    wall floor(point(0,0,0), vector(0,0,1), "floor");
    writeWallVtp("VTK_Bounce/wall.vtp", floor);
    
    // Particle (Center at H0 + L/2)
    point initialPos(0, 0, H0 + 0.5*L);
    particle p = createCubeParticle(initialPos, L, m);
    
    // Contact System
    gjkContact contactDetector;
    contactForce forceComputer(k1, k2, kt, gammaN, gammaT, mu);
    
    // Time Registry
    timeRegistry timeReg(dt, endTime, "VTK_Bounce");
    timeReg.setWriteInterval(0.01);
    timeReg.addParticle(p, "cube");
    timeReg.addStaticObject("floor", "VTK_Bounce/wall.vtp");

    // Data Logging
    std::ofstream dataFile("bounce_data.csv");
    dataFile << "Time,Z,Vz,ForceN,Overlap\n";

    // 3. Simulation Loop
    Info << "\nStarting simulation..." << endl;
    
    scalar t = 0.0;
    int step = 0;
    
    while (t < endTime + SMALL)
    {
        // A. Clear Forces
        p.clearForceAndTorque();
        
        // B. Apply Gravity
        vector Fg(0, 0, -m * g);
        p.applyForce(Fg, p.position());
        
        // C. Detect Wall Contact
        ContactInfo contact = contactDetector.detectContactParticleWall(p, floor, crustWidth);
        
        scalar Fn = 0.0;
        scalar overlap = 0.0;
        
        // Track active contacts for history cleanup
        List<string> activeKeys;

        if (contact.isContact)
        {
            // D. Compute Contact Force
            forceComputer.computeForceParticleWall(p, floor, contact, dt);
            
            // Add key to active list
            // Key format from contactForce.C: "W_" + w.name() + "_" + name(p.id())
            activeKeys.append("W_" + floor.name() + "_" + name(p.id()));
            
            // Extract normal force magnitude for logging
            vector F_contact = p.force() - Fg;
            Fn = F_contact.z();
            overlap = contact.depth;
            
            // if (step % 1000 == 0) Info << "Contact detected at t=" << t << " Depth=" << overlap << " Fn=" << Fn << endl;
        }
        
        // Clean up history for contacts that are no longer active
        forceComputer.cleanHistory(activeKeys);
        
        // E. Update Particle
        p.update(dt);
        
        // F. Output
        if (step % 100 == 0) // Log every 1ms
        {
            dataFile << t << "," << p.position().z() - 0.5*L << "," << p.velocity().z() << "," << Fn << "," << overlap << "\n";
        }
        
        timeReg.setCurrentTime(t);
        timeReg.advanceTime(); // Writes VTK if interval met
        
        if (step % 10000 == 0)
        {
             Info << "Time = " << t << " s, Z = " << p.position().z() - 0.5*L << " m" << endl;
        }

        t += dt;
        step++;
    }
    
    dataFile.close();
    timeReg.writeTimeSeries();
    
    Info << "\nSimulation completed. Data written to bounce_data.csv and VTK_Bounce/" << endl;

    return 0;
}
